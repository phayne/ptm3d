// horizon_classifier_scene.cpp -- see horizon_classifier_scene.h.
// Compiled without -ffast-math / FP contraction, as horizon_classifier.cpp (the curtain
// vertices and the recentering are part of the scene the hash certifies).
#include "horizon_classifier_scene.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace ptm {
namespace hzc {

namespace {
inline double norm(const V3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

double minCentroidR(const std::vector<V3>& V, const std::vector<F3>& F) {
    double m = 1e300;
    for (auto& f : F) {
        V3 c = {(V[f[0]][0] + V[f[1]][0] + V[f[2]][0]) / 3, (V[f[0]][1] + V[f[1]][1] + V[f[2]][1]) / 3,
                (V[f[0]][2] + V[f[1]][2] + V[f[2]][2]) / 3};
        m = std::min(m, norm(c));
    }
    return m;
}
double minVertR(const std::vector<V3>& V) { double m = 1e300; for (auto& v : V) m = std::min(m, norm(v)); return m; }

inline uint64_t edgeKey(int a, int b) {
    const uint32_t lo = a < b ? a : b, hi = a < b ? b : a;
    return ((uint64_t)lo << 32) | hi;
}

// Seam curtain: every boundary edge above radius rb drops a quad to rb (PTM-785).
void curtain(const std::vector<V3>& V, const std::vector<F3>& F, double rb, std::vector<V3>& cv, std::vector<F3>& cf) {
    std::unordered_map<uint64_t, uint32_t> cnt; cnt.reserve(F.size() * 3);
    for (auto& f : F) { ++cnt[edgeKey(f[0], f[1])]; ++cnt[edgeKey(f[1], f[2])]; ++cnt[edgeKey(f[2], f[0])]; }
    for (auto& f : F) for (int k = 0; k < 3; ++k) {
        const int a = f[k], b = f[(k + 1) % 3];
        if (cnt[edgeKey(a, b)] != 1) continue;
        const V3& A = V[a]; const V3& B = V[b]; const double ra = norm(A), rbb = norm(B);
        if (!(ra > rb && rbb > rb)) continue;
        const double sa = rb / ra, sb = rb / rbb; const int i0 = (int)cv.size();
        cv.push_back(A); cv.push_back(B); cv.push_back({B[0] * sb, B[1] * sb, B[2] * sb}); cv.push_back({A[0] * sa, A[1] * sa, A[2] * sa});
        cf.push_back({i0, i0 + 1, i0 + 2}); cf.push_back({i0, i0 + 2, i0 + 3});
    }
}

void addGeom(RTCDevice dev, RTCScene sc, const std::vector<V3>& V, const std::vector<F3>& F, const V3& off, Hash64& h) {
    RTCGeometry g = rtcNewGeometry(dev, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vb = (float*)rtcSetNewGeometryBuffer(g, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), V.size());
    for (size_t i = 0; i < V.size(); ++i) for (int k = 0; k < 3; ++k) vb[3 * i + k] = (float)(V[i][k] - off[k]);
    unsigned* ib = (unsigned*)rtcSetNewGeometryBuffer(g, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), F.size());
    for (size_t i = 0; i < F.size(); ++i) for (int k = 0; k < 3; ++k) ib[3 * i + k] = (unsigned)F[i][k];
    // The hash certifies what is TRACED: the float32 buffers, not the source doubles.
    const uint64_t nv = V.size(), nf = F.size();
    h.pod(nv); h.bytes(vb, 12 * nv); h.pod(nf); h.bytes(ib, 12 * nf);
    rtcCommitGeometry(g); rtcAttachGeometry(sc, g); rtcReleaseGeometry(g);
}

struct Ctx { RTCRayQueryContext base; unsigned src; };
void selfFilter(const RTCFilterFunctionNArguments* a) {
    const Ctx* q = reinterpret_cast<const Ctx*>(a->context);
    for (unsigned k = 0; k < a->N; ++k) {
        if (a->valid[k] != -1) continue;
        if (RTCHitN_geomID(a->hit, a->N, k) == 0 && RTCHitN_primID(a->hit, a->N, k) == q->src) a->valid[k] = 0;
    }
}
}  // namespace

ClassifierScene::ClassifierScene(const std::vector<V3>& Vc, const std::vector<F3>& Fc,
                                 const std::vector<V3>& Vr, const std::vector<F3>& Fr,
                                 std::vector<int8_t> sign, std::string rule)
    : Vc_(Vc), Fc_(Fc), sign_(std::move(sign)), orientation_rule_(std::move(rule)) {
    if (sign_.size() != Fc_.size()) throw std::runtime_error("ClassifierScene: orientation sign count != core faces");
    const bool ring = !Fr.empty();
    nring_ = Fr.size();
    facts_.R_occ = std::min(minCentroidR(Vc, Fc), ring ? minCentroidR(Vr, Fr) : 1e300);
    {
        std::unordered_map<uint64_t, int> c;
        for (auto& f : Fc) { ++c[edgeKey(f[0], f[1])]; ++c[edgeKey(f[1], f[2])]; ++c[edgeKey(f[2], f[0])]; }
        for (auto& kv : c) nbound_ += (kv.second == 1);
    }
    facts_.sphere_on = nbound_ > 0;
    std::vector<V3> Vk; std::vector<F3> Fk;
    if (ring) {
        const double rb = std::min({facts_.R_occ, minVertR(Vc), minVertR(Vr)}) - 1.0;
        curtain(Vc, Fc, rb, Vk, Fk); curtain(Vr, Fr, rb, Vk, Fk);
    }
    ncurtain_ = Fk.size();
    V3 bmin = {1e300, 1e300, 1e300}, bmax = {-1e300, -1e300, -1e300};
    for (auto* VV : {&Vc, &Vr}) for (auto& v : *VV) for (int k = 0; k < 3; ++k) {
        bmin[k] = std::min(bmin[k], v[k]); bmax[k] = std::max(bmax[k], v[k]);
    }
    off_ = {0.5 * (bmin[0] + bmax[0]), 0.5 * (bmin[1] + bmax[1]), 0.5 * (bmin[2] + bmax[2])};
    Hash64 h;
    h.str("ptm-hzc-scene");
    h.pod(off_); h.pod(facts_.R_occ); const uint8_t so = facts_.sphere_on ? 1 : 0; h.pod(so);
    dev_ = rtcNewDevice("");
    if (!dev_) throw std::runtime_error("ClassifierScene: rtcNewDevice failed");
    sc_ = rtcNewScene(dev_);
    rtcSetSceneFlags(sc_, RTC_SCENE_FLAG_ROBUST | RTC_SCENE_FLAG_COMPACT | RTC_SCENE_FLAG_FILTER_FUNCTION_IN_ARGUMENTS);
    rtcSetSceneBuildQuality(sc_, RTC_BUILD_QUALITY_HIGH);
    addGeom(dev_, sc_, Vc, Fc, off_, h);
    if (ring) { addGeom(dev_, sc_, Vr, Fr, off_, h); if (!Fk.empty()) addGeom(dev_, sc_, Vk, Fk, off_, h); }
    rtcCommitScene(sc_);
    scene_hash_ = h.h;
}

ClassifierScene::~ClassifierScene() {
    if (sc_) rtcReleaseScene(sc_);
    if (dev_) rtcReleaseDevice(dev_);
}

FaceFrame ClassifierScene::frame(size_t i) const {
    const F3& f = Fc_[i];
    return FaceFrame::make(Vc_[f[0]], Vc_[f[1]], Vc_[f[2]], sign_[i]);
}

int ClassifierScene::occluded(size_t core_face, const FaceFrame& f, const V3& d) const {
    const V3 org = {f.fc[0] - off_[0], f.fc[1] - off_[1], f.fc[2] - off_[2]};
    Ctx q; rtcInitRayQueryContext(&q.base); q.src = (unsigned)core_face;
    RTCRay r{};
    r.org_x = (float)org[0]; r.org_y = (float)org[1]; r.org_z = (float)org[2];
    r.dir_x = (float)d[0]; r.dir_y = (float)d[1]; r.dir_z = (float)d[2];
    r.tnear = 0.f; r.tfar = 1e30f; r.mask = 0xFFFFFFFFu;
    RTCOccludedArguments oa; rtcInitOccludedArguments(&oa);
    oa.flags = (RTCRayQueryFlags)(RTC_RAY_QUERY_FLAG_INCOHERENT | RTC_RAY_QUERY_FLAG_INVOKE_ARGUMENT_FILTER);
    oa.context = &q.base; oa.filter = selfFilter;
    rtcOccluded1(sc_, &r, &oa);
    return r.tfar < 0.f ? 1 : 0;
}

const char* ClassifierScene::predicateDescriptor() {
    return "embree4 rtcOccluded1 ROBUST|COMPACT|FILTER_IN_ARGS quality=HIGH; origin=face centroid (float32, "
           "recentered on core+ring bbox center); tnear=0 tfar=1e30; any-hit filter rejects geom 0 prim == source "
           "face only; R_occ sphere (min centroid radius core,ring) on open cores; seam curtain to min radius - 1 m";
}

Key ClassifierScene::key(const Params& p) const {
    return Key{scene_hash_, predicateHash(predicateDescriptor(), orientation_rule_, sign_, p)};
}

}  // namespace hzc
}  // namespace ptm
