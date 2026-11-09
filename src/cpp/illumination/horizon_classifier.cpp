// horizon_classifier.cpp -- see horizon_classifier.h.
//
// ⛔ Compiled WITHOUT -ffast-math and with FP contraction off (illumination/CMakeLists.txt).
// The classification compares disk extents against bracket ends that were found by the
// same arithmetic; a fused multiply-add or a reassociated sum moves those by an ulp and
// flips boundary cases. IEEE-exact arithmetic keeps the library bit-identical to the
// -O2 measurement tool it was ported from, and identical across x86 and arm64.
#include "horizon_classifier.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ptm {
namespace hzc {

namespace {
// The exact expressions of the measurement tool (hcommon.h). Do not "simplify".
inline double norm(const V3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
inline double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline V3 cross(const V3& a, const V3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline V3 unit(V3 a) { const double l = norm(a); return {a[0] / l, a[1] / l, a[2] / l}; }
constexpr double D2R = M_PI / 180;

inline V3 dirAE(const FaceFrame& f, double az, double el) {
    const double ce = std::cos(el), se = std::sin(el), ca = std::cos(az), sa = std::sin(az);
    return {ce * (ca * f.north[0] + sa * f.east[0]) + se * f.up[0],
            ce * (ca * f.north[1] + sa * f.east[1]) + se * f.up[1],
            ce * (ca * f.north[2] + sa * f.east[2]) + se * f.up[2]};
}

// Columns the disk (plus margin) overlaps, plus one more on each side; j may be out of
// [0, ncol) and is wrapped by the caller. Whole ring when the disk covers the zenith.
inline void colRange(const EpochGeom& q, double MRG, double DC, int ncol, int& j0, int& j1) {
    const double ce = std::cos(q.el);
    const double w = (ce <= std::sin(q.alpha + MRG)) ? M_PI : std::asin(std::sin(q.alpha + MRG) / ce);
    if (w >= M_PI * 0.999) { j0 = 0; j1 = ncol - 1; return; }
    j0 = (int)std::floor((q.az - w) / (DC * D2R)) - 1;
    j1 = (int)std::floor((q.az + w) / (DC * D2R)) + 2;
}
inline int wrap(int j, int n) { return ((j % n) + n) % n; }

// 0 clear, 1 blocked, 2 wildcard (at/below the face plane). Sphere test as blockedDir.
inline int status(const FaceFrame& f, const V3& d, const SceneFacts& sf, const OccludedFn& occ, uint64_t* nray) {
    if (!(dot(f.fn, d) > 0.0)) return 2;
    if (sf.sphere_on && dot(f.fc, d) < 0.0) {
        const double b = dot(f.fc, d), c = f.rf * f.rf - sf.R_occ * sf.R_occ, disc = b * b - c;
        if (disc >= 0) { const double q = std::sqrt(disc); if (-b - q > 0 || -b + q > 0) return 1; }
    }
    if (nray) ++*nray;
    return occ(d);
}

// Conservative float32: never certify more than the double envelope did.
inline float roundUp(double v) {
    float f = static_cast<float>(v);
    if (static_cast<double>(f) < v) f = std::nextafter(f, std::numeric_limits<float>::infinity());
    return f;
}
inline float roundDown(double v) {
    float f = static_cast<float>(v);
    if (static_cast<double>(f) > v) f = std::nextafter(f, -std::numeric_limits<float>::infinity());
    return f;
}
}  // namespace

int Params::ncol() const { return (int)std::lround(360.0 / dc_deg); }

FaceFrame FaceFrame::make(const V3& A, const V3& B, const V3& C, int8_t sign) {
    FaceFrame f;
    f.fc = {(A[0] + B[0] + C[0]) / 3, (A[1] + B[1] + C[1]) / 3, (A[2] + B[2] + C[2]) / 3};
    f.fn = unit(cross({B[0] - A[0], B[1] - A[1], B[2] - A[2]}, {C[0] - A[0], C[1] - A[1], C[2] - A[2]}));
    if (sign < 0) f.fn = {-f.fn[0], -f.fn[1], -f.fn[2]};
    f.rf = norm(f.fc);
    f.up = {f.fc[0] / f.rf, f.fc[1] / f.rf, f.fc[2] / f.rf};
    f.east = unit(cross({0, 0, 1}, f.up));
    f.north = cross(f.up, f.east);
    return f;
}

EpochGeom epochGeom(const FaceFrame& f, const Epoch& e, const SceneFacts& sf, const Params& p) {
    const double MRG = p.margin_deg * M_PI / 180;
    const double RE = p.earth_radius_km;
    const V3& fc = f.fc;
    const V3 sm = {e.sun_m[0] - fc[0], e.sun_m[1] - fc[1], e.sun_m[2] - fc[2]};
    const double sd = norm(sm);
    const V3 sv = {sm[0] / sd, sm[1] / sd, sm[2] / sd};
    double sa = 695700.0 / (sd / 1000); if (sa > 1) sa = 1;
    const double alpha = std::asin(sa);
    EpochGeom q; q.s = sv; q.alpha = alpha; q.c0 = dot(f.fn, sv);
    q.el = std::asin(std::max(-1.0, std::min(1.0, dot(sv, f.up))));
    q.az = std::atan2(dot(sv, f.east), dot(sv, f.north)); if (q.az < 0) q.az += 2 * M_PI;
    q.k = UNSET;
    if (q.c0 <= -std::sin(alpha)) { q.k = FREE_PLANE; return q; }
    if (sf.sphere_on && f.rf > sf.R_occ) {
        const double cs = std::max(-1.0, std::min(1.0, -dot(f.up, sv)));
        if (std::acos(cs) + alpha < std::asin(sf.R_occ / f.rf) - 1e-9) { q.k = FREE_CONE; return q; }
    }
    const V3 pc = {e.earth_m[0] - fc[0], e.earth_m[1] - fc[1], e.earth_m[2] - fc[2]};
    const double d = norm(pc);
    const double along = dot(pc, sv);
    if (d <= RE * 1000) { q.k = ECLIPSE; return q; }
    if (d < sd && along > 0) {
        const double sep = std::acos(std::max(-1.0, std::min(1.0, along / d)));
        if (sep < std::asin(RE * 1000 / d) + alpha + MRG) { q.k = ECLIPSE; return q; }
    }
    return q;
}

const Column* FaceEnvelope::find(int col) const {
    auto it = std::lower_bound(cols.begin(), cols.end(), col,
                               [](const Column& c, int v) { return (int)c.col < v; });
    return (it != cols.end() && (int)it->col == col) ? &*it : nullptr;
}

FaceEnvelope buildEnvelope(const FaceFrame& f, const std::vector<EpochGeom>& g, const SceneFacts& sf,
                           const Params& p, const OccludedFn& occluded, uint64_t* nrays) {
    const double DC = p.dc_deg, DEL = p.del_deg, MRG = p.margin_deg * M_PI / 180;
    const int NCOL = p.ncol();
    if (NCOL < 4 || NCOL > 65536)   // Column::col is uint16_t
        throw std::invalid_argument("horizon classifier: column width " + std::to_string(DC) +
                                    " deg gives " + std::to_string(NCOL) + " columns (need 4..65536)");
    std::vector<double> blo(NCOL, 1e9), bhi(NCOL, -1e9);
    for (const EpochGeom& q : g) {
        if (q.k != UNSET) continue;
        int j0, j1; colRange(q, MRG, DC, NCOL, j0, j1);
        for (int j = j0; j <= j1; ++j) {
            const int jj = wrap(j, NCOL);
            blo[jj] = std::min(blo[jj], q.el - q.alpha - 2 * MRG);
            bhi[jj] = std::max(bhi[jj], q.el + q.alpha + 2 * MRG);
        }
    }
    FaceEnvelope env;
    std::vector<double> ev; std::vector<int> st;
    uint64_t nray = 0;
    for (int j = 0; j < NCOL; ++j) {
        if (blo[j] > bhi[j]) continue;
        const double az = j * DC * D2R;
        const double lo = blo[j], hi = std::min(bhi[j], kElevationCapDeg * D2R);
        const int n = std::max(2, (int)std::ceil((hi - lo) / (DEL * D2R)) + 1);
        ev.resize(n); st.resize(n);
        for (int k = 0; k < n; ++k) { ev[k] = lo + (hi - lo) * k / (n - 1); st[k] = status(f, dirAE(f, az, ev[k]), sf, occluded, &nray); }
        Column c; c.col = (uint16_t)j; c.lo = lo; c.hi = hi;
        int lastB = -1, firstC = -1;
        for (int k = 0; k < n; ++k) { if (st[k] == 1) lastB = k; if (st[k] == 0 && firstC < 0) firstC = k; }
        if (lastB >= 0 && firstC >= 0 && firstC < lastB) { c.multi = 1; env.cols.push_back(c); continue; }
        double tb = -1e9, bc = 1e9;
        if (lastB >= 0) tb = (lastB + 1 < n) ? ev[lastB + 1] : 1e9;   // above the last blocked sample: not yet known clear
        if (firstC >= 0) bc = (firstC > 0) ? ev[firstC - 1] : -1e9;
        if (lastB >= 0 && firstC == lastB + 1) {   // adjacent blocked->clear: bisect the bracket
            double a = ev[lastB], b = ev[firstC];
            for (int it = 0; it < p.nbis; ++it) {
                const double m = 0.5 * (a + b);
                const int r = status(f, dirAE(f, az, m), sf, occluded, &nray);
                if (r == 1) a = m; else if (r == 0) b = m; else break;
            }
            tb = b; bc = a;
        }
        c.topB = tb; c.botC = bc;
        env.cols.push_back(c);
    }
    if (nrays) *nrays += nray;
    return env;
}

uint8_t classify(const FaceEnvelope& env, const EpochGeom& q, const Params& p, uint8_t* reason) {
    uint8_t rs = R_PRECULL;
    uint8_t k = q.k;
    if (k == UNSET) {
        const double MRG = p.margin_deg * M_PI / 180;
        const int NCOL = p.ncol();
        int j0, j1; colRange(q, MRG, p.dc_deg, NCOL, j0, j1);
        const double bot = q.el - q.alpha - MRG, top = q.el + q.alpha + MRG;
        // Fall-through first: a column the build never visited, or a disk outside the band
        // the build sampled, cannot be certified by this envelope.
        for (int j = j0; j <= j1 && rs == R_PRECULL; ++j) {
            const Column* c = env.find(wrap(j, NCOL));
            if (!c) rs = R_UNVISITED;
            else if (bot < c->lo || top > c->hi) rs = R_OUT_OF_LADDER;
        }
        if (rs != R_PRECULL) {
            k = STRADDLE;
        } else {
            bool clr = true, blk = true, multi = false;
            for (int j = j0; j <= j1 && (clr || blk); ++j) {
                const Column* c = env.find(wrap(j, NCOL));
                if (c->multi) { clr = blk = false; multi = true; break; }
                if (!(bot >= c->topB)) clr = false;
                if (!(top <= c->botC)) blk = false;
            }
            k = clr ? CLEAR : (blk ? BLOCKED : STRADDLE);
            rs = (k == STRADDLE) ? (multi ? R_MULTI : R_BRACKET) : R_CERTIFIED;
            if (p.plane_guard && std::fabs(q.c0) < std::sin(q.alpha + MRG)) {
                if (k != STRADDLE || rs == R_BRACKET) rs = R_PLANE_GUARD;
                k = STRADDLE;
            }
        }
    }
    if (reason) *reason = rs;
    return k;
}

void ClassCensus::merge(const ClassCensus& o) {
    for (int i = 0; i < 6; ++i) by_class[i] += o.by_class[i];
    for (int i = 0; i < R_COUNT; ++i) by_reason[i] += o.by_reason[i];
}
uint64_t ClassCensus::total() const { uint64_t t = 0; for (auto v : by_class) t += v; return t; }
std::string ClassCensus::summary() const {
    const double t = std::max<double>(1.0, (double)total());
    char b[640];
    std::snprintf(b, sizeof b,
                  "face-epochs %llu: free-plane %.4f%% free-cone %.4f%% CLEAR %.4f%% BLOCKED %.4f%% STRADDLE %.4f%% "
                  "ECLIPSE %.4f%% | straddle reasons: bracket %llu multi %llu plane-guard %llu FALL-THROUGH "
                  "unvisited %llu out-of-ladder %llu",
                  (unsigned long long)total(), 100 * by_class[0] / t, 100 * by_class[1] / t, 100 * by_class[2] / t,
                  100 * by_class[3] / t, 100 * by_class[4] / t, 100 * by_class[5] / t,
                  (unsigned long long)by_reason[R_BRACKET], (unsigned long long)by_reason[R_MULTI],
                  (unsigned long long)by_reason[R_PLANE_GUARD], (unsigned long long)by_reason[R_UNVISITED],
                  (unsigned long long)by_reason[R_OUT_OF_LADDER]);
    return b;
}

bool auditPick(uint32_t face, double et, double rate) {
    uint64_t eb; std::memcpy(&eb, &et, 8);
    uint64_t h = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)face * 0xBF58476D1CE4E5B9ULL) ^ (eb * 0x94D049BB133111EBULL);
    h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 27; h *= 0x94D049BB133111EBULL; h ^= h >> 31;
    return (double)(h >> 11) * (1.0 / 9007199254740992.0) < rate;
}

void AuditCensus::add(uint8_t cls, double exact_frac, double tol) {
    ++n;
    const double classified = (cls == CLEAR) ? 1.0 : 0.0;
    if (cls == CLEAR && exact_frac < 1.0 - tol) ++wrong_clear;
    if (cls == BLOCKED && exact_frac > tol) ++wrong_blocked;
    max_abs_frac = std::max(max_abs_frac, std::fabs(classified - exact_frac));
}
void AuditCensus::merge(const AuditCensus& o) {
    n += o.n; wrong_clear += o.wrong_clear; wrong_blocked += o.wrong_blocked;
    max_abs_frac = std::max(max_abs_frac, o.max_abs_frac);
}
std::string AuditCensus::summary() const {
    char b[256];
    std::snprintf(b, sizeof b, "audit %llu certified face-epochs: wrong CLEAR %llu, wrong BLOCKED %llu, max |frac err| %.3g%s",
                  (unsigned long long)n, (unsigned long long)wrong_clear, (unsigned long long)wrong_blocked,
                  max_abs_frac, clean() ? "" : "  <-- CLASSIFIER DISAGREES WITH THE CHORD");
    return b;
}

void Hash64::bytes(const void* p, size_t n) {
    const unsigned char* c = static_cast<const unsigned char*>(p);
    for (size_t i = 0; i < n; ++i) { h ^= c[i]; h *= 1099511628211ULL; }
}

uint64_t predicateHash(const std::string& desc, const std::string& orientation_rule,
                       const std::vector<int8_t>& sign, const Params& p) {
    Hash64 h;
    h.str("ptm-hzc-predicate");
    h.pod(kRuleVersion);
    h.str(desc);
    h.str(orientation_rule);
    uint64_t n = sign.size(); h.pod(n); h.bytes(sign.data(), sign.size());
    h.pod(p.dc_deg); h.pod(p.del_deg); h.pod(p.nbis); h.pod(p.margin_deg);
    const uint8_t pg = p.plane_guard ? 1 : 0; h.pod(pg); h.pod(p.earth_radius_km);
    h.pod(kElevationCapDeg);
    return h.h;
}

// ---- .hzc ---------------------------------------------------------------------------
// little-endian:
//   char[8] "PTMHZC01"; u32 version (1); u32 rule version; u64 scene hash; u64 predicate hash;
//   f64 dc, del, margin_deg, earth_radius_km; i32 nbis; u8 plane_guard; u8[3] 0;
//   u64 nfaces_scene; u64 nfaces stored; u64 ncols total; u64 note length; char note[];
//   u32 face[nf]; u64 colptr[nf+1];
//   per column: u16 col; u8 multi; u8 0; f32 lo (rounded up), hi (down), topB (up), botC (down);
//   u64 FNV-1a of every byte above (the payload checksum).
namespace {
constexpr char kMagic[8] = {'P', 'T', 'M', 'H', 'Z', 'C', '0', '1'};
constexpr uint32_t kFileVersion = 1;
#pragma pack(push, 1)
struct ColRec { uint16_t col; uint8_t multi; uint8_t pad; float lo, hi, topB, botC; };
#pragma pack(pop)
static_assert(sizeof(ColRec) == 20, "ColRec layout");

struct HashingWriter {
    std::ofstream o; Hash64 h;
    void w(const void* p, size_t n) { o.write(static_cast<const char*>(p), (std::streamsize)n); h.bytes(p, n); }
    template <class T> void pod(const T& v) { w(&v, sizeof(T)); }
};
struct HashingReader {
    std::ifstream i; Hash64 h;
    bool r(void* p, size_t n) { i.read(static_cast<char*>(p), (std::streamsize)n); if (!i) return false; h.bytes(p, n); return true; }
    template <class T> bool pod(T& v) { return r(&v, sizeof(T)); }
};
}  // namespace

std::string writeHzc(const std::string& path, const HzcFile& f) {
    if (f.faces.size() != f.env.size()) return "faces/envelopes size mismatch";
    HashingWriter W; W.o.open(path, std::ios::binary | std::ios::trunc);
    if (!W.o) return "cannot open " + path + " for writing";
    W.w(kMagic, 8); W.pod(kFileVersion); W.pod(kRuleVersion);
    W.pod(f.key.scene); W.pod(f.key.predicate);
    W.pod(f.params.dc_deg); W.pod(f.params.del_deg); W.pod(f.params.margin_deg); W.pod(f.params.earth_radius_km);
    const int32_t nb = f.params.nbis; W.pod(nb);
    const uint8_t pg[4] = {(uint8_t)(f.params.plane_guard ? 1 : 0), 0, 0, 0}; W.w(pg, 4);
    const uint64_t nfs = f.nfaces_scene, nf = f.faces.size();
    uint64_t nc = 0; for (const auto& e : f.env) nc += e.cols.size();
    W.pod(nfs); W.pod(nf); W.pod(nc);
    const uint64_t nl = f.build_note.size(); W.pod(nl); W.w(f.build_note.data(), nl);
    W.w(f.faces.data(), 4 * nf);
    std::vector<uint64_t> ptr(nf + 1, 0);
    for (size_t i = 0; i < nf; ++i) ptr[i + 1] = ptr[i] + f.env[i].cols.size();
    W.w(ptr.data(), 8 * (nf + 1));
    std::vector<ColRec> buf;
    for (const auto& e : f.env) {
        buf.resize(e.cols.size());
        for (size_t k = 0; k < e.cols.size(); ++k) {
            const Column& c = e.cols[k];
            buf[k] = {c.col, c.multi, 0, roundUp(c.lo), roundDown(c.hi), roundUp(c.topB), roundDown(c.botC)};
        }
        W.w(buf.data(), sizeof(ColRec) * buf.size());
    }
    const uint64_t cs = W.h.h; W.o.write(reinterpret_cast<const char*>(&cs), 8);
    W.o.close();
    if (!W.o) return "write failed: " + path;
    return "";
}

std::string readHzc(const std::string& path, HzcFile& f, const Key* expect) {
    HashingReader R; R.i.open(path, std::ios::binary);
    if (!R.i) return "cannot open " + path;
    char mg[8]; uint32_t fv = 0, rv = 0;
    if (!R.r(mg, 8) || std::memcmp(mg, kMagic, 8) != 0) return path + ": not a .hzc file (bad magic)";
    if (!R.pod(fv) || fv != kFileVersion) return path + ": unsupported .hzc file version " + std::to_string(fv);
    if (!R.pod(rv) || rv != kRuleVersion) return path + ": built with classifier rule v" + std::to_string(rv) +
                                                 ", this build is v" + std::to_string(kRuleVersion);
    if (!R.pod(f.key.scene) || !R.pod(f.key.predicate)) return path + ": truncated header";
    if (expect && *expect != f.key) {
        char b[256];
        std::snprintf(b, sizeof b, ": key mismatch (file scene %016llx predicate %016llx, run scene %016llx predicate %016llx)",
                      (unsigned long long)f.key.scene, (unsigned long long)f.key.predicate,
                      (unsigned long long)expect->scene, (unsigned long long)expect->predicate);
        return path + b + " -- this envelope was built for a different scene or ray predicate";
    }
    int32_t nb = 0; uint8_t pg[4];
    if (!R.pod(f.params.dc_deg) || !R.pod(f.params.del_deg) || !R.pod(f.params.margin_deg) ||
        !R.pod(f.params.earth_radius_km) || !R.pod(nb) || !R.r(pg, 4)) return path + ": truncated header";
    f.params.nbis = nb; f.params.plane_guard = pg[0] != 0;
    uint64_t nfs = 0, nf = 0, nc = 0, nl = 0;
    if (!R.pod(nfs) || !R.pod(nf) || !R.pod(nc) || !R.pod(nl) || nl > (1u << 20)) return path + ": truncated header";
    f.nfaces_scene = nfs;
    f.build_note.assign(nl, '\0'); if (nl && !R.r(&f.build_note[0], nl)) return path + ": truncated note";
    f.faces.resize(nf); if (!R.r(f.faces.data(), 4 * nf)) return path + ": truncated face list";
    std::vector<uint64_t> ptr(nf + 1);
    if (!R.r(ptr.data(), 8 * (nf + 1)) || ptr[0] != 0 || ptr[nf] != nc) return path + ": bad column pointers";
    f.env.assign(nf, FaceEnvelope());
    std::vector<ColRec> buf;
    for (size_t i = 0; i < nf; ++i) {
        if (ptr[i + 1] < ptr[i]) return path + ": bad column pointers";
        buf.resize(ptr[i + 1] - ptr[i]);
        if (!R.r(buf.data(), sizeof(ColRec) * buf.size())) return path + ": truncated columns";
        auto& cols = f.env[i].cols; cols.resize(buf.size());
        for (size_t k = 0; k < buf.size(); ++k) {
            const ColRec& r = buf[k];
            cols[k].col = r.col; cols[k].multi = r.multi;
            cols[k].lo = r.lo; cols[k].hi = r.hi; cols[k].topB = r.topB; cols[k].botC = r.botC;
            if (k && cols[k].col <= cols[k - 1].col) return path + ": columns not ascending";
            if ((int)cols[k].col >= f.params.ncol()) return path + ": column index beyond the column count";
        }
    }
    const uint64_t want = R.h.h; uint64_t cs = 0;
    R.i.read(reinterpret_cast<char*>(&cs), 8);
    if (!R.i || cs != want) return path + ": payload checksum mismatch (corrupt or truncated file)";
    return "";
}

}  // namespace hzc
}  // namespace ptm
