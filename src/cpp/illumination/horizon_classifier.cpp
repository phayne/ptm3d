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
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
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

// ---- 16-bit elevation codes ------------------------------------------------------------
namespace {
constexpr double kQStep = M_PI / 65532.0;          // code 1 = -pi/2, code 65533 = +pi/2
constexpr double kQLo = -M_PI / 2;
}  // namespace

double qDecode(uint16_t c) {
    if (c == kQNegInf) return -std::numeric_limits<double>::infinity();
    if (c >= kQPosInf - 1) return std::numeric_limits<double>::infinity();   // 65534 is unused; treat as +inf
    return kQLo + (double)(c - 1) * kQStep;
}
uint16_t qUp(double e) {
    if (!(e > -1e8)) return kQNegInf;               // the -1e9 sentinel (and -inf)
    if (e > M_PI / 2) return kQPosInf;               // beyond the representable band, incl. +1e9
    double x = std::ceil((e - kQLo) / kQStep);
    long c = (long)x + 1; if (c < 1) c = 1; if (c > kQMaxFinite) c = kQMaxFinite;
    while (c < kQMaxFinite && qDecode((uint16_t)c) < e) ++c;     // guard the division's rounding
    while (c > 1 && qDecode((uint16_t)(c - 1)) >= e) --c;
    return qDecode((uint16_t)c) >= e ? (uint16_t)c : kQPosInf;
}
uint16_t qDown(double e) {
    if (!(e < 1e8)) return kQPosInf;                // the +1e9 sentinel (and +inf)
    if (e < -M_PI / 2) return kQNegInf;
    double x = std::floor((e - kQLo) / kQStep);
    long c = (long)x + 1; if (c < 1) c = 1; if (c > kQMaxFinite) c = kQMaxFinite;
    while (c > 1 && qDecode((uint16_t)c) > e) --c;
    while (c < kQMaxFinite && qDecode((uint16_t)(c + 1)) <= e) ++c;
    return qDecode((uint16_t)c) <= e ? (uint16_t)c : kQNegInf;
}

// ---- .hzc v2 ------------------------------------------------------------------------------
// little-endian:
//   header (128 B): char[8] "PTMHZC02"; u32 file version (2); u32 rule version; u64 scene; u64 predicate;
//     f64 dc, del, margin_deg, earth_radius_km; i32 nbis; u8 plane_guard; u8[3] 0; u64 nfaces_scene;
//     u64 note length; zero pad to 128;
//   note (padded with zeros to a multiple of 8);
//   data: per stored face, ascending: u16 nruns; per run: u16 start, u16 count, then count x
//     {u16 lo (up), hi (down), topB (up), botC (down)}; MULTI = topB +inf, botC -inf;
//     zero pad to a multiple of 8;
//   index: u32 face[nf] (zero pad to 8); u64 off[nf + 1] (bytes from the data start);
//   footer (48 B): u64 nfaces stored; u64 columns; u64 data offset; u64 index offset;
//     u64 checksum; char[8] "PTMHZCEN".
// checksum = FNV-1a 64 over little-endian u64 WORDS of every byte before the footer (the
// file length before the footer is a multiple of 8 by construction).
namespace {
constexpr char kMagic2[8] = {'P', 'T', 'M', 'H', 'Z', 'C', '0', '2'};
constexpr char kEnd[8] = {'P', 'T', 'M', 'H', 'Z', 'C', 'E', 'N'};
constexpr uint32_t kFileVersion = 2;
constexpr size_t kHeader = 128, kFooter = 48;
constexpr uint64_t kFnvOff = 1469598103934665603ULL, kFnvPrime = 1099511628211ULL;
inline uint64_t wordHash(uint64_t h, const unsigned char* p, size_t nwords) {
    for (size_t i = 0; i < nwords; ++i) { uint64_t w; std::memcpy(&w, p + 8 * i, 8); h ^= w; h *= kFnvPrime; }
    return h;
}
inline size_t pad8(size_t n) { return (8 - (n & 7)) & 7; }
}  // namespace

HzcWriter::~HzcWriter() { if (f_) std::fclose(static_cast<FILE*>(f_)); }

std::string HzcWriter::w(const void* p, size_t n) {
    if (n && std::fwrite(p, 1, n, static_cast<FILE*>(f_)) != n) return "write failed: " + path_;
    const unsigned char* c = static_cast<const unsigned char*>(p);
    size_t i = 0;
    if (!tail_.empty()) {
        while (i < n && tail_.size() < 8) tail_.push_back(c[i++]);
        if (tail_.size() == 8) { hash_ = wordHash(hash_, tail_.data(), 1); tail_.clear(); }
    }
    const size_t nw = (n - i) / 8;
    hash_ = wordHash(hash_, c + i, nw); i += 8 * nw;
    while (i < n) tail_.push_back(c[i++]);
    pos_ += n;
    return "";
}

std::string HzcWriter::open(const std::string& path, const Key& key, const Params& p, uint64_t nfaces_scene,
                            const std::string& note) {
    path_ = path; hash_ = kFnvOff; pos_ = 0; ncols_ = 0; tail_.clear(); faces_.clear(); offs_.clear();
    f_ = std::fopen(path.c_str(), "wb");
    if (!f_) return "cannot open " + path + " for writing";
    unsigned char h[kHeader] = {0}; size_t o = 0;
    auto put = [&](const void* v, size_t n) { std::memcpy(h + o, v, n); o += n; };
    const uint32_t fv = kFileVersion, rv = kRuleVersion; const int32_t nb = p.nbis;
    const uint8_t pg[4] = {(uint8_t)(p.plane_guard ? 1 : 0), 0, 0, 0}; const uint64_t nl = note.size();
    put(kMagic2, 8); put(&fv, 4); put(&rv, 4); put(&key.scene, 8); put(&key.predicate, 8);
    put(&p.dc_deg, 8); put(&p.del_deg, 8); put(&p.margin_deg, 8); put(&p.earth_radius_km, 8); put(&nb, 4); put(pg, 4);
    put(&nfaces_scene, 8); put(&nl, 8);
    std::string e = w(h, kHeader); if (!e.empty()) return e;
    e = w(note.data(), note.size()); if (!e.empty()) return e;
    static const unsigned char z[8] = {0};
    e = w(z, pad8(note.size())); if (!e.empty()) return e;
    data_off_ = pos_;
    offs_.push_back(0);
    return "";
}

std::string HzcWriter::append(uint32_t face, const FaceEnvelope& env) {
    if (!f_) return "append on a writer that is not open";
    if (!faces_.empty() && face <= faces_.back()) return "faces must be appended in strictly ascending order";
    buf_.clear(); buf_.push_back(0);
    uint16_t nruns = 0;
    for (size_t k = 0; k < env.cols.size();) {
        size_t m = k + 1;
        while (m < env.cols.size() && env.cols[m].col == env.cols[m - 1].col + 1 && m - k < 65535) ++m;
        buf_.push_back(env.cols[k].col); buf_.push_back((uint16_t)(m - k)); ++nruns;
        for (size_t j = k; j < m; ++j) {
            const Column& c = env.cols[j];
            buf_.push_back(qUp(c.lo)); buf_.push_back(qDown(c.hi));
            if (c.multi) { buf_.push_back(kQPosInf); buf_.push_back(kQNegInf); }
            else { buf_.push_back(qUp(c.topB)); buf_.push_back(qDown(c.botC)); }
        }
        k = m;
        if (nruns == 65535) return "too many column runs on one face";
    }
    buf_[0] = nruns;
    std::string e = w(buf_.data(), 2 * buf_.size()); if (!e.empty()) return e;
    ncols_ += env.cols.size();
    faces_.push_back(face); offs_.push_back(pos_ - data_off_);
    return "";
}

std::string HzcWriter::close() {
    if (!f_) return "close on a writer that is not open";
    static const unsigned char z[8] = {0};
    std::string e = w(z, pad8(pos_)); if (!e.empty()) return e;
    const uint64_t index_off = pos_, nf = faces_.size();
    e = w(faces_.data(), 4 * nf); if (!e.empty()) return e;
    e = w(z, pad8(pos_)); if (!e.empty()) return e;
    e = w(offs_.data(), 8 * offs_.size()); if (!e.empty()) return e;
    if (!tail_.empty()) return "internal: unaligned length before the footer";
    unsigned char ft[kFooter]; size_t o = 0;
    auto put = [&](const void* v, size_t n) { std::memcpy(ft + o, v, n); o += n; };
    put(&nf, 8); put(&ncols_, 8); put(&data_off_, 8); put(&index_off, 8); put(&hash_, 8); put(kEnd, 8);
    FILE* f = static_cast<FILE*>(f_); f_ = nullptr;
    const bool ok = std::fwrite(ft, 1, kFooter, f) == kFooter;
    const bool closed = std::fclose(f) == 0;
    if (!ok || !closed) return "write failed: " + path_;
    pos_ += kFooter;
    return "";
}

HzcReader::~HzcReader() { if (base_) munmap(const_cast<unsigned char*>(base_), size_); }

std::string HzcReader::open(const std::string& path, const Key* expect, bool verify_checksum) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) return "cannot open " + path;
    struct stat st; if (fstat(fd, &st) != 0) { ::close(fd); return "cannot stat " + path; }
    size_ = (size_t)st.st_size;
    if (size_ < kHeader + kFooter) { ::close(fd); return path + ": too short for a .hzc file"; }
    void* m = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd, 0); ::close(fd);
    if (m == MAP_FAILED) return "cannot map " + path;
    base_ = static_cast<const unsigned char*>(m);
    auto rd = [&](size_t off, void* v, size_t n) { std::memcpy(v, base_ + off, n); };
    if (std::memcmp(base_, kMagic2, 8) != 0) {
        if (std::memcmp(base_, "PTMHZC01", 8) == 0) return path + ": .hzc file version 1 (float32) is no longer read; rebuild it";
        return path + ": not a .hzc file (bad magic)";
    }
    uint32_t fv, rv; rd(8, &fv, 4); rd(12, &rv, 4);
    if (fv != kFileVersion) return path + ": unsupported .hzc file version " + std::to_string(fv);
    if (rv != kRuleVersion) return path + ": built with classifier rule v" + std::to_string(rv) + ", this build is v" + std::to_string(kRuleVersion);
    rd(16, &key_.scene, 8); rd(24, &key_.predicate, 8);
    if (expect && *expect != key_) {
        char b[256];
        std::snprintf(b, sizeof b, ": key mismatch (file scene %016llx predicate %016llx, run scene %016llx predicate %016llx)",
                      (unsigned long long)key_.scene, (unsigned long long)key_.predicate,
                      (unsigned long long)expect->scene, (unsigned long long)expect->predicate);
        return path + b + " -- this envelope was built for a different scene or ray predicate";
    }
    int32_t nb; uint8_t pg[4]; uint64_t nl;
    rd(32, &params_.dc_deg, 8); rd(40, &params_.del_deg, 8); rd(48, &params_.margin_deg, 8); rd(56, &params_.earth_radius_km, 8);
    rd(64, &nb, 4); rd(68, pg, 4); rd(72, &nfaces_scene_, 8); rd(80, &nl, 8);
    params_.nbis = nb; params_.plane_guard = pg[0] != 0;
    if (nl > (1u << 20) || kHeader + nl > size_ - kFooter) return path + ": bad note length";
    note_.assign(reinterpret_cast<const char*>(base_ + kHeader), nl);
    const size_t fo = size_ - kFooter;
    if (std::memcmp(base_ + fo + 40, kEnd, 8) != 0) return path + ": missing footer (truncated file?)";
    uint64_t ncol, cs; rd(fo, &nf_, 8); rd(fo + 8, &ncol, 8); rd(fo + 16, &data_off_, 8); rd(fo + 24, &index_off_, 8); rd(fo + 32, &cs, 8);
    const uint64_t idx_bytes = 4 * nf_ + pad8(4 * nf_) + 8 * (nf_ + 1);
    if (data_off_ != kHeader + nl + pad8(nl) || index_off_ < data_off_ || index_off_ + idx_bytes != fo || (fo & 7))
        return path + ": inconsistent section offsets";
    if (verify_checksum && wordHash(kFnvOff, base_, fo / 8) != cs)
        return path + ": checksum mismatch (corrupt or truncated file)";
    uint64_t last = 0;
    for (uint64_t i = 0; i < nf_; ++i) {
        const uint32_t fc = faceAt(i);
        if (i && fc <= last) return path + ": face index not ascending";
        last = fc;
    }
    return "";
}

uint32_t HzcReader::faceAt(uint64_t i) const { uint32_t f; std::memcpy(&f, base_ + index_off_ + 4 * i, 4); return f; }

bool HzcReader::get(uint32_t face, FaceEnvelope& out) const {
    out.cols.clear();
    uint64_t lo = 0, hi = nf_;
    while (lo < hi) { const uint64_t mid = (lo + hi) / 2; if (faceAt(mid) < face) lo = mid + 1; else hi = mid; }
    if (lo == nf_ || faceAt(lo) != face) return false;
    const size_t po = index_off_ + 4 * nf_ + pad8(4 * nf_);
    uint64_t a, b; std::memcpy(&a, base_ + po + 8 * lo, 8); std::memcpy(&b, base_ + po + 8 * (lo + 1), 8);
    const unsigned char* p = base_ + data_off_ + a; const unsigned char* e = base_ + data_off_ + b;
    auto u16 = [&]() { uint16_t v; std::memcpy(&v, p, 2); p += 2; return v; };
    if (p + 2 > e) return false;
    const uint16_t nruns = u16();
    for (uint16_t r = 0; r < nruns; ++r) {
        if (p + 4 > e) { out.cols.clear(); return false; }
        const uint16_t start = u16(), count = u16();
        if (p + 8 * (size_t)count > e) { out.cols.clear(); return false; }
        for (uint16_t j = 0; j < count; ++j) {
            Column c; c.col = (uint16_t)(start + j);
            const uint16_t ql = u16(), qh = u16(), qt = u16(), qb = u16();
            c.lo = qDecode(ql); c.hi = qDecode(qh);
            if (qt == kQPosInf && qb == kQNegInf) { c.multi = 1; c.topB = 1e9; c.botC = -1e9; }
            else { c.topB = qDecode(qt); c.botC = qDecode(qb); }
            out.cols.push_back(c);
        }
    }
    return true;
}

std::string writeHzc(const std::string& path, const HzcFile& f) {
    if (f.faces.size() != f.env.size()) return "faces/envelopes size mismatch";
    HzcWriter W; std::string e = W.open(path, f.key, f.params, f.nfaces_scene, f.build_note);
    for (size_t i = 0; e.empty() && i < f.faces.size(); ++i) e = W.append(f.faces[i], f.env[i]);
    if (e.empty()) e = W.close();
    return e;
}

std::string readHzc(const std::string& path, HzcFile& f, const Key* expect) {
    HzcReader R; std::string e = R.open(path, expect, true);
    if (!e.empty()) return e;
    f.key = R.key(); f.params = R.params(); f.nfaces_scene = R.numFacesScene(); f.build_note = R.note();
    f.faces.resize(R.numFaces()); f.env.assign(R.numFaces(), FaceEnvelope());
    for (uint64_t i = 0; i < R.numFaces(); ++i) {
        f.faces[i] = R.faceAt(i);
        if (!R.get(f.faces[i], f.env[i])) return path + ": corrupt face record " + std::to_string(f.faces[i]);
    }
    return "";
}

}  // namespace hzc
}  // namespace ptm
