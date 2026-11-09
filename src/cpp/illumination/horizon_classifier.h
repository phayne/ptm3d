// horizon_classifier.h -- per-face SOLAR horizon envelope used as a CLASSIFIER for the
// direct-solar term of the streamed history marcher (Fermi build spec 2026-11-09, M0).
//
// Library form of the 2026-11-08 measurement tool hzcls.cpp (spin-up study kit; not in the repo)
// (rule v2 = v1 + the face-plane guard). NOTHING in the simulator calls this yet: it is
// not linked into ptm_thermal_sim, so it cannot change a default.
//
// What it decides, per (face, epoch): whether the solar disk is certainly fully CLEAR,
// certainly fully BLOCKED, or must be integrated exactly (STRADDLE -> the engine chord).
// It is a screen, not an approximation: anything it cannot certify falls through to the
// exact path, and every fall-through is counted (ClassCensus). It never fills in a value.
//
// Envelope (per face; face-local frame up = unit(fc)): every column of DC degrees of
// azimuth that the disk (plus margin) visits over the BUILD epochs gets an elevation
// ladder over the union of the disk's elevation extents there [el - alpha - 2m,
// el + alpha + 2m] (top clamped at 89.9 deg), step DEL, then NBIS halvings of an adjacent
// blocked->clear bracket. Directions at/below the face plane are wildcards. A column with
// a clear sample below a blocked one is MULTI (never certifies).
//
// Classification of an epoch over every column the disk (plus margin) overlaps, plus one
// more on each side:
//   * FALL-THROUGH (-> STRADDLE, counted): a needed column was never built (unvisited),
//     or the disk's [bottom, top] leaves the column's built elevation band. This is what
//     makes reuse of one .hzc across windows fail-safe: epochs outside the build set can
//     only lose efficiency, never correctness.
//   * any MULTI column -> STRADDLE;
//   * CLEAR iff the disk bottom (el - alpha - m) >= the top of the blocked bracket in
//     every column; BLOCKED iff the disk top (el + alpha + m) <= the bottom of the clear
//     bracket in every column; else STRADDLE;
//   * face-plane guard (v2): a disk crossing the face plane (|c0| < sin(alpha + m))
//     always STRADDLEs (MEASURED necessary on n4990: every v1 dropout/false clear had it).
// Pre-culls (as the engine chord rule): whole disk below the face plane (FREE_PLANE),
// whole disk inside the R_occ sphere's occluding cone (FREE_CONE, open meshes),
// Earth-disk overlap (ECLIPSE -> the engine's own eclipse path).
//
// Persistence (.hzc): one file per site, keyed by TWO hashes, never by a face or
// triangle count (PTM-800: a BVH cache keyed on triangle count served a stale scene):
//   * scene hash     -- the float32 vertex and uint32 index buffers actually traced, in
//                       attach order, plus the recentering offset, R_occ and the sphere flag;
//   * predicate hash -- rule version, ray predicate descriptor, the applied orientation
//                       rule and every per-face orientation sign, and the build params.
// A reader with a different key gets an error and must fail loud.
// Stored values are float32 rounded CONSERVATIVELY (blocked tops and band bottoms up,
// clear bottoms and band tops down), so a round trip can only add STRADDLEs.
//
// Audit: a deterministic 0.1% sample of certified face-epochs (selected by a hash of the
// face index and the epoch's ET bits, so it is window-independent) is meant to be
// re-integrated with the exact chord by the caller and booked into AuditCensus.
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace ptm {
namespace hzc {

using V3 = std::array<double, 3>;

enum Class : uint8_t { FREE_PLANE = 0, FREE_CONE = 1, CLEAR = 2, BLOCKED = 3, STRADDLE = 4, ECLIPSE = 5,
                       UNSET = 255 };

constexpr uint32_t kRuleVersion = 2;            // v2 = v1 + face-plane guard
constexpr double   kElevationCapDeg = 89.9;     // top of any ladder
constexpr double   kAuditRate = 1e-3;           // 0.1% of certified face-epochs

struct Params {
    double dc_deg = 0.125;      // column width (azimuth)
    double del_deg = 0.05;      // ladder step (elevation)
    int    nbis = 8;            // bracket halvings
    double margin_deg = 0.02;   // disk margin m
    bool   plane_guard = true;  // rule v2
    double earth_radius_km = 6378.137;
    int ncol() const;
};

// Face geometry in the classifier's frame. fn is the oriented unit normal.
struct FaceFrame {
    V3 fc{}, fn{}, up{}, east{}, north{};
    double rf = 0.0;
    static FaceFrame make(const V3& A, const V3& B, const V3& C, int8_t orientation_sign);
};

// Per-epoch disk geometry for one face. k = pre-cull class or UNSET (needs the envelope).
struct EpochGeom {
    double az = 0, el = 0, alpha = 0, c0 = 0;
    V3 s{};
    uint8_t k = UNSET;
};

// Sun and Earth centers in the mesh frame, meters (the engine's SPICE positions x 1000).
struct Epoch {
    double et = 0.0;
    V3 sun_m{}, earth_m{};
};

// Scene facts the pre-culls need.
struct SceneFacts {
    bool sphere_on = false;     // open core: the R_occ sphere pre-check is on
    double R_occ = 0.0;         // meters
};

EpochGeom epochGeom(const FaceFrame& f, const Epoch& e, const SceneFacts& sf, const Params& p);

// Ray predicate for direction d (unit) from face f: 0 clear, 1 blocked. Only called for
// directions above the face plane and outside the R_occ sphere (the library tests those).
using OccludedFn = std::function<int(const V3& d)>;

struct Column {
    uint16_t col = 0;
    uint8_t multi = 0;
    double lo = 0, hi = 0;          // built elevation band (radians)
    double topB = -1e9, botC = 1e9; // top of the blocked bracket, bottom of the clear bracket
};

struct FaceEnvelope {
    std::vector<Column> cols;       // sorted by col
    const Column* find(int col) const;
};

// Builds the envelope from the epochs' geometry. Returns Embree rays cast (counted
// exactly as the standalone tool counts them) through *nrays.
FaceEnvelope buildEnvelope(const FaceFrame& f, const std::vector<EpochGeom>& g, const SceneFacts& sf,
                           const Params& p, const OccludedFn& occluded, uint64_t* nrays);

// Why a face-epoch got its class. Every STRADDLE has exactly one reason.
enum Reason : uint8_t {
    R_PRECULL = 0,          // FREE_PLANE / FREE_CONE / ECLIPSE
    R_CERTIFIED = 1,        // CLEAR / BLOCKED
    R_BRACKET = 2,          // disk overlaps a bracket (genuine straddle)
    R_MULTI = 3,            // a multi-valued column
    R_PLANE_GUARD = 4,      // disk crosses the face plane (v2)
    R_UNVISITED = 5,        // FALL-THROUGH: column never built
    R_OUT_OF_LADDER = 6,    // FALL-THROUGH: disk outside the column's built band
    R_COUNT = 7
};

struct ClassCensus {
    uint64_t by_class[6] = {0, 0, 0, 0, 0, 0};
    uint64_t by_reason[R_COUNT] = {0, 0, 0, 0, 0, 0, 0};
    void add(uint8_t cls, uint8_t reason) { ++by_class[cls]; ++by_reason[reason]; }
    void merge(const ClassCensus& o);
    uint64_t total() const;
    uint64_t fallthrough() const { return by_reason[R_UNVISITED] + by_reason[R_OUT_OF_LADDER]; }
    std::string summary() const;   // one line
};

uint8_t classify(const FaceEnvelope& env, const EpochGeom& q, const Params& p, uint8_t* reason);

// ---- audit ---------------------------------------------------------------------------
// Deterministic selection: same (face, et) -> same answer on every machine and window.
bool auditPick(uint32_t face, double et, double rate = kAuditRate);

struct AuditCensus {
    uint64_t n = 0;              // audited certified face-epochs
    uint64_t wrong_clear = 0;    // classified CLEAR, exact visible fraction < 1 - tol
    uint64_t wrong_blocked = 0;  // classified BLOCKED, exact visible fraction > tol
    double max_abs_frac = 0.0;   // max |classified - exact| as a fraction of the full disk
    // exact_frac = exact chord value / full-disk value (0 when the full disk is 0).
    void add(uint8_t cls, double exact_frac, double tol = 1e-6);
    void merge(const AuditCensus& o);
    bool clean() const { return wrong_clear == 0 && wrong_blocked == 0; }
    std::string summary() const;
};

// ---- hashing and the .hzc file ----------------------------------------------------------
struct Hash64 {   // FNV-1a, 64-bit; byte order is the host's (little-endian on both ISAs we run)
    uint64_t h = 1469598103934665603ULL;
    void bytes(const void* p, size_t n);
    template <class T> void pod(const T& v) { bytes(&v, sizeof(T)); }
    void str(const std::string& s) { uint64_t n = s.size(); pod(n); bytes(s.data(), s.size()); }
};

struct Key {
    uint64_t scene = 0, predicate = 0;
    bool operator==(const Key& o) const { return scene == o.scene && predicate == o.predicate; }
    bool operator!=(const Key& o) const { return !(*this == o); }
};

// predicate descriptor: the ray predicate's own description (from the scene builder).
uint64_t predicateHash(const std::string& predicate_descriptor, const std::string& orientation_rule,
                       const std::vector<int8_t>& orientation_sign, const Params& p);

struct HzcFile {
    Key key;
    Params params;
    uint64_t nfaces_scene = 0;            // core faces in the scene (informational only)
    std::vector<uint32_t> faces;          // stored face indices, ascending
    std::vector<FaceEnvelope> env;        // parallel to faces
    std::string build_note;               // free text (epoch span etc.)
};

// Values are rounded conservatively to float32 on write. Returns "" or an error message.
std::string writeHzc(const std::string& path, const HzcFile& f);
// Reads and verifies magic, version, payload checksum and (if expect != nullptr) the key.
// Returns "" or an error message; the caller must fail loud on any error.
std::string readHzc(const std::string& path, HzcFile& f, const Key* expect);

}  // namespace hzc
}  // namespace ptm
