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
// Stored values (file v2) are 16-bit elevation codes rounded CONSERVATIVELY -- toward STRADDLE:
// blocked-bracket tops and band bottoms up, clear-bracket bottoms and band tops down -- so a
// round trip can only send more face-epochs to the chord, never certify one the double
// envelope did not (test_horizon_classifier A6, ptm_hzcls HZ_ROUNDTRIP). Code 0 = -inf,
// 65535 = +inf, finite 1..65533 in steps of pi/65532 rad (0.0027 deg) from -pi/2. A MULTI
// column is stored as (topB, botC) = (+inf, -inf), a pair no other column can produce
// (it would need a blocked top sample above a clear bottom sample). Columns are stored as
// runs (start, count), 8 bytes per column + 4 per run. The writer streams face by face
// (bounded RAM); the reader maps the file and decodes one face on demand.
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
    // The disk being classified (M2, 2026-11-09): the Sun by default; the Earth cache uses the Earth's
    // radius and no eclipse pre-cull (Epoch::sun_m then carries the EARTH's center). Defaults leave every
    // solar key and file unchanged.
    double disk_radius_km = 695700.0;
    bool eclipse_cull = true;
    bool isSunDisk() const { return disk_radius_km == 695700.0 && eclipse_cull; }
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

// Disk (Sun, or the Earth for the Earth cache) and Earth centers in the mesh frame, meters (the engine's SPICE
// positions x 1000). The eclipse pre-cull (Sun disk only) reads earth_m.
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

// ---- 16-bit elevation codes ------------------------------------------------------------
constexpr uint16_t kQNegInf = 0, kQPosInf = 65535, kQMaxFinite = 65533;
double qDecode(uint16_t c);           // -inf / +inf for the sentinels
uint16_t qUp(double e);               // smallest code with qDecode(code) >= e
uint16_t qDown(double e);             // largest code with qDecode(code) <= e

// ---- the .hzc file (v2) --------------------------------------------------------------------
// Streams envelopes to disk in ascending face order; RAM = the face index (12 B/face).
class HzcWriter {
public:
    // Returns "" or an error message (the caller must fail loud).
    std::string open(const std::string& path, const Key& key, const Params& p, uint64_t nfaces_scene,
                     const std::string& note);
    std::string append(uint32_t face, const FaceEnvelope& env);   // faces strictly ascending
    std::string close();
    ~HzcWriter();
    uint64_t bytesWritten() const { return pos_; }
    uint64_t columns() const { return ncols_; }
private:
    std::string w(const void* p, size_t n);
    void* f_ = nullptr;        // FILE*
    std::string path_;
    uint64_t pos_ = 0, data_off_ = 0, ncols_ = 0;
    uint64_t hash_ = 0;        // word hash state, see the .cpp
    std::vector<unsigned char> tail_;   // < 8 bytes not yet hashed
    std::vector<uint32_t> faces_;
    std::vector<uint64_t> offs_;
    std::vector<uint16_t> buf_;
};

// Maps a .hzc read-only. get() decodes one face's columns into doubles.
class HzcReader {
public:
    // Verifies magic, version, rule version, footer, offsets, the key (if expect), and -- when
    // verify_checksum -- the whole-file checksum (one pass over the file). "" or an error.
    std::string open(const std::string& path, const Key* expect, bool verify_checksum = true);
    ~HzcReader();
    const Key& key() const { return key_; }
    const Params& params() const { return params_; }
    uint64_t numFaces() const { return nf_; }
    uint64_t numFacesScene() const { return nfaces_scene_; }
    uint32_t faceAt(uint64_t i) const;
    const std::string& note() const { return note_; }
    // false if the face is not stored (the caller treats every epoch of it as fall-through).
    bool get(uint32_t face, FaceEnvelope& out) const;
private:
    const unsigned char* base_ = nullptr; size_t size_ = 0;
    Key key_; Params params_; std::string note_;
    uint64_t nf_ = 0, nfaces_scene_ = 0, data_off_ = 0, index_off_ = 0;
};

// Whole-file convenience for small envelopes (tests, tools).
struct HzcFile {
    Key key;
    Params params;
    uint64_t nfaces_scene = 0;            // core faces in the scene (informational only)
    std::vector<uint32_t> faces;          // stored face indices, ascending
    std::vector<FaceEnvelope> env;        // parallel to faces
    std::string build_note;               // free text (epoch span etc.)
};
std::string writeHzc(const std::string& path, const HzcFile& f);
std::string readHzc(const std::string& path, HzcFile& f, const Key* expect);

}  // namespace hzc
}  // namespace ptm
