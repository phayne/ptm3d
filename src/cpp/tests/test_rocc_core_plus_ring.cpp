// PTM-780 regression test: the open-mesh body sphere R_occ must be the minimum
// face-centroid radius over the core AND the occluder ring, not the core alone.
//
// Geometry (meters, Moon-sized body, site at the north pole):
//   - CORE: a small plateau whose vertices lie on the sphere R + 2000 m, plus one
//     test facet at its +x edge tilted 10 deg toward +x (the Sun's azimuth).
//   - RING: a flat annulus on the sphere R (2 km below the core), 2-20 km from
//     the site. Open patch, so the analytic body-sphere pre-check is active.
//
// Sun directions from the test facet (depression below the local horizontal,
// toward +x), one per time step:
//   t0  -5.0 deg (Sun 5 deg ABOVE the horizontal): lit control, sets the scale.
//   t1  +0.5 deg: the ray's closest approach is ~66 m below the facet, so it
//       dips under the core-only sphere (R + ~2000 m) but stays ~1.9 km above
//       the ring and every other triangle. Modeled terrain lets the Sun in.
//       PRE-FIX: culled by the core-only sphere -> 0 W/m^2 (the PTM-780 false
//       shadow). POST-FIX: R_occ = ring minimum ~R -> lit, flux ~ S*sin(9.5 deg).
//   t2  +3.0 deg: closest approach ~380 m below R, about 91 km out, well past
//       the 20 km ring. The ring does not block it; only the sphere does, and
//       it must keep doing so after the fix (the sphere still stands in for
//       the body beyond the modeled terrain). Dark both before and after.
//
// Arms: {with ring, without ring} x {batch path (calculateIlluminationOptimized,
// the kernel production uses), single-face path (calculateIllumination)}.
// Without a ring, R_occ comes from the core in both builds, so t1 must stay
// dark there: the fix acts only through the ring.
//
// SEAM case (the companion seal). Same core, but the ring is FLUSH with the core
// (radius R + 2000 m) and starts 600 m out, leaving a 300 m seam gap past the
// 300 m plateau edge, like the registered AGT rings; one ring sector 15 km away
// sits 5 km lower, so min(core, ring) is
// far below the core. With the Sun 2 deg below the horizontal the facet's ray
// dips into the gap ~140 m out and then runs UNDER the flush ring sheet, never
// touching a triangle; its closest approach (R + ~945 m) stays above the lowered
// sphere. Physically it is underground: it must be dark.
//   pre-fix engine:      dark (the core-only sphere catches it)
//   R_occ fix alone:     LIT  (the leak this case exists to catch)
//   R_occ + seam curtain: dark (the ring's inner-loop curtain blocks it)
//
// No SPICE, no GPU: positions come from a precomputed file written next to the
// binary's working directory (removed on exit). Exit 0 = PASS, 1 = FAIL.
#include "illumination/illumination_model.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <unistd.h>

using ptm::IlluminationModel;

namespace {

using V3 = std::array<double, 3>;
using F3 = std::array<int, 3>;

constexpr double kR       = 1737400.0;     // body radius (m)
constexpr double kCoreUp  = 2000.0;        // plateau height above R (m)
constexpr double kAU_km   = 1.495978707e8;
constexpr double kDeg     = M_PI / 180.0;
constexpr double kTiltDeg = 10.0;          // test facet tilt toward +x

int g_fail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++g_fail;
}

V3 onSphere(double x, double y, double r) {
    return {x, y, std::sqrt(r * r - x * x - y * y)};
}

// Plateau grid (on radius R + kCoreUp) plus one tilted test facet. Returns the
// index of the test facet in `faces`.
int buildCore(std::vector<V3>& V, std::vector<F3>& F) {
    const double rc = kR + kCoreUp;
    const int n = 4;                 // 4 x 4 cells, 150 m each, x,y in [-300, 300]
    const double h = 150.0, x0 = -300.0;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
            V.push_back(onSphere(x0 + i * h, x0 + j * h, rc));
    auto id = [&](int i, int j) { return j * (n + 1) + i; };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            F.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            F.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    // Test facet: a 20 m triangle centered at (330, 0), 5 m above the plateau
    // sphere, with normal (sin t, 0, cos t). The in-plane axes are
    // u = (cos t, 0, -sin t) and v = (0, 1, 0); counter-clockwise about n.
    const double t = kTiltDeg * kDeg;
    const V3 c = {330.0, 0.0, std::sqrt(rc * rc - 330.0 * 330.0) + 5.0};
    const V3 u = {std::cos(t), 0.0, -std::sin(t)};
    const V3 v = {0.0, 1.0, 0.0};
    auto P = [&](double a, double b) {
        return V3{c[0] + a * u[0] + b * v[0], c[1] + a * u[1] + b * v[1],
                  c[2] + a * u[2] + b * v[2]};
    };
    const int b0 = static_cast<int>(V.size());
    V.push_back(P(-10.0, -10.0));
    V.push_back(P(10.0, 0.0));
    V.push_back(P(-10.0, 10.0));
    F.push_back({b0, b0 + 1, b0 + 2});
    return static_cast<int>(F.size()) - 1;
}

// Annulus on radius `rsph`, 36 sectors, at the given ground radii. If pit_depth
// > 0, the vertices of the sector around azimuth 180 deg beyond 10 km are lowered
// by pit_depth (a distant low basin that sets min(core, ring)).
void buildRingAt(std::vector<V3>& V, std::vector<F3>& F, const std::vector<double>& radii,
                 double rsph, double pit_depth) {
    const int na = 36;
    for (double rr : radii)
        for (int k = 0; k < na; ++k) {
            const double a = 2.0 * M_PI * k / na;
            const bool pit = pit_depth > 0.0 && rr >= 10e3 && (k >= 16 && k <= 20);
            V.push_back(onSphere(rr * std::cos(a), rr * std::sin(a), rsph - (pit ? pit_depth : 0.0)));
        }
    auto id = [&](int ir, int k) { return ir * na + (k % na); };
    for (int ir = 0; ir + 1 < static_cast<int>(radii.size()); ++ir)
        for (int k = 0; k < na; ++k) {
            F.push_back({id(ir, k), id(ir + 1, k), id(ir + 1, k + 1)});
            F.push_back({id(ir, k), id(ir + 1, k + 1), id(ir, k + 1)});
        }
}

// Flat annulus on radius R, 2-20 km from the pole (2 km below the core).
void buildRing(std::vector<V3>& V, std::vector<F3>& F) {
    buildRingAt(V, F, {2e3, 4e3, 7e3, 10e3, 14e3, 20e3}, kR, 0.0);
}

// Flush annulus on the core sphere, 0.6-20 km (a 300 m seam gap past the plateau
// edge), with a 5 km-deep basin sector beyond 10 km.
void buildSeamRing(std::vector<V3>& V, std::vector<F3>& F) {
    buildRingAt(V, F, {600.0, 1e3, 2e3, 4e3, 7e3, 10e3, 14e3, 20e3}, kR + kCoreUp, 5000.0);
}

double minCentroidRadius(const std::vector<V3>& V, const std::vector<F3>& F) {
    double m = 1e300;
    for (const auto& f : F) {
        double c[3];
        for (int k = 0; k < 3; ++k) c[k] = (V[f[0]][k] + V[f[1]][k] + V[f[2]][k]) / 3.0;
        m = std::min(m, std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]));
    }
    return m;
}

// Sun directions: depression d below the horizontal toward +x.
const std::vector<double> kDepressionDeg = {-5.0, 0.5, 3.0};

// Seam case: t0 lit control (Sun 5 deg up), t1 the underground leak (2 deg down).
const std::vector<double> kSeamDepressionDeg = {-5.0, 2.0};

std::string writePositions(const std::vector<double>& dep, const char* tag) {
    const std::string path = std::string("test_rocc_core_plus_ring_positions_") + tag + "_" +
                             std::to_string(static_cast<long>(getpid())) + ".csv";
    std::ofstream o(path);
    o << "# idx,et,sun_x,sun_y,sun_z,primary_x,primary_y,primary_z (km)\n";
    o.precision(17);
    for (size_t t = 0; t < dep.size(); ++t) {
        const double d = dep[t] * kDeg;
        // Primary placed off to +y, nowhere near any Sun line of sight.
        o << t << "," << 1.0e9 + 3600.0 * t << "," << kAU_km * std::cos(d) << ",0,"
          << -kAU_km * std::sin(d) << ",0,384400,0\n";
    }
    return path;
}

// Compile-time detection of the post-fix getters, so the SAME source builds
// against the pre-fix engine and fails there on behavior, not on compilation.
template <typename T, typename = void>
struct HasRoccGetter : std::false_type {};
template <typename T>
struct HasRoccGetter<T, std::void_t<decltype(std::declval<const T&>().bodyOcclusionRadius()),
                                    decltype(std::declval<const T&>().bodyOcclusionRadiusSource())>>
    : std::true_type {};

template <typename M>
void checkRoccGetter(const M& m, double expect, const char* expect_src, const std::string& arm) {
    if constexpr (HasRoccGetter<M>::value) {
        const double got = m.bodyOcclusionRadius();
        check(std::fabs(got - expect) < 1e-6,
              arm + ": R_occ = " + std::to_string(got) + " m, expected " + std::to_string(expect));
        check(m.bodyOcclusionRadiusSource() == expect_src,
              arm + ": R_occ source = '" + m.bodyOcclusionRadiusSource() + "', expected '" +
                  expect_src + "'");
    } else {
        std::printf("  [INFO] %s: no R_occ getter in this build (pre-fix engine)\n", arm.c_str());
    }
}

struct ArmResult { std::vector<double> flux; bool ok = false; };

ArmResult runArm(bool with_ring, bool batch, const std::string& pos_file, size_t nt,
                 const std::vector<V3>& coreV, const std::vector<F3>& coreF, int test_face,
                 const std::vector<V3>& ringV, const std::vector<F3>& ringF,
                 double expect_rocc, const char* expect_src, const std::string& arm) {
    IlluminationModel m(false);
    m.setShowProgress(false);
    m.setIncludePrimaryReflection(false);
    m.setTimeSteps(static_cast<int>(nt));
    if (with_ring) m.setOccluderGeometry(ringV, ringF);
    ArmResult r;
    if (!m.loadPositionsFromFile(pos_file)) return r;
    const bool ok = batch ? m.calculateIlluminationOptimized(coreV, coreF, "", 1.0)
                          : m.calculateIllumination(coreV, coreF, "", 1.0);
    if (!ok) return r;
    r.ok = true;
    for (size_t t = 0; t < nt; ++t)
        r.flux.push_back(batch ? m.getDirectSolarValue(test_face, t)
                               : m.getIlluminationValue(test_face, t));
    checkRoccGetter(m, expect_rocc, expect_src, arm);
    return r;
}

}  // namespace

int main() {
    std::vector<V3> coreV, ringV;
    std::vector<F3> coreF, ringF;
    const int test_face = buildCore(coreV, coreF);
    buildRing(ringV, ringF);
    const double rmin_core = minCentroidRadius(coreV, coreF);
    const double rmin_ring = minCentroidRadius(ringV, ringF);
    std::printf("test_rocc_core_plus_ring: core %zu faces (min centroid r = R%+.2f m), "
                "ring %zu faces (min centroid r = R%+.2f m)\n",
                coreF.size(), rmin_core - kR, ringF.size(), rmin_ring - kR);

    const std::string pos = writePositions(kDepressionDeg, "base");
    // Expected flux ratio t1/t0 when lit: mu1/mu0 = sin(10-0.5)/sin(10+5).
    const double ratio_lit = std::sin((kTiltDeg - 0.5) * kDeg) / std::sin((kTiltDeg + 5.0) * kDeg);

    for (int batch = 1; batch >= 0; --batch) {
        const std::string path = batch ? "batch path" : "single-face path";
        for (int ring = 1; ring >= 0; --ring) {
            const std::string arm = path + (ring ? ", core + ring" : ", core only");
            std::printf("%s\n", arm.c_str());
            ArmResult r = runArm(ring, batch, pos, kDepressionDeg.size(), coreV, coreF, test_face, ringV, ringF,
                                 ring ? rmin_ring : rmin_core, ring ? "ring" : "core", arm);
            check(r.ok, arm + ": illumination ran");
            if (!r.ok) continue;
            std::printf("    flux t0 (Sun 5 deg up)   = %.4f W/m^2\n", r.flux[0]);
            std::printf("    flux t1 (0.5 deg down)   = %.4f W/m^2\n", r.flux[1]);
            std::printf("    flux t2 (3 deg down)     = %.4f W/m^2\n", r.flux[2]);
            check(r.flux[0] > 100.0, arm + ": t0 lit control is lit");
            if (ring) {
                const double ratio = r.flux[1] / r.flux[0];
                check(std::fabs(ratio / ratio_lit - 1.0) < 0.02,
                      arm + ": t1 clears all modeled terrain -> lit, flux ratio " +
                          std::to_string(ratio) + " vs mu ratio " + std::to_string(ratio_lit) +
                          " (PTM-780: pre-fix engine gives 0)");
            } else {
                check(r.flux[1] == 0.0,
                      arm + ": t1 without a ring stays dark (sphere from the core, unchanged)");
            }
            check(r.flux[2] == 0.0,
                  arm + ": t2 passes below R beyond the ring -> the sphere still culls it");
        }
    }
    std::remove(pos.c_str());

    // ---- SEAM case -------------------------------------------------------
    std::vector<V3> seamV;
    std::vector<F3> seamF;
    buildSeamRing(seamV, seamF);
    const double rmin_seam = std::min(rmin_core, minCentroidRadius(seamV, seamF));
    std::printf("seam case: flush ring %zu faces from 600 m (300 m gap past the 300 m plateau edge), "
                "min centroid r over core + ring = R%+.2f m\n", seamF.size(), rmin_seam - kR);
    const std::string pos_seam = writePositions(kSeamDepressionDeg, "seam");
    for (int batch = 1; batch >= 0; --batch) {
        const std::string arm = std::string(batch ? "batch path" : "single-face path") +
                                ", seam (flush ring + gap)";
        std::printf("%s\n", arm.c_str());
        ArmResult r = runArm(true, batch, pos_seam, kSeamDepressionDeg.size(), coreV, coreF,
                             test_face, seamV, seamF, rmin_seam, "ring", arm);
        check(r.ok, arm + ": illumination ran");
        if (!r.ok) continue;
        std::printf("    flux t0 (Sun 5 deg up)   = %.4f W/m^2\n", r.flux[0]);
        std::printf("    flux t1 (2 deg down)     = %.4f W/m^2\n", r.flux[1]);
        check(r.flux[0] > 100.0, arm + ": t0 lit control is lit");
        check(r.flux[1] == 0.0,
              arm + ": t1 dips into the seam gap and runs under the ring -> dark "
                    "(R_occ fix without the seam curtain lights it)");
    }
    std::remove(pos_seam.c_str());

    std::printf("test_rocc_core_plus_ring: %s (%d failure%s)\n", g_fail ? "FAILED" : "PASSED",
                g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
