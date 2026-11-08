// PTM-819 regression test: primary-disk (Earthshine) rays must take the same
// open-mesh body-sphere cull (R_occ) as the direct Sun.
//
// Geometry (meters, Moon-sized body, site at the north pole, no occluder ring):
//   - CORE: a flat 600 m plateau on the sphere R (open patch, so the analytic
//     body pre-check is active and R_occ ~ R).
//   - facet A: 20 m triangle 5 m above the plateau, tilted 10 deg toward +x.
//   - facet B: same, tilted 3.1 deg toward +x (a grazing facet; see below).
// Earth (radius 6371 km, 384,400 km away, angular radius ~0.95 deg) is placed
// toward +x at a depression d below the local horizontal, one epoch per d:
//   t0  d = -5.0  (whole disk ABOVE the limb): lit control. Must be BIT-IDENTICAL
//                 with the cull on and off.
//   t1  d = +3.0  (whole disk 2.05..3.95 deg BELOW the limb). Facets A and B are
//                 tilted far enough to "see" it above their own planes, and the
//                 tiny patch has no terrain beyond, so a bare Embree ray escapes.
//                 PRE-FIX (cull off): Earth IR and reflected light leak in.
//                 POST-FIX: exactly 0 on both facets.
//   t2  d = +0.5  (disk straddles the limb, top 0.45 deg above): partial. With
//                 the cull, 0 < flux < the uncut flux.
// Facet B at t1 is priced by the per-face budget (armed with the Moon's own
// emissivity/floor) at the NO-RAY rung (it could absorb the whole disk within
// 0.163 K); the test checks that it IS that rung, so the rung-0 branch of the
// fix is exercised, not just the ray branch.
//
// Arms: cull on (default) / off (illumination.primary_disk_body_cull: false).
// The same source builds against a pre-fix engine (no setter): there the "on"
// arm IS the pre-fix engine and fails t1. No SPICE, no GPU: positions come from
// a precomputed file (removed on exit). Exit 0 = PASS, 1 = FAIL.
#include "illumination/illumination_model.h"
#include "illumination/primary_disk.h"

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

constexpr double kR      = 1737400.0;      // body radius (m)
constexpr double kAU_km  = 1.495978707e8;
constexpr double kEarthD = 384400.0;       // km
constexpr double kDeg    = M_PI / 180.0;
constexpr double kEps    = 0.95;           // Moon emissivity (budget arming)
constexpr double kQb     = 0.018;          // Moon deep-heat floor (budget arming)

int g_fail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++g_fail;
}

V3 onSphere(double x, double y, double r) { return {x, y, std::sqrt(r * r - x * x - y * y)}; }

// Tilted 20 m test facet centered at (x, y) 5 m above the plateau; returns its index.
int addFacet(std::vector<V3>& V, std::vector<F3>& F, double x, double y, double tilt_deg) {
    const double t = tilt_deg * kDeg;
    const V3 c = {x, y, std::sqrt(kR * kR - x * x - y * y) + 5.0};
    const V3 u = {std::cos(t), 0.0, -std::sin(t)};
    const V3 v = {0.0, 1.0, 0.0};
    auto P = [&](double a, double b) {
        return V3{c[0] + a * u[0] + b * v[0], c[1] + a * u[1] + b * v[1], c[2] + a * u[2] + b * v[2]};
    };
    const int b0 = static_cast<int>(V.size());
    V.push_back(P(-10.0, -10.0));
    V.push_back(P(10.0, 0.0));
    V.push_back(P(-10.0, 10.0));
    F.push_back({b0, b0 + 1, b0 + 2});
    return static_cast<int>(F.size()) - 1;
}

void buildCore(std::vector<V3>& V, std::vector<F3>& F, int& fa, int& fb) {
    const int n = 4;
    const double h = 150.0, x0 = -300.0;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) V.push_back(onSphere(x0 + i * h, x0 + j * h, kR));
    auto id = [&](int i, int j) { return j * (n + 1) + i; };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            F.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            F.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    fa = addFacet(V, F, 330.0, 0.0, 10.0);
    fb = addFacet(V, F, 330.0, 100.0, 3.1);
}

V3 centroidKm(const std::vector<V3>& V, const F3& f) {
    V3 c{};
    for (int k = 0; k < 3; ++k) c[k] = (V[f[0]][k] + V[f[1]][k] + V[f[2]][k]) / 3000.0;
    return c;
}

const std::vector<double> kDepressionDeg = {-5.0, 3.0, 0.5};

// Earth toward +x at depression d below the horizontal at `site_km` (local up ~ +z
// at the pole); the Sun 30 deg up toward +y, far from every Earth line of sight.
std::vector<std::pair<V3, V3>> positions(const V3& site_km) {
    std::vector<std::pair<V3, V3>> out;
    for (double dd : kDepressionDeg) {
        const double d = dd * kDeg;
        const V3 e = {site_km[0] + kEarthD * std::cos(d), site_km[1], site_km[2] - kEarthD * std::sin(d)};
        const V3 s = {0.0, kAU_km * std::cos(30.0 * kDeg), kAU_km * std::sin(30.0 * kDeg)};
        out.push_back({s, e});
    }
    return out;
}

std::string writePositions(const std::vector<std::pair<V3, V3>>& sp) {
    const std::string path = "test_primary_disk_body_cull_positions_" +
                             std::to_string(static_cast<long>(getpid())) + ".csv";
    std::ofstream o(path);
    o << "# idx,et,sun_x,sun_y,sun_z,primary_x,primary_y,primary_z (km)\n";
    o.precision(17);
    for (size_t t = 0; t < sp.size(); ++t)
        o << t << "," << 1.0e9 + 3600.0 * t << "," << sp[t].first[0] << "," << sp[t].first[1] << ","
          << sp[t].first[2] << "," << sp[t].second[0] << "," << sp[t].second[1] << "," << sp[t].second[2]
          << "\n";
    return path;
}

template <typename T, typename = void>
struct HasCullSetter : std::false_type {};
template <typename T>
struct HasCullSetter<T, std::void_t<decltype(std::declval<T&>().setPrimaryDiskBodyCull(true))>>
    : std::true_type {};

// Sets the cull if this engine has the setter; returns false on a pre-fix engine.
template <typename M>
bool setCull(M& m, bool on) {
    if constexpr (HasCullSetter<M>::value) {
        m.setPrimaryDiskBodyCull(on);
        return true;
    } else {
        (void)m; (void)on;
        return false;
    }
}

struct Arm { std::vector<double> ir_a, rf_a, ir_b, rf_b; bool ok = false; bool rung0_b = false; };

Arm runArm(bool cull, const std::string& pos, const std::vector<V3>& V, const std::vector<F3>& F,
           int fa, int fb) {
    IlluminationModel m(false);
    m.setShowProgress(false);
    m.setPrimaryBody("EARTH");
    m.setIncludePrimaryReflection(true);
    m.setPrimaryDiskBudget(kEps, kQb);
    if (!setCull(m, cull) && !cull)
        std::printf("  [INFO] pre-fix engine: no cull setter; the 'off' arm equals the 'on' arm\n");
    m.setTimeSteps(static_cast<int>(kDepressionDeg.size()));
    Arm r;
    if (!m.loadPositionsFromFile(pos)) return r;
    if (!m.calculateIlluminationOptimized(V, F, "", 1.0)) return r;
    r.ok = true;
    for (size_t t = 0; t < kDepressionDeg.size(); ++t) {
        r.ir_a.push_back(m.getPrimaryIRValue(fa, t));
        r.rf_a.push_back(m.getPrimaryReflectedValue(fa, t));
        r.ir_b.push_back(m.getPrimaryIRValue(fb, t));
        r.rf_b.push_back(m.getPrimaryReflectedValue(fb, t));
    }
    // Facet B at t1: is the full-disk flux inside the no-ray rung's budget?
    {
        const auto sp = positions(centroidKm(V, F[fb]));
        const V3 fc = centroidKm(V, F[fb]);
        const V3 e = sp[1].second;
        const V3 rel = {e[0] - fc[0], e[1] - fc[1], e[2] - fc[2]};
        const double dist = std::sqrt(rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]);
        const V3 dir = {rel[0] / dist, rel[1] / dist, rel[2] / dist};
        const double t = 3.1 * kDeg;
        const double c0 = std::sin(t) * dir[0] + std::cos(t) * dir[2];
        const double alpha = std::asin(6371.0 / dist);
        const double A = ptm::disk::diskMeanCosAboveHorizon(c0, alpha);
        const double f_full = m.probePrimaryHeating(A, dist) +
                              m.probePrimaryReflected(A, dir, sp[1].first, e, dist);
        const double dv = ptm::disk::primaryDiskAllowedVisibilityError(f_full, kQb, kEps * 5.670374419e-8);
        r.rung0_b = (ptm::disk::primaryDiskSamplesForAllowedError(dv) == 0);
        std::printf("    facet B at t1: c0 = %.5f, A = %.5f, F_full = %.3e W/m^2, budget rung %s\n", c0, A,
                    f_full, r.rung0_b ? "0 (no rays)" : ">0");
    }
    return r;
}

}  // namespace

int main() {
    std::vector<V3> V;
    std::vector<F3> F;
    int fa = -1, fb = -1;
    buildCore(V, F, fa, fb);
    // The Earth geometry is referenced to facet A; facet B is 100 m away, a
    // 1.5e-5 deg parallax at 384,400 km.
    const std::string pos = writePositions(positions(centroidKm(V, F[fa])));
    std::printf("test_primary_disk_body_cull: %zu core faces (open patch), Earth depressions -5 / +3 / +0.5 deg\n",
                F.size());

    Arm on, off;
    std::printf("arm: cull ON (default)\n");
    on = runArm(true, pos, V, F, fa, fb);
    std::printf("arm: cull OFF (pre-fix behavior)\n");
    off = runArm(false, pos, V, F, fa, fb);
    std::remove(pos.c_str());
    check(on.ok && off.ok, "both arms ran");
    if (!(on.ok && off.ok)) return 1;

    for (size_t t = 0; t < kDepressionDeg.size(); ++t)
        std::printf("    t%zu d=%+.1f  A: IR on %.6e off %.6e | refl on %.6e off %.6e    "
                    "B: IR on %.6e off %.6e | refl on %.6e off %.6e\n",
                    t, kDepressionDeg[t], on.ir_a[t], off.ir_a[t], on.rf_a[t], off.rf_a[t], on.ir_b[t],
                    off.ir_b[t], on.rf_b[t], off.rf_b[t]);

    std::printf("    t0 on - off: A IR %.3e refl %.3e | B IR %.3e refl %.3e\n", on.ir_a[0] - off.ir_a[0],
                on.rf_a[0] - off.rf_a[0], on.ir_b[0] - off.ir_b[0], on.rf_b[0] - off.rf_b[0]);
    std::printf("    t0 exact (cross-build identity): on A IR %a refl %a B IR %a refl %a | off A refl %a\n",
                on.ir_a[0], on.rf_a[0], on.ir_b[0], on.rf_b[0], off.rf_a[0]);
    // t0: above the limb -> unchanged, bit for bit.
    check(off.ir_a[0] > 0.0 && off.rf_a[0] > 0.0, "t0 control: Earth above the limb lights facet A");
    check(on.ir_a[0] == off.ir_a[0] && on.rf_a[0] == off.rf_a[0] && on.ir_b[0] == off.ir_b[0] &&
              on.rf_b[0] == off.rf_b[0],
          "t0: Earth above the limb -> cull on == off, bit-identical (both facets, IR + reflected)");
    // t1: whole disk below the limb.
    if (HasCullSetter<IlluminationModel>::value)
        check(off.ir_a[1] > 0.0 && off.ir_b[1] > 0.0,
              "t1 cull OFF reproduces the PTM-819 leak (Earth IR through the body on A and B)");
    check(on.ir_a[1] == 0.0 && on.rf_a[1] == 0.0,
          "t1: Earth wholly below the limb -> facet A gets exactly 0 IR and 0 reflected (probe/sample rays)");
    check(on.rung0_b, "t1: facet B is priced at the no-ray budget rung (the rung-0 branch is exercised)");
    check(on.ir_b[1] == 0.0 && on.rf_b[1] == 0.0,
          "t1: Earth wholly below the limb -> facet B (no-ray rung) gets exactly 0 IR and 0 reflected");
    // t2: straddling the limb -> partial.
    check(on.ir_a[2] > 0.0 && on.ir_a[2] < off.ir_a[2],
          "t2: disk straddles the limb -> 0 < IR(on) < IR(off) on facet A");
    std::printf("test_primary_disk_body_cull: %s (%d failure%s)\n", g_fail ? "FAILED" : "PASSED", g_fail,
                g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
