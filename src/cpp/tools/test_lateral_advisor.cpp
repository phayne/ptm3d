// test_lateral_advisor -- the lateral-conduction advisor (thermal/lateral_advisor.h).
//
//   [1] D*(T) reproduces the Moon3D paper's quoted values for lunar surface regolith
//       (Hayne et al., in revision, Section 2): ~0.9 mm sunlit, 6.5 mm at 110 K, ~2 cm at 70 K,
//       ~4 cm at 50 K. Also pins the computed values, so a change in the formula or the inputs
//       shows up here as a number, not as a changed sentence in a log.
//   [2] The advisor's D* equals the engine's own ptm::lateral::conductionLength (the function
//       the lateral edge-activation rule calls), so the two cannot drift apart.
//   [3] D* decreases monotonically with temperature (40-400 K).
//   [4] Verdicts on synthetic meshes: a Moon3D-like millimeter mesh -> RECOMMENDED; a
//       Haworth-like 250 m mesh -> not needed; a 1 m shape-from-shading-like mesh -> not needed
//       (with the rock caveat); a meter mesh with a 0.5%-area patch of 2 cm facets -> BORDERLINE;
//       an empty mesh -> UNMEASURED (warns, never silent).
//   [5] Formatting: the lateral-ON variants, and the silence key named in every block.
//
// Links libthermal only for [2]; everything else is the header under test.
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "thermal/lateral_advisor.h"
#include "thermal/lateral_conduction.h"

namespace la = ptm::lateral;

static int g_fail = 0, g_pass = 0;
static void check(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    (ok ? g_pass : g_fail)++;
}

// Moon, Hayne et al. 2017 Table A1, as planetary_properties.h getMoonProperties() sets it.
static la::SurfaceMaterial lunarRegolith() {
    la::SurfaceMaterial m;
    m.ks = 7.4e-4;
    m.rhos = 1100.0;
    m.cp0 = -3.6125;
    m.cpCoeff = {2.7431, 2.3616e-3, -1.2340e-5, 8.9093e-9};
    m.emissivity = 0.95;
    m.chi = 2.7;
    m.R350 = 2.7 / (350.0 * 350.0 * 350.0);
    m.model = thermal::ConductivityModel::Hayne2017;
    m.day_s = 29.530588 * 24 * 3600;
    return m;
}

// A flat square patch of side L tiled by right isosceles triangles with legs s
// (mean edge = (2 + sqrt 2) s / 3 = 1.138 s), in meters.
struct Grid {
    std::vector<std::array<double, 3>> V;
    std::vector<std::array<int, 3>> F;
};
static void addGrid(Grid& g, double L, double s, double x0 = 0.0) {
    const int n = static_cast<int>(std::lround(L / s));
    const int base = static_cast<int>(g.V.size());
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) g.V.push_back({x0 + i * s, j * s, 0.0});
    auto id = [&](int i, int j) { return base + j * (n + 1) + i; };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            g.F.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            g.F.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
}

static bool near(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); }

static void printBlock(const std::vector<std::string>& lines) {
    for (const auto& l : lines) std::printf("      %s\n", l.c_str());
}

int main() {
    const la::SurfaceMaterial moon = lunarRegolith();

    std::printf("[1] D*(T) on lunar surface regolith vs the Moon3D paper\n");
    const double d300 = la::criticalLength(moon, 300.0), d380 = la::criticalLength(moon, 380.0);
    const double d110 = la::criticalLength(moon, 110.0), d70 = la::criticalLength(moon, 70.0);
    const double d50 = la::criticalLength(moon, 50.0);
    std::printf("      D*(380 K) = %.4f mm, D*(300 K) = %.4f mm, D*(110 K) = %.4f mm, D*(70 K) = %.3f cm, "
                "D*(50 K) = %.3f cm\n", d380 * 1e3, d300 * 1e3, d110 * 1e3, d70 * 1e2, d50 * 1e2);
    check(std::lround(d300 * 1e4) == 9, "sunlit (300 K) rounds to the paper's 0.9 mm");
    check(std::lround(d110 * 1e4) == 65, "110 K rounds to the paper's 6.5 mm");
    check(std::lround(d70 * 1e2) == 2, "70 K rounds to the paper's 2 cm");
    check(std::lround(d50 * 1e2) == 4, "50 K rounds to the paper's 4 cm");
    check(d380 < d300 && d380 > 0.5e-3, "380 K (lunar noon) is below the 300 K value but still sub-millimeter-scale");
    // Regression pins (independent float64 Python evaluation of the same formula and inputs,
    // 2026-11-02): 0.8650237 mm, 6.528925 mm, 2.049285 cm, 4.353098 cm.
    check(near(d300, 8.650237e-4, 1e-6) && near(d110, 6.528925e-3, 1e-6) && near(d70, 2.049285e-2, 1e-6) &&
              near(d50, 4.353098e-2, 1e-6),
          "pinned values to 1e-6: 0.8650 mm (300 K), 6.5289 mm (110 K), 2.0493 cm (70 K), 4.3531 cm (50 K)");

    std::printf("[2] the advisor's D* equals the engine's ptm::lateral::conductionLength\n");
    double worst = 0.0;
    for (double T = 30.0; T <= 420.0; T += 5.0) {
        const double k = la::surfaceConductivity(moon, T);
        const double cp = la::surfaceHeatCapacity(moon, T);
        const double e = la::conductionLength(k, moon.rhos, cp, moon.emissivity, T, 2.0 * M_PI / moon.day_s);
        worst = std::max(worst, std::fabs(la::criticalLength(moon, T) / e - 1.0));
    }
    std::printf("      worst relative difference 30-420 K: %.3g\n", worst);
    check(worst < 1e-12, "agrees to round-off (< 1e-12) at every temperature");

    std::printf("[3] D* decreases with temperature\n");
    bool mono = true;
    for (double T = 40.0; T < 400.0; T += 1.0)
        if (!(la::criticalLength(moon, T + 1.0) < la::criticalLength(moon, T))) mono = false;
    check(mono, "strictly decreasing, 40-400 K");

    std::printf("[4] verdicts on synthetic meshes\n");
    la::AdvisorContext off;
    {
        Grid g; addGrid(g, 1.3, 2.2e-3);   // Moon3D-like: 1.3 m patch, ~2.5 mm mean edge
        const auto r = la::evaluateLateralAdvice(la::measureFacetSizes(g.V, g.F, 1.0), moon);
        std::printf("    Moon3D-like (1.3 m patch, %zu faces):\n", g.F.size());
        printBlock(la::formatLateralAdvice(r, moon, off));
        check(r.verdict == la::AdvisorVerdict::Recommended, "Moon3D-like millimeter mesh -> RECOMMENDED");
        check(r.frac_cold > 0.99 && r.frac_ice > 0.99, "every facet is below D*(110 K) and D*(50 K)");
    }
    {
        Grid g; addGrid(g, 50000.0, 220.0);   // Haworth-like: 50 km, ~250 m facets
        const auto r = la::evaluateLateralAdvice(la::measureFacetSizes(g.V, g.F, 1.0), moon);
        std::printf("    Haworth-like (50 km, %zu faces):\n", g.F.size());
        printBlock(la::formatLateralAdvice(r, moon, off));
        check(r.verdict == la::AdvisorVerdict::NotNeeded, "Haworth-like 250 m mesh -> not needed");
    }
    {
        Grid g; addGrid(g, 300.0, 0.88);   // SfS-like: 300 m, ~1 m facets
        const auto r = la::evaluateLateralAdvice(la::measureFacetSizes(g.V, g.F, 1.0), moon);
        std::printf("    SfS-like (300 m, %zu faces):\n", g.F.size());
        const auto lines = la::formatLateralAdvice(r, moon, off);
        printBlock(lines);
        check(r.verdict == la::AdvisorVerdict::NotNeeded, "1 m SfS-like mesh -> not needed");
        bool rock = false;
        for (const auto& l : lines) rock = rock || l.find("In rock") != std::string::npos;
        check(rock, "the 1 m verdict carries the rock caveat (median facet < 100 D*(50 K))");
    }
    {
        // 1 m facets over 20 m x 20 m, plus a 1.42 m x 1.42 m patch of 2 cm facets beside it
        // (0.5% of the area): below D*(50 K), above 0.1%, below the 1% recommend line.
        Grid g; addGrid(g, 20.0, 0.88); addGrid(g, 1.42, 0.0176, 25.0);
        const auto r = la::evaluateLateralAdvice(la::measureFacetSizes(g.V, g.F, 1.0), moon);
        std::printf("    mixed (1 m mesh with a 0.5%%-area patch of 2 cm facets, %zu faces):\n", g.F.size());
        printBlock(la::formatLateralAdvice(r, moon, off));
        check(r.verdict == la::AdvisorVerdict::Borderline, "a 0.5% fine-facet patch -> BORDERLINE");
    }
    {
        Grid g;
        const auto r = la::evaluateLateralAdvice(la::measureFacetSizes(g.V, g.F, 1.0), moon);
        const auto lines = la::formatLateralAdvice(r, moon, off);
        printBlock(lines);
        check(r.verdict == la::AdvisorVerdict::Unmeasured && lines.size() == 1 &&
                  lines[0].find("did NOT run") != std::string::npos,
              "empty mesh -> UNMEASURED, and it says so");
    }
    {
        // Kilometer units: the same 250 m grid written in km must give the same answer at to_m = 1000.
        Grid gm; addGrid(gm, 50000.0, 220.0);
        Grid gk; addGrid(gk, 50.0, 0.22);
        const auto rm = la::evaluateLateralAdvice(la::measureFacetSizes(gm.V, gm.F, 1.0), moon);
        const auto rk = la::evaluateLateralAdvice(la::measureFacetSizes(gk.V, gk.F, 1000.0), moon);
        check(rk.verdict == rm.verdict && near(rk.size_median, rm.size_median, 1e-9),
              "a mesh in kilometers with to_m = 1000 matches the same mesh in meters");
    }

    std::printf("[5] formatting\n");
    {
        la::AdvisorContext on; on.lateral_on = true;
        Grid fine; addGrid(fine, 0.3, 2.2e-3);
        Grid coarse; addGrid(coarse, 5000.0, 220.0);
        const auto lf = la::formatLateralAdvice(
            la::evaluateLateralAdvice(la::measureFacetSizes(fine.V, fine.F, 1.0), moon), moon, on);
        const auto lc = la::formatLateralAdvice(
            la::evaluateLateralAdvice(la::measureFacetSizes(coarse.V, coarse.F, 1.0), moon), moon, on);
        std::printf("    lateral ON, fine mesh:\n"); printBlock({lf[0]});
        std::printf("    lateral ON, coarse mesh:\n"); printBlock({lc[0]});
        check(lf[0].find("is ON, and appropriate") != std::string::npos, "ON + fine mesh: 'ON, and appropriate'");
        check(lc[0].find("ON but unnecessary") != std::string::npos, "ON + coarse mesh: 'ON but unnecessary'");
        bool silence = false, prefix = true;
        for (const auto& l : lf) {
            silence = silence || l.find("thermal.lateral_advisor: false") != std::string::npos;
            prefix = prefix && l.rfind("[lateral-advisor] ", 0) == 0;
        }
        check(silence && prefix, "every line carries the [lateral-advisor] prefix; the block names the silence key");
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    std::printf("%s\n", g_fail ? "FAIL" : "PASS: all lateral-advisor checks");
    return g_fail ? 1 : 0;
}
