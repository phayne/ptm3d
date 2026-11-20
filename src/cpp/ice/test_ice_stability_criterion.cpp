// test_ice_stability_criterion.cpp - analytic gates for ice::computeIceStabilityDepths.
//
// No mesh, no solver, no SPICE: synthetic accumulators whose depth answer is
// known in closed form, so a slip in the criterion itself (a unit error in J0,
// a wrong prefactor, a surface test that only fires on underflow) fails here.
//
//   A  The threshold in physical units: the equivalent surface temperature of
//      J0 = 3e-14 (1 m of ice per Gyr) is 114.4 K and of 1 kg m^-2 Gyr^-1
//      (Schorghofer & Taylor 2007; Paige et al. 2010) 101.3 K, and the engine's
//      free-sublimation mass flux at that temperature returns J0.
//   B  Surface criterion: an isothermal column 0.05 K BELOW the equivalent
//      temperature is stable at the surface (z* = 0); 0.05 K ABOVE it is not,
//      and its depth is the closed form z* = mu l E / (2 J0) (about l/2).
//      Before 2026-11 the colder column got z1/512 (the spurious micrometers).
//   C  S&T (2007) Eq. 11 at depth, isothermal column: <E> constant, so
//      z* = mu l E(T) / (2 J0) exactly; T is chosen so that z* = 0.10 m.
//   D  S&T Eq. 11 with an exponentially damped <E(z)> = E0 exp(-z/H) (the
//      log-linear interpolation is exact for it): z* solves
//      z = (mu l / 2 J0) E0 exp(-z/H), solved here independently by Newton.
//   E  Censored: J > J0 at every node -> z* = z_max.
//   F  Underflow: <E(0)> = 0 -> surface (the legacy rule still holds).
//   G  Invalid J0 or l (0, negative, NaN, inf) -> empty result, never a depth.
// Tolerance for B-D: the solver brackets between nodes and bisects 8 times,
// so |z - z*| <= (z_i - z_{i-1}) / 2^9 (plus float rounding).

#include "ice_table.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace {

int g_fail = 0;

void check(bool ok, const char* name, const char* fmt, double a, double b) {
    std::printf("  %-3s ", name);
    std::printf(fmt, a, b);
    std::printf("  %s\n", ok ? "PASS" : "FAIL");
    if (!ok) ++g_fail;
}

// Geometric grid like the engine's: z0 = 0, z1 = 3.3 mm, 23 nodes to ~0.47 m.
std::vector<double> grid() {
    std::vector<double> z(23, 0.0);
    double dz = 3.3e-3;
    for (int i = 1; i < 23; ++i) { z[i] = z[i - 1] + dz; dz *= 1.17; }
    return z;
}

// One face; <E(z_i)> = E[i] (number flux), n_samples = 1.
ice::IceStabilityAccumulator column(const std::vector<double>& E) {
    const std::vector<double> z = grid();
    ice::IceStabilityAccumulator acc;
    acc.allocate(1, static_cast<int>(z.size()));
    acc.z = z;
    acc.n_samples = 1;
    for (std::size_t i = 0; i < z.size(); ++i) acc.sum_E[i] = static_cast<float>(E[i]);
    return acc;
}

double solve1(const std::vector<double>& E, double J0) {
    auto out = ice::computeIceStabilityDepths(column(E), J0, 75e-6, -1.0, false);
    return out.size() == 1 ? out[0] : std::numeric_limits<double>::quiet_NaN();
}

// Cell width around depth zs on the test grid, for the bisection tolerance.
double cell(double zs) {
    const std::vector<double> z = grid();
    for (std::size_t i = 1; i < z.size(); ++i)
        if (zs <= z[i]) return z[i] - z[i - 1];
    return z.back() - z[z.size() - 2];
}

std::vector<double> isothermal(double T) {
    const double E = ice::sublimationMassFlux(T) / ice::kWaterMoleculeMassKg;
    return std::vector<double>(grid().size(), E);
}

}  // namespace

int main() {
    const double MU = ice::kWaterMoleculeMassKg, ELL = 75e-6;
    const double J_code = 3.0e-14, J_lit = 1.0 / ice::kSecondsPerGyr;
    std::printf("test_ice_stability_criterion\n");
    std::printf("  %s\n  %s\n", ice::iceThresholdSummary(J_code).c_str(),
                ice::iceThresholdSummary(J_lit).c_str());

    // A
    const double T_code = ice::surfaceStabilityTemperature(J_code);
    const double T_lit  = ice::surfaceStabilityTemperature(J_lit);
    check(std::fabs(T_code - 114.4) < 0.05 && std::fabs(T_lit - 101.3) < 0.05, "A",
          "T_eq(3e-14) = %.3f K, T_eq(1 kg/m2/Gyr) = %.3f K", T_code, T_lit);
    check(std::fabs(ice::sublimationMassFlux(T_code) / J_code - 1) < 1e-9 &&
          std::fabs(ice::sublimationMassFlux(T_lit) / J_lit - 1) < 1e-9, "A2",
          "mu E(T_eq) / J0 - 1 = %.1e, %.1e",
          ice::sublimationMassFlux(T_code) / J_code - 1, ice::sublimationMassFlux(T_lit) / J_lit - 1);

    for (double J0 : {J_code, J_lit}) {
        const double Teq = ice::surfaceStabilityTemperature(J0);
        // B
        const double z_cold = solve1(isothermal(Teq - 0.05), J0);
        const double z_warm = solve1(isothermal(Teq + 0.05), J0);
        const double zs_warm = 0.5 * MU * ELL * (isothermal(Teq + 0.05)[0]) / J0;
        check(z_cold == 0.0, "B", "T_eq-0.05 K: z* = %.3g m (want 0), J0 = %.3g", z_cold, J0);
        check(z_warm > 0.0 && std::fabs(z_warm - zs_warm) <= cell(zs_warm) / 512 * 1.01, "B2",
              "T_eq+0.05 K: z* = %.4g m, closed form %.4g m", z_warm, zs_warm);
        // C
        const double target = 0.10;
        const double T_c = ice::surfaceStabilityTemperature(J0 * 2.0 * target / ELL);
        const double z_c = solve1(isothermal(T_c), J0);
        check(std::fabs(z_c - target) <= cell(target) / 512 * 1.01 + 1e-7, "C",
              "isothermal at %.3f K: z* = %.6f m (want 0.100000)", T_c, z_c);
        // D: E0 at a warm surface, H = 3 cm
        const double H = 0.03, E0 = 400.0 * J0 / MU;   // mu E0 = 400 J0
        std::vector<double> Ed;
        for (double z : grid()) Ed.push_back(E0 * std::exp(-z / H));
        double zs = 0.05;                                // Newton on g(z) = z - P e^{-z/H}
        const double P = 0.5 * MU * ELL * E0 / J0;
        for (int k = 0; k < 100; ++k) zs -= (zs - P * std::exp(-zs / H)) / (1 + P / H * std::exp(-zs / H));
        const double z_d = solve1(Ed, J0);
        check(std::fabs(z_d - zs) <= cell(zs) / 512 * 1.01 + 1e-7, "D",
              "E0 e^{-z/H}: z* = %.6f m, Eq. 11 root %.6f m", z_d, zs);
        // E
        const double z_e = solve1(std::vector<double>(grid().size(), 1e6 * J0 / MU), J0);
        check(z_e == static_cast<float>(grid().back()), "E", "censored: z* = %.6f m, z_max %.6f m",
              z_e, grid().back());
    }
    // F
    check(solve1(std::vector<double>(grid().size(), 0.0), J_code) == 0.0, "F",
          "<E(0)> = 0: z* = %.1f (want 0)%.0f", 0.0, 0.0);
    // G
    int n_empty = 0;
    for (double J : {0.0, -3e-14, std::nan(""), std::numeric_limits<double>::infinity()})
        n_empty += ice::computeIceStabilityDepths(column(isothermal(110.0)), J, ELL, -1, false).empty();
    n_empty += ice::computeIceStabilityDepths(column(isothermal(110.0)), J_code, 0.0, -1, false).empty();
    check(n_empty == 5, "G", "invalid J0 / l refused: %.0f of %.0f", n_empty, 5.0);

    if (g_fail) { std::printf("[test_ice_stability_criterion] FAIL (%d)\n", g_fail); return 1; }
    std::printf("[test_ice_stability_criterion] PASS\n");
    return 0;
}
