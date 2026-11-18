// test_ice_retention.cpp - analytic gates for the ice retention criterion
// (ice_retention.h). Synthetic columns only: no mesh, no SPICE, no GPU; about
// a second.
//
//   A  isothermal column, constant m and ell: z(tau) = (ell/2)[(1 + 4 mu E tau/(m ell))^1/2 - 1],
//      at 105-130 K, relative 1e-6 (the solver is exact for log-linear <E>),
//      against the closed form and the reference table to its printed digits;
//      z(4 tau)/z(tau) within 1e-3 of 2 at 125 K; a 0.30 m column at 130 K is
//      "none within the column".
//   B  two layers (m 500 -> 400 kg m^-3, ell 100 -> 30 um at 5 cm), 125 K:
//      the closed form at 0.1, 1 and 3 Gyr, relative 1e-6, and the reference
//      values 4.1546 / 11.0975 / 17.1920 cm.
//   C  <E> = 0 at every node: retained in layer 0 with the underflow class,
//      lifetimes at FLT_MAX, no NaN anywhere.
//   D  grid convergence on a damped-sinusoid column with pore-filling ice on
//      nested geometric grids of 12, 23, 45 and 89 nodes to 0.4706 m: the error
//      falls by at least 3 per halving and is below 0.5% at 23 nodes; the
//      reference depths match 3.302 / 8.016 / 0.797 cm.
//   E  layer index and class bookkeeping (layer boundaries at node midpoints).
//   F  <E> falling by 325, 1e4 and 8e4 across the first node interval: top-layer
//      lifetime, ice-table depth in that interval and column lifetime against
//      an independent Simpson reference, relative 1e-9.
//   G  the grain density comes from the body: the Moon and Mercury their own
//      (silicate, 3100 kg m^-3) or the configured value; a generic body only
//      when configured; Europa (ice regolith) refused even when configured;
//      no grain density at all is refused by the parameter check.

#include "ice_retention.h"
#include "ice_table.h"
#include "../thermal/planetary_properties.h"   // header-only body table

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
constexpr double MU = ice::kWaterMoleculeMassKg;
constexpr double GYR = ice::kSecondsPerGyr;
int g_fail = 0;

void gate(const char* name, bool ok, const char* fmt_detail, double a = 0, double b = 0, double c = 0) {
    char d[256];
    std::snprintf(d, sizeof d, fmt_detail, a, b, c);
    std::printf("  %-3s %-86s %s\n", name, d, ok ? "PASS" : "FAIL");
    if (!ok) ++g_fail;
}

double rate(double T) { return ice::sublimationMassFlux(T) / MU; }    // molecules m^-2 s^-1

std::vector<double> geometric(int n_intervals, double z_bot, double ratio) {
    std::vector<double> z(n_intervals + 1, 0.0);
    const double dz0 = z_bot * (ratio - 1.0) / (std::pow(ratio, n_intervals) - 1.0);
    double dz = dz0;
    for (int i = 1; i <= n_intervals; ++i) { z[i] = z[i - 1] + dz; dz *= ratio; }
    z[n_intervals] = z_bot;
    return z;
}

ice::IceRetentionParams params(double age) {
    ice::IceRetentionParams p;
    p.age_gyr = age;
    p.grain_density_kgm3 = 3100.0;
    p.rho_surface_kgm3 = 1100.0; p.rho_deep_kgm3 = 1800.0; p.density_scale_m = 0.06;
    return p;
}

double solve_const(double T, double m, double ell, double age, const std::vector<double>& z,
                   ice::IceRetentionFace* out = nullptr) {
    const auto col = ice::makeIceColumnFromFunctions(z, [m](double) { return m; },
                                                     [ell](double) { return ell; });
    std::vector<double> E(z.size(), rate(T));
    auto r = ice::solveIceRetentionColumn(z, E.data(), col, params(age));
    if (out) *out = r;
    return r.z_tau_m;
}
} // namespace

int main() {
    std::printf("[test_ice_retention] ice retention criterion, analytic gates\n");
    // The engine's production grid shape: 23 nodes to 0.4706 m, z_1 = 3.308 mm.
    const std::vector<double> zprod = geometric(22, 0.4706301768531821, 1.14201);

    // ---- A: isothermal ------------------------------------------------------
    {
        const double m = 400.0, ell = 75e-6;
        const double Ts[] = {105, 110, 115, 120, 125, 130};
        const double ref_mm[] = {1.201, 4.583, 15.335, 46.212, 127.324, 324.270};
        const double ref_tcol[] = {1.4e5, 1.0e4, 935, 103, 13.6, 2.1};
        double worst = 0, worst_tab = 0, worst_tcol = 0;
        for (int i = 0; i < 6; ++i) {
            const double muE = rate(Ts[i]) * MU;
            const double exact = 0.5 * ell * (std::sqrt(1.0 + 4.0 * muE * GYR / (m * ell)) - 1.0);
            ice::IceRetentionFace fr;
            const double zt = solve_const(Ts[i], m, ell, 1.0, zprod, &fr);
            worst = std::fmax(worst, std::fabs(zt / exact - 1.0));
            // printed digits: 4-6 significant; half a unit in the last place
            const double half_ulp = 0.5e-3 / ref_mm[i];
            worst_tab = std::fmax(worst_tab, std::fabs(zt * 1e3 - ref_mm[i]) / ref_mm[i] / half_ulp);
            const double zb = zprod.back();
            const double tcol = m * (zb + zb * zb / ell) / muE / GYR;
            worst_tcol = std::fmax(worst_tcol, std::fabs(fr.t_column_s / GYR / tcol - 1.0));
            const double tab_rel = std::fabs(tcol / ref_tcol[i] - 1.0);
            if (tab_rel > 0.05) { std::printf("    t(0.47 m) at %g K: %.4g Gyr vs table %.3g\n", Ts[i], tcol, ref_tcol[i]); ++g_fail; }
            std::printf("    %3.0f K  mu E = %.4g  z_tau = %9.4f mm (exact %9.4f, table %8.3f)  t(z_bot) = %.4g Gyr\n",
                        Ts[i], muE, zt * 1e3, exact * 1e3, ref_mm[i], fr.t_column_s / GYR);
        }
        gate("A1", worst < 1e-6, "isothermal z_tau vs closed form, 105-130 K: worst relative %.2e (tol 1e-6)", worst);
        gate("A2", worst_tab <= 1.0, "against the reference table to its printed digits: worst %.2f half-units", worst_tab);
        gate("A3", worst_tcol < 1e-9, "column lifetime t(z_bot) vs closed form: worst relative %.2e", worst_tcol);
        const double z1 = solve_const(125, m, ell, 1.0, zprod), z4 = solve_const(125, m, ell, 4.0, zprod);
        gate("A4", std::fabs(z4 / z1 - 2.0) < 1e-3, "square-root retreat at 125 K: z(4 tau)/z(tau) = %.6f", z4 / z1);
        std::vector<double> z30 = geometric(20, 0.30, 1.14201);
        ice::IceRetentionFace fr;
        solve_const(130, m, ell, 1.0, z30, &fr);
        gate("A5", fr.cls == ice::kNoneWithinColumn && fr.layer == -1 && fr.z_tau_m == 0.30,
             "130 K, z_bot = 0.30 m: none within the column (class %.0f, layer %.0f)", fr.cls, fr.layer);
    }

    // ---- B: two layers --------------------------------------------------------
    {
        const double a = 0.05, m1 = 500, l1 = 100e-6, m2 = 400, l2 = 30e-6, T = 125;
        std::vector<double> z;                       // nodes every cm, so 5 cm is a node
        for (int i = 0; i <= 50; ++i) z.push_back(0.01 * i);
        const auto col = ice::makeIceColumnFromFunctions(
            z, [=](double s) { return s < a ? m1 : m2; }, [=](double s) { return s < a ? l1 : l2; });
        const double muE = rate(T) * MU;
        const double ta = m1 * (a + a * a / l1) / muE;
        const double ages[] = {0.1, 1.0, 3.0}, ref_cm[] = {4.1546, 11.0975, 17.1920};
        double worst = 0, worst_tab = 0;
        std::vector<double> E(z.size(), rate(T));
        for (int i = 0; i < 3; ++i) {
            const double t = ages[i] * GYR;
            double exact;
            if (t <= ta) exact = 0.5 * l1 * (std::sqrt(1.0 + 4.0 * muE * t / (m1 * l1)) - 1.0);
            else {
                const double b = 1.0 + 2.0 * a / l1, c = (t - ta) * muE / m2;
                exact = a + 0.5 * l2 * (-b + std::sqrt(b * b + 4.0 * c / l2));
            }
            const auto r = ice::solveIceRetentionColumn(z, E.data(), col, params(ages[i]));
            worst = std::fmax(worst, std::fabs(r.z_tau_m / exact - 1.0));
            worst_tab = std::fmax(worst_tab, std::fabs(r.z_tau_m * 100 - ref_cm[i]) / 0.5e-4);
            std::printf("    %.1f Gyr  z_tau = %.5f cm (exact %.5f, table %.4f)\n", ages[i],
                        r.z_tau_m * 100, exact * 100, ref_cm[i]);
        }
        gate("B1", std::fabs(ta / GYR - 0.1448) < 0.5e-4, "t(a) = %.5f Gyr (table 0.1448)", ta / GYR);
        gate("B2", worst < 1e-6, "two-layer z_tau vs closed form: worst relative %.2e (tol 1e-6)", worst);
        gate("B3", worst_tab <= 1.0, "against the reference values to their printed digits: worst %.2f half-units", worst_tab);
    }

    // ---- C: underflow ---------------------------------------------------------
    {
        const int F = 3, L = static_cast<int>(zprod.size());
        std::vector<float> sum(static_cast<size_t>(F) * L, 0.0f);
        const auto r = ice::computeIceRetention(sum.data(), F, L, 10, zprod, params(1.0), nullptr, false);
        bool ok = r.error.empty() && r.n_underflow == F;
        for (int f = 0; f < F && ok; ++f)
            ok = r.cls[f] == ice::kRetainedUnderflow && r.layer[f] == 0.0f && r.depth_m[f] >= 0.0f &&
                 r.depth_m[f] < 1e-9f && r.toplayer_gyr[f] == FLT_MAX && r.column_gyr[f] == FLT_MAX &&
                 r.exposed_loss_m[f] == 0.0f;
        gate("C1", ok, "<E> = 0 everywhere: class 3, layer 0, lifetimes FLT_MAX, no NaN (%.0f faces)", (double)r.n_underflow);
        std::vector<float> bad(static_cast<size_t>(L), 1.0f);
        bad[3] = std::nanf("");
        const auto rb = ice::computeIceRetention(bad.data(), 1, L, 1, zprod, params(1.0), nullptr, false);
        gate("C2", rb.error.empty() && rb.n_not_evaluated == 1 && rb.cls[0] == 255.0f && std::isnan(rb.depth_m[0]),
             "a non-finite <E>: not evaluated (class 255, depth NaN), never a class", 0);
    }

    // ---- D: grid convergence --------------------------------------------------
    {
        struct Case { double Tm, A, delta, ref_cm; } cases[] = {
            {112, 25, 0.03, 3.302}, {118, 40, 0.04, 8.016}, {108, 10, 0.03, 0.797}};
        const double zb = 0.4706301768531821, r23 = 1.14201;
        const double q = std::pow(r23, 0.25);
        const std::vector<double> g89 = geometric(88, zb, q);
        auto decimate = [](const std::vector<double>& g, int k) {
            std::vector<double> o; for (size_t i = 0; i < g.size(); i += k) o.push_back(g[i]); return o; };
        const std::vector<double> gref = geometric(88 * 32, zb, std::pow(q, 1.0 / 32));
        const std::vector<std::vector<double>> grids = {decimate(g89, 8), decimate(g89, 4), decimate(g89, 2), g89};
        ice::IceRetentionParams p = params(1.0);
        bool ok_rate = true, ok_23 = true, ok_ref = true;
        for (const auto& c : cases) {
            auto meanE = [&](double z) {             // time mean of E over one cycle (spectral)
                const int N = 2048; double s = 0;
                for (int k = 0; k < N; ++k)
                    s += rate(c.Tm + c.A * std::exp(-z / c.delta) * std::cos(2 * M_PI * (k + 0.5) / N));
                return s / N;
            };
            auto zt = [&](const std::vector<double>& g) {
                std::vector<double> E(g.size()); for (size_t i = 0; i < g.size(); ++i) E[i] = meanE(g[i]);
                const auto col = ice::makeIceColumn(p, g, p.density_scale_m);
                return ice::solveIceRetentionColumn(g, E.data(), col, p).z_tau_m;
            };
            const double zr = zt(gref);
            double err[4];
            for (int i = 0; i < 4; ++i) err[i] = zt(grids[i]) / zr - 1.0;
            std::printf("    T_m %3.0f K, A %2.0f K, delta %2.0f cm: ref %.4f cm (table %.3f); "
                        "error %+.3f%% %+.3f%% %+.3f%% %+.3f%% at 12/23/45/89 nodes (z_1 %.2f/%.2f/%.2f/%.2f mm)\n",
                        c.Tm, c.A, c.delta * 100, zr * 100, c.ref_cm, err[0] * 100, err[1] * 100,
                        err[2] * 100, err[3] * 100, grids[0][1] * 1e3, grids[1][1] * 1e3,
                        grids[2][1] * 1e3, grids[3][1] * 1e3);
            for (int i = 0; i < 3; ++i) ok_rate &= std::fabs(err[i]) >= 3.0 * std::fabs(err[i + 1]);
            ok_23 &= std::fabs(err[1]) < 0.005;
            ok_ref &= std::fabs(zr * 100 - c.ref_cm) <= 0.5e-3 + 1e-9;
        }
        gate("D1", ok_rate, "grid convergence: error falls by >= 3 per halving of the spacing", 0);
        gate("D2", ok_23, "error below 0.5%% on the 23-node production grid", 0);
        gate("D3", ok_ref, "reference depths match the specification's table to its printed digits", 0);
    }

    // ---- E: layer bookkeeping ------------------------------------------------
    {
        // Constant inputs; choose ages whose ice table lands just above and
        // just below the layer-0/1 boundary (z_1/2) and the layer-1/2 boundary.
        const double m = 400, ell = 75e-6, T = 112;
        const double muE = rate(T) * MU;
        auto age_for = [&](double z) { return m * (z + z * z / ell) / muE / GYR; };
        const double b0 = 0.5 * zprod[1], b1 = 0.5 * (zprod[1] + zprod[2]);
        struct { double z; int layer; } probes[] = {{0.99 * b0, 0}, {1.01 * b0, 1}, {0.99 * b1, 1}, {1.01 * b1, 2}};
        bool ok = true;
        for (const auto& pr : probes) {
            ice::IceRetentionFace fr;
            solve_const(T, m, ell, age_for(pr.z), zprod, &fr);
            ok &= fr.layer == pr.layer && fr.cls == (pr.layer == 0 ? ice::kRetainedTopLayer : ice::kRetainedAtDepth);
            ok &= std::fabs(fr.t_top_s / GYR / age_for(b0) - 1.0) < 1e-9;
        }
        gate("E1", ok, "layer index from control volumes; top-layer lifetime = t(z_1/2)", 0);
        // A weight fraction that does not fit in the pores is refused.
        ice::IceRetentionParams pw = params(1.0);
        pw.inventory = ice::IceInventory::WeightFraction; pw.weight_fraction = 0.6;
        std::string err;
        const auto col = ice::makeIceColumn(pw, zprod, 0.06, &err);
        gate("E2", col.zq.empty() && !err.empty(), "60 wt%% ice exceeds the pore space: refused", 0);
    }

    // ---- F: steep rate drop across the first interval ----------------------
    // <E> falls by 325, 1e4 and 8e4 between the surface node and z_1 (the
    // range met on cold-trap faces), then by 2 per node. Constant m and ell,
    // so the integrand of (4) is m (1 + 2 z/ell) / (mu <E(z)>) with log-linear
    // <E> and nothing else. Reference: composite Simpson on 2^14 panels per
    // node interval and bisection on it, written here independently of the
    // solver. A quadrature that is exact only for slowly varying <E> (the
    // trapezoid rule on 16 sub-steps is 1 to 4% high here) fails F1 and F2.
    {
        const double m = 500.0, ell = 75e-6;
        const std::vector<double>& z = zprod;
        const auto col = ice::makeIceColumnFromFunctions(z, [m](double) { return m; },
                                                         [ell](double) { return ell; });
        double worst_top = 0, worst_z = 0, worst_col = 0;
        for (double drop : {325.0, 1.0e4, 8.0e4}) {
            std::vector<double> E(z.size());
            E[0] = rate(112.0);
            E[1] = E[0] / drop;
            for (size_t j = 2; j < z.size(); ++j) E[j] = E[j - 1] / 2.0;
            auto g = [&](double x, int i) {   // integrand on node interval i
                const double s = (x - z[i]) / (z[i + 1] - z[i]);
                const double lE = std::log(E[i]) + s * (std::log(E[i + 1]) - std::log(E[i]));
                return m * (1.0 + 2.0 * x / ell) / (MU * std::exp(lE));
            };
            auto simpson = [&](int i, double a, double b) {
                const int N = 1 << 14;
                const double h = (b - a) / N;
                double s = g(a, i) + g(b, i);
                for (int k = 1; k < N; ++k) s += (k % 2 ? 4.0 : 2.0) * g(a + k * h, i);
                return s * h / 3.0;
            };
            const double t_top = simpson(0, 0.0, 0.5 * z[1]);
            double t_col = 0;
            for (size_t i = 0; i + 1 < z.size(); ++i) t_col += simpson(static_cast<int>(i), z[i], z[i + 1]);
            const double t_z1 = simpson(0, 0.0, z[1]);
            // Ages that put the ice table at 20%, 50% (inside layer 0) and 80%
            // of the steep interval.
            for (double frac : {0.2, 0.5, 0.8}) {
                const double zt_ref_target = frac * z[1];
                const double age_s = simpson(0, 0.0, zt_ref_target);
                double lo = 0.0, hi = z[1];
                for (int it = 0; it < 60; ++it) {
                    const double mid = 0.5 * (lo + hi);
                    (simpson(0, 0.0, mid) < age_s ? lo : hi) = mid;
                }
                const double zt_ref = 0.5 * (lo + hi);
                ice::IceRetentionParams pa = params(age_s / GYR);
                const auto fr = ice::solveIceRetentionColumn(z, E.data(), col, pa);
                worst_z = std::fmax(worst_z, std::fabs(fr.z_tau_m / zt_ref - 1.0));
                if (frac == 0.5) {
                    worst_top = std::fmax(worst_top, std::fabs(fr.t_top_s / t_top - 1.0));
                    worst_col = std::fmax(worst_col, std::fabs(fr.t_column_s / t_col - 1.0));
                    std::printf("    drop %6.0f: t(z_1/2) = %.6e Gyr (reference %.6e), t(z_1) %.4e Gyr\n",
                                drop, fr.t_top_s / GYR, t_top / GYR, t_z1 / GYR);
                }
            }
        }
        gate("F1", worst_top < 1e-9, "steep <E> drop (325-8e4 across z_1): top-layer lifetime vs reference, worst %.2e", worst_top);
        gate("F2", worst_z < 1e-9, "steep <E> drop: ice-table depth inside the steep interval, worst relative %.2e", worst_z);
        gate("F3", worst_col < 1e-9, "steep <E> drop: column lifetime vs reference, worst relative %.2e", worst_col);
    }

    // ---- G: the grain density is a property of the body ---------------------
    {
        auto resolve = [](const std::string& body, const std::string& mat, double body_rho,
                          double cfg, double& out, std::string& src) {
            out = -1; src.clear();
            return ice::resolveRetentionGrainDensity(body, mat, body_rho, cfg, out, src);
        };
        const auto moon = thermal::PlanetaryData::getPropertiesByName("Moon");
        const auto merc = thermal::PlanetaryData::getPropertiesByName("Mercury");
        const auto eur  = thermal::PlanetaryData::getPropertiesByName("Europa");
        double rg; std::string src;
        bool ok_moon = moon.regolith_material == "silicate" && moon.grain_density == 3100.0 &&
                       resolve("Moon", moon.regolith_material, moon.grain_density, -1, rg, src).empty() &&
                       rg == 3100.0 && src.find("Moon default") != std::string::npos;
        ok_moon &= resolve("Moon", moon.regolith_material, moon.grain_density, 2900, rg, src).empty() &&
                   rg == 2900.0 && src.find("set in") != std::string::npos;
        ok_moon &= merc.regolith_material == "silicate" &&
                   resolve("Mercury", merc.regolith_material, merc.grain_density, -1, rg, src).empty() &&
                   rg == 3100.0;
        gate("G1", ok_moon, "Moon, Mercury: the body's grain density (3100), or the configured one", 0);
        // A generic body (asteroid): no material, no grain density.
        const std::string e_ast = resolve("101955 Bennu", "", 0.0, -1, rg, src);
        const bool ok_ast = !e_ast.empty() && e_ast.find("grain_density_kgm3") != std::string::npos &&
                            resolve("101955 Bennu", "", 0.0, 3300, rg, src).empty() && rg == 3300.0;
        gate("G2", ok_ast, "generic body: refused until grain_density_kgm3 is set, then that value", 0);
        const std::string e_eur = resolve("Europa", eur.regolith_material, eur.grain_density, -1, rg, src);
        const std::string e_eur2 = resolve("Europa", eur.regolith_material, eur.grain_density, 3100, rg, src);
        gate("G3", eur.regolith_material == "water_ice" && !e_eur.empty() && !e_eur2.empty() &&
                   e_eur.find("criterion: flux") != std::string::npos,
             "Europa (ice regolith): refused, also with a configured grain density", 0);
        ice::IceRetentionParams pn = params(1.0);
        pn.grain_density_kgm3 = 0.0;
        gate("G4", ice::validateIceRetentionParams(pn).find("grain density is not set") != std::string::npos,
             "no grain density: the parameter check refuses (no silent default)", 0);
    }

    if (g_fail) { std::printf("[test_ice_retention] FAIL (%d)\n", g_fail); return 1; }
    std::printf("[test_ice_retention] PASS\n");
    return 0;
}
