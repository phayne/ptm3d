// test_cadence_advisor -- the output-cadence advisor (thermal/cadence_advisor.h). Header-only under test.
//
//   [1] sunGeometry recovers omega = (2 pi / P) cos(dec), theta = 2 asin(R_sun / d) and the distance from a
//       synthetic sun track at three declinations; the Moon's theta/omega ~ 1.05 h (Fermi 2026-11-04).
//   [2] evaluateCadence: the recommended step is 0.4 theta/omega whatever the occluder angles; the 5%
//       peak-flux limit 0.3 theta/omega; "not needed" only when no face is occluded; the indicator uses
//       1 frame for dt >= theta/omega and 1-2 frames below it.
//   [3] litFrameFraction and coldTrapOverstatementBound (numerator over ALL faces).
//   [4] missedFlashDeltaT: the conductive form, the radiative cap, and Fermi's Europa estimate (~16 K at
//       35.5 min, ~6 K at 6.7 min).
//   [5] formatting: verdicts, the silence key, the prefix, the UNMEASURED warning, the post-run line.
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "thermal/cadence_advisor.h"

namespace cd = ptm::cadence;
static int g_fail = 0, g_pass = 0;
static void check(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    (ok ? g_pass : g_fail)++;
}
static bool near(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); }

static void track(double P_s, double dec, double d_km, double dt, int n, std::vector<std::array<double, 3>>& s,
                  std::vector<double>& t) {
    s.clear(); t.clear();
    for (int k = 0; k < n; ++k) {
        const double ph = 2.0 * M_PI * (k * dt) / P_s;
        s.push_back({d_km * std::cos(dec) * std::cos(ph), d_km * std::cos(dec) * std::sin(ph), d_km * std::sin(dec)});
        t.push_back(k * dt);
    }
}

int main() {
    const double Rsun = 695700.0, AU = 1.495978707e8, P = 29.530588 * 86400.0;
    std::printf("[1] sun geometry from a synthetic track\n");
    std::vector<std::array<double, 3>> s; std::vector<double> t;
    for (double dec_deg : {0.0, 1.54, 25.0}) {
        track(P, dec_deg * M_PI / 180.0, AU, 7200.0, 360, s, t);
        const auto g = cd::sunGeometry(s, t, 0, s.size(), Rsun, P);
        const double om = 2.0 * M_PI / P * std::cos(dec_deg * M_PI / 180.0);
        char b[160];
        std::snprintf(b, sizeof b, "dec %.2f deg: omega %.6g vs (2pi/P)cos(dec) %.6g; theta %.5f deg; r %.6f AU",
                      dec_deg, g.omega, om, g.theta * 180 / M_PI, g.r_au);
        check(g.ok && near(g.omega, om, 1e-4) && near(g.theta, 2.0 * std::asin(Rsun / AU), 1e-9) &&
                  near(g.dt_out, 7200.0, 1e-12) && near(g.r_au, 1.0, 1e-12), b);
    }
    track(P, 0.0, AU, 7200.0, 10, s, t);
    {
        const auto g = cd::sunGeometry(s, t, 0, s.size(), Rsun, P);
        const double h = g.theta / g.omega / 3600.0;
        check(std::fabs(h - 1.05) < 0.01, "Moon theta/omega = " + std::to_string(h) + " h (Fermi: 1.05 h)");
        const auto g1 = cd::sunGeometry(s, t, 0, 1, Rsun, P);
        check(g1.ok && g1.omega_from_day && near(g1.omega, 2.0 * M_PI / P, 1e-12),
              "one output step: omega falls back to 2 pi / day, and says so");
    }

    std::printf("[2] recommendation and verdict\n");
    track(P, 0.0, AU, 7200.0, 10, s, t);
    const auto g = cd::sunGeometry(s, t, 0, s.size(), Rsun, P);   // 2 h step, theta/omega 1.05 h
    {
        std::vector<double> a(100, 1.0);
        std::vector<float> small(100, 0.0f), big(100, static_cast<float>(10.0 * M_PI / 180.0));
        const auto r1 = cd::evaluateCadence(g, small, a), r2 = cd::evaluateCadence(g, big, a);
        check(near(r1.t_rec, 0.4 * g.theta / g.omega, 1e-12) && r1.t_rec == r2.t_rec,
              "recommended step = 0.4 theta/omega, independent of the occluder angles");
        check(near(r1.t_peak5, 0.3 * g.theta / g.omega, 1e-12), "5% peak-flux limit = 0.3 theta/omega");
        check(cd::verdictName(r1) == "COARSE" && cd::verdictName(r2) == "COARSE", "a 2 h step is COARSE on the Moon");
        check(r1.max_lit == 1, "dt (2 h) >= theta/omega (1.05 h): the indicator counts one frame");
        const std::vector<float> none(5, std::numeric_limits<float>::infinity());
        const auto rn = cd::evaluateCadence(g, none, std::vector<double>(5, 1.0));
        check(cd::verdictName(rn) == "not needed" && rn.frac_resolved == 0.0, "no occluded face: not needed");
        const auto rp = cd::evaluateCadence(g, {}, std::vector<double>(5, 1.0));
        check(!rp.probed && cd::verdictName(rp) == "COARSE", "no probe: flashes assumed possible");
        const std::vector<float> ring(5, static_cast<float>(cd::kDeltaNonCore));
        check(cd::evaluateCadence(g, ring, std::vector<double>(5, 1.0)).frac_resolved == 1.0,
              "a block by a non-core geometry (pi) counts as occluded");
    }
    {
        track(P, 0.0, AU, 1500.0, 10, s, t);   // 25 min step
        const auto g25 = cd::sunGeometry(s, t, 0, s.size(), Rsun, P);
        const auto r = cd::evaluateCadence(g25, std::vector<float>(3, 0.01f), std::vector<double>(3, 1.0));
        check(cd::verdictName(r) == "OK" && r.max_lit == 2, "25 min step: OK, and the indicator counts 1-2 frames");
        cd::SunGeometry g2 = g25; g2.dt_out = 1.015 * r.t_rec;
        cd::SunGeometry g3 = g25; g3.dt_out = 1.03 * r.t_rec;
        check(cd::verdictName(cd::evaluateCadence(g2, {}, {})) == "OK" &&
                  cd::verdictName(cd::evaluateCadence(g3, {}, {})) == "COARSE",
              "verdict tolerance: 1.5% over the recommendation is OK, 3% is COARSE");
    }

    std::printf("[3] indicator fractions\n");
    {
        const std::vector<int> n = {0, 1, 2, 3, 1};
        const std::vector<double> a = {1, 2, 3, 4, 10};
        check(near(cd::litFrameFraction(n, a, 1), 12.0 / 20.0, 1e-15), "one frame: 12/20");
        check(near(cd::litFrameFraction(n, a, 2), 15.0 / 20.0, 1e-15), "one or two frames: 15/20");
        const std::vector<float> tmax = {40, 300, 100, 50, 300};   // cold (<110 K): faces 0, 2, 3 = area 8
        double ct = 0;
        const double r = cd::coldTrapOverstatementBound(n, a, tmax, 110.0, 1, &ct);
        check(near(r, 12.0 / 8.0, 1e-15) && near(ct, 8.0 / 20.0, 1e-15),
              "R_ct: numerator over ALL faces (12) / cold-trap area (8)");
        check(cd::coldTrapOverstatementBound(n, a, tmax, 10.0, 1) < 0, "no area below the threshold: -1");
    }

    std::printf("[4] missed-flash temperature rise\n");
    {
        const double F = 1361.0 * (1.0 - 0.6) / (5.2 * 5.2);   // Europa, A = 0.6
        const double d35 = cd::missedFlashDeltaT(F, 65.0, 0.9, 110.0, 35.5 * 60.0);
        const double d67 = cd::missedFlashDeltaT(F, 65.0, 0.9, 110.0, 6.72 * 60.0);
        check(std::fabs(d35 - 16.0) < 1.5 && std::fabs(d67 - 7.0) < 1.5,
              "Europa: " + std::to_string(d35) + " K at 35.5 min, " + std::to_string(d67) +
                  " K at 6.7 min (Fermi: ~16 K, ~6 K)");
        const double cap = F / (4.0 * 0.9 * 5.670374419e-8 * 110.0 * 110.0 * 110.0);
        check(near(cd::missedFlashDeltaT(F, 65.0, 0.9, 110.0, 1e9), cap, 1e-12), "long step: the radiative cap");
    }

    std::printf("[5] formatting\n");
    {
        auto r = cd::evaluateCadence(g, std::vector<float>(10, 0.01f), std::vector<double>(10, 1.0));
        r.indicator_frac = 0.0123; r.tau110 = 3340; r.tau54 = 107000; r.F_abs = 1198; r.dT_run = 4000; r.dT_rec = 12;
        const auto L = cd::formatCadenceAdvice(r, 64);
        for (const auto& l : L) std::printf("      %s\n", l.c_str());
        bool prefix = true, silence = false, over = false;
        for (const auto& l : L) {
            prefix = prefix && l.rfind("[cadence-advisor] ", 0) == 0;
            silence = silence || l.find("illumination.cadence_advisor: false") != std::string::npos;
            over = over || l.find("more than 100 K") != std::string::npos;
        }
        check(prefix && silence && over && L[0].find("COARSE: the output step (2 h)") != std::string::npos,
              "prefix, verdict, '> 100 K' wording, silence key");
        const auto Lb = cd::formatCadenceAdvice(cd::evaluateCadence(cd::SunGeometry{}, {}, {}), 0);
        check(Lb.size() == 1 && Lb[0].find("did NOT run") != std::string::npos, "unmeasured sun: says so");
        const std::string pr = cd::formatPostRun(2, 0.01, true, 0.0124, 0.35, -1, 0.0);
        std::printf("      %s\n", pr.c_str());
        check(pr.find("one or two output frames") != std::string::npos && pr.find("1.24% at 110 K") != std::string::npos &&
                  pr.find("n/a (no area below the threshold) at 54 K") != std::string::npos,
              "post-run line: frame count, R_ct, n/a");
    }
    std::printf("\n%d passed, %d failed\n%s\n", g_pass, g_fail, g_fail ? "FAIL" : "PASS: all cadence-advisor checks");
    return g_fail ? 1 : 0;
}
