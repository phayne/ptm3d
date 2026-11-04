// test_eq_seasonal -- the two-timescale (seasonal) seed correction, eqSeedSeasonalApply().
//
// No arguments: self-contained checks.
//   S1  constant season (every phase identical to the base series): the correction
//       is identically zero, so the seed comes back unchanged (|dT| < 1e-6 K).
//   S2  a warm season BEFORE the seed epoch raises the deep column and lowers
//       nothing at depth below the annual mean; refusal paths report a reason.
// With a file argument: parity harness. The file (little-endian float64) holds
//   nz, z[nz], kc[nz], rho[nz], K, NS, period_s, cycle_s, tau_s, flux[K*NS],
//   n_eq, span_eq_s, Feq[n_eq]
// and the program prints the base seed and the corrected seed, one node per line,
// so an independent implementation can be compared on IDENTICAL inputs.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "../thermal/eq_seed.h"
#include "../thermal/thermal_model.h"

namespace {
int g_fail = 0;
void check(bool ok, const std::string& what, const std::string& detail = "") {
    std::printf("  [%s] %s%s%s\n", ok ? "PASS" : "FAIL", what.c_str(), detail.empty() ? "" : " -- ", detail.c_str());
    if (!ok) ++g_fail;
}
const double kCp[4] = {2.7431, 2.3616e-3, -1.2340e-5, 8.9093e-9};   // Moon cp(T) coefficients, cp0 = -3.6125

struct Run {
    std::vector<double> z, kc, rho, T0, T1;
    thermal::EqSeedResult base;
    thermal::EqSeedSeasonalResult sr;
};

Run run(const std::vector<double>& z, const std::vector<double>& kc, const std::vector<double>& rho,
        const std::vector<double>& Feq, double span_eq, int K, int NS, double period, double cycle,
        double tau, const std::vector<float>& flux) {
    Run r; r.z = z; r.kc = kc; r.rho = rho;
    thermal::EqSeedForcing fc;
    fc.flux = Feq.data(); fc.n_samples = static_cast<int>(Feq.size()); fc.span_s = span_eq; fc.cycle_s = cycle;
    fc.cp0 = -3.6125; fc.cp_coeff = kCp; fc.n_cp_coeff = 4;
    thermal::EqSeedInput in;
    double m = 0.0; for (double v : Feq) m += v; m /= Feq.size();
    in.mean_absorbed_flux = m; in.Qb = 0.018; in.emissivity = 0.95; in.sigma = thermal::SIGMA_SB;
    in.z = r.z.data(); in.kc = r.kc.data(); in.rho = r.rho.data(); in.n = static_cast<int>(z.size());
    in.R350 = thermal::R350(2.7); in.forcing = &fc;
    r.T0.resize(z.size());
    r.base = thermal::eqSeedColumn(in, r.T0.data());
    r.T1 = r.T0;
    thermal::EqSeedSeasonal s;
    s.K = K; s.NS = NS; s.period_s = period; s.cycle_s = cycle; s.tau_seed_s = tau; s.flux = flux.data();
    r.sr = thermal::eqSeedSeasonalApply(in, r.base, s, r.T1.data());
    return r;
}

void moonColumn(std::vector<double>& z, std::vector<double>& kc, std::vector<double>& rho, int n = 40, double zmax = 0.4706) {
    z.resize(n); kc.resize(n); rho.resize(n);
    const double rr = std::pow(zmax / 1e-3, 1.0 / (n - 2));
    double dz = 1e-3; z[0] = 0.0;
    for (int i = 1; i < n; ++i) { z[i] = z[i - 1] + dz; dz *= rr; }
    const double sc = zmax / z[n - 1];
    for (int i = 0; i < n; ++i) {
        z[i] *= sc; const double df = std::exp(-z[i] / 0.06);
        kc[i] = 3.4e-3 + (7.4e-4 - 3.4e-3) * df; rho[i] = 1800.0 + (1100.0 - 1800.0) * df;
    }
}
}  // namespace

int main(int argc, char** argv) {
    const double PS = 29.530588 * 86400.0, PA = 346.62 * 86400.0;
    if (argc > 1) {
        std::ifstream f(argv[1], std::ios::binary);
        auto rd = [&]() { double v = 0.0; f.read(reinterpret_cast<char*>(&v), 8); return v; };
        const int nz = static_cast<int>(rd());
        std::vector<double> z(nz), kc(nz), rho(nz);
        for (auto& v : z) v = rd(); for (auto& v : kc) v = rd(); for (auto& v : rho) v = rd();
        const int K = static_cast<int>(rd()), NS = static_cast<int>(rd());
        const double period = rd(), cycle = rd(), tau = rd();
        std::vector<float> flux(static_cast<std::size_t>(K) * NS);
        for (auto& v : flux) v = static_cast<float>(rd());
        const int neq = static_cast<int>(rd()); const double span = rd();
        std::vector<double> Feq(neq); for (auto& v : Feq) v = rd();
        if (!f) { std::fprintf(stderr, "FATAL: short parity file %s\n", argv[1]); return 2; }
        const Run r = run(z, kc, rho, Feq, span, K, NS, period, cycle, tau, flux);
        std::printf("applied %d why '%s' g_mean %.9f g_base %.9f anchor %s\n", r.sr.applied ? 1 : 0, r.sr.why,
                    r.sr.g_mean, r.sr.g_base, thermal::eqSeedAnchorStatusName(r.base.anchor_status));
        for (int j = 0; j < nz; ++j) std::printf("%.9f %.9f %.9f\n", z[j], r.T0[j], r.T1[j]);
        return r.sr.applied ? 0 : 1;
    }
    std::printf("=== two-timescale seasonal seed correction ===\n");
    std::vector<double> z, kc, rho; moonColumn(z, kc, rho);
    const int K = 12, NS = 24;
    std::vector<double> day(NS);
    for (int i = 0; i < NS; ++i) { const double ph = 2.0 * M_PI * (i + 0.5) / NS; day[i] = std::max(0.0, 300.0 * std::sin(ph)); }
    {   // S1: constant season, base on the very same series -> zero correction
        std::vector<float> flux(static_cast<std::size_t>(K) * NS);
        for (int k = 0; k < K; ++k) for (int i = 0; i < NS; ++i) flux[k * NS + i] = static_cast<float>(day[i]);
        std::vector<double> Feq(day.begin(), day.end());
        for (auto& v : Feq) v = static_cast<double>(static_cast<float>(v));
        const Run r = run(z, kc, rho, Feq, PS, K, NS, PA, PS, 0.3 * PA, flux);
        double dmax = 0.0; for (std::size_t j = 0; j < z.size(); ++j) dmax = std::max(dmax, std::fabs(r.T1[j] - r.T0[j]));
        check(r.sr.applied, "S1 applied");
        char b[96]; std::snprintf(b, sizeof b, "max |dT| = %.3e K", dmax);
        check(dmax < 1e-6, "S1 constant season leaves the seed unchanged", b);
    }
    {   // S2: bright season centered half a year before the seed epoch, dark around it
        std::vector<float> flux(static_cast<std::size_t>(K) * NS);
        for (int k = 0; k < K; ++k) {
            const double sc = 0.5 + 0.5 * std::cos(2.0 * M_PI * (k + 0.5) / K);   // max at phase ~0
            for (int i = 0; i < NS; ++i) flux[k * NS + i] = static_cast<float>(sc * day[i]);
        }
        std::vector<double> Feq(NS); for (int i = 0; i < NS; ++i) Feq[i] = 0.0 * day[i] + 1e-3;   // dark eq window
        const double tau = 0.5 * PA;                                                       // seed in the dark season
        const Run r = run(z, kc, rho, Feq, PS, K, NS, PA, PS, tau, flux);
        check(r.sr.applied, "S2 applied", r.sr.why);
        char b[128]; std::snprintf(b, sizeof b, "deep %.2f -> %.2f K", r.T0.back(), r.T1.back());
        check(r.T1.back() > r.T0.back() + 5.0, "S2 dark-season seed after a bright season: deep column warmer than window-periodic", b);
        thermal::EqSeedSeasonal s0; s0.K = 2;
        thermal::EqSeedInput bad; std::vector<double> T = r.T0;
        const auto rr = thermal::eqSeedSeasonalApply(bad, r.base, s0, T.data());
        check(!rr.applied && std::string(rr.why).size() > 0, "S2 refusal reports a reason", rr.why);
    }
    {   // S3: state + apply-state == one-shot apply; the min-correction gate leaves the column untouched
        std::vector<float> flux(static_cast<std::size_t>(K) * NS);
        for (int k = 0; k < K; ++k) {
            const double sc = 0.5 + 0.5 * std::cos(2.0 * M_PI * (k + 0.5) / K);
            for (int i = 0; i < NS; ++i) flux[k * NS + i] = static_cast<float>(sc * day[i]);
        }
        std::vector<double> Feq(day.begin(), day.end());
        const Run r = run(z, kc, rho, Feq, PS, K, NS, PA, PS, 0.2 * PA, flux);
        thermal::EqSeedForcing fc; fc.flux = Feq.data(); fc.n_samples = NS; fc.span_s = PS; fc.cycle_s = PS;
        fc.cp0 = -3.6125; fc.cp_coeff = kCp; fc.n_cp_coeff = 4;
        thermal::EqSeedInput in; in.mean_absorbed_flux = 0.0; in.Qb = 0.018; in.emissivity = 0.95; in.sigma = thermal::SIGMA_SB;
        in.z = z.data(); in.kc = kc.data(); in.rho = rho.data(); in.n = static_cast<int>(z.size()); in.R350 = thermal::R350(2.7); in.forcing = &fc;
        thermal::EqSeedSeasonal s; s.K = K; s.NS = NS; s.period_s = PA; s.cycle_s = PS; s.tau_seed_s = 0.2 * PA; s.flux = flux.data();
        std::vector<double> H(z.size()), T = r.T0;
        const auto st = thermal::eqSeedSeasonalState(in, s, H.data());
        const auto ap = thermal::eqSeedSeasonalApplyState(in, r.base, H.data(), 0.0, T.data());
        double dmax = 0.0; for (std::size_t j = 0; j < z.size(); ++j) dmax = std::max(dmax, std::fabs(T[j] - r.T1[j]));
        char b[96]; std::snprintf(b, sizeof b, "max |state path - one-shot| = %.3e K", dmax);
        check(st.applied && ap.applied && dmax < 1e-9, "S3 state + apply-state reproduces the one-shot correction", b);
        std::vector<double> T2 = r.T0;
        const auto gate = thermal::eqSeedSeasonalApplyState(in, r.base, H.data(), 1e6, T2.data());
        double d2 = 0.0; for (std::size_t j = 0; j < z.size(); ++j) d2 = std::max(d2, std::fabs(T2[j] - r.T0[j]));
        check(!gate.applied && d2 == 0.0 && std::string(gate.why) == "below min correction", "S3 min-correction gate leaves the column untouched", gate.why);
    }
    std::printf(g_fail ? "=== %d FAILURE(S) ===\n" : "=== OK (0 failures) ===\n", g_fail);
    return g_fail ? 1 : 0;
}
