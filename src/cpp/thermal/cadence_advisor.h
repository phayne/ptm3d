// cadence_advisor.h -- is the output step short enough to see a run's sunlit flashes?
//
// A log-and-provenance-only advisory. It never changes a setting or an output.
//
// Frame-sampled statistics (T_max, direct_solar_max) see the Sun only at the output epochs, so a face
// lit for less than one output step can be recorded as colder than it is. On rough terrain such
// flashes come from two kinds of event: the Sun crossing a gap between occluders, which lasts about
// (delta + theta)/omega, and the Sun's disk grazing an occluding edge, which can be arbitrarily short.
// Here omega is the Sun's apparent angular rate across the sky ((2 pi / P_syn) cos(subsolar latitude),
// measured from the run's own sun directions), theta the solar disk's angular diameter, and delta the
// gap's angular width. Because grazing events have no lower limit on duration, no output step catches
// every flash, and the angular size of the mesh's occluders does NOT predict how short its flashes are
// (Fermi, n4990 cadence ladder, 2026-11-04: Spearman -0.20 between the occluder angle and the true event
// duration). The advisor therefore recommends a fixed fraction of the disk-crossing time,
//     dt_out <= 0.4 theta / omega          (about 25 min on the Moon),
// and gives the 5% peak-flux limit 0.3 theta/omega.
//
// The occluder probe (IlluminationModel::probeBlockerAngles, one shadow ray per face per sampled sun
// position) is used only to decide whether flashes can occur at all ("not needed" when no face has an
// occluder between it and the Sun) and to identify at-risk faces for PTM-818.
//
// Indicators, from the frames the run already has:
//   N1 = area lit in exactly one output frame. For dt >= theta/omega it was never below the area whose
//        sunlit events the frames miss entirely, over all phases of the n4990 ladder (ratio min 1.00,
//        median 1.7-3.2). Below theta/omega it can under-count (ratio down to 0.34), because every full
//        transit is then caught at least twice; there the advisor uses N1 + N2 (lit in one or two
//        frames; ratio min 0.94).
//   R_ct(thr) = (that area, over ALL faces) / (area with frame-sampled T_max < thr): the estimated upper
//        bound on the fractional overstatement of the cold-trap area at threshold thr. n4990 (one phase),
//        110 K: 4.74 / 1.24 / 1.26% at 5.9 / 2 / 1 h (the 1 h step is below theta/omega: one or two frames;
//        one frame alone gives 0.70%), against true overstatements of 2.5 / 0.71-0.87 / 0.38-0.54%.
//   dT_miss = min[(2 F_abs / Gamma) sqrt(dt / pi), F_abs / (4 eps sigma T^3)], F_abs = S (1 - A) / r_AU^2:
//        the most a missed full-sun flash of duration dt can raise the surface temperature (semi-infinite
//        conductor, normal incidence), capped by the linearized radiative balance.
//
// TODO(PTM-818): an optional direct count of missed events -- re-trace only the at-risk faces (an
// occluder, lit in <= 2 output frames) at a step of theta/(4 omega) over the output window and compare
// with the frames. Not built: on the 509k Moon3D mesh that is ~2,700 sun positions per lunation times
// ~10^5 faces, minutes of ray tracing, so it needs its own switch and cost estimate.
//
// Header-only; the pure functions here are what test_cadence_advisor exercises.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace ptm {
namespace cadence {

constexpr double kRecC = 0.4;                 // recommended step = kRecC theta/omega
constexpr double kPeak5W = 0.3;               // dt <= kPeak5W theta/omega captures the peak flux to 5%
constexpr double kLitThresholdWm2 = 1.0;      // a face is "lit" in a frame when direct solar > 1 W/m^2
constexpr int    kProbeSunSamples = 64;       // sun positions sampled by the occluder probe
constexpr double kSolarConstant1AU = 1361.0;  // W/m^2
constexpr double kDeltaNonCore = M_PI;        // probe value for a block by a non-core geometry (size not measured)
constexpr double kVerdictTolerance = 1.02;    // a step within 2% of the recommendation is OK (25.0 vs 24.98 min)

struct SunGeometry {
    bool   ok = false;
    double omega = 0.0;      // rad/s, median over output steps
    double theta = 0.0;      // rad, median solar-disk angular diameter
    double dt_out = 0.0;     // s, median output step
    double r_au = 0.0;       // median Sun distance [AU]
    bool   omega_from_day = false;   // fewer than 2 output steps: omega = 2 pi / day (no declination)
};

// sun_km: body-frame sun positions [km] at the run's epochs; times_s: the epochs [s];
// [i0, i1): the output window. sun_radius_km: the solar radius. day_s: fallback period.
inline SunGeometry sunGeometry(const std::vector<std::array<double, 3>>& sun_km,
                               const std::vector<double>& times_s, size_t i0, size_t i1,
                               double sun_radius_km, double day_s) {
    SunGeometry g;
    i1 = std::min({i1, sun_km.size(), times_s.size()});
    if (i0 >= i1) return g;
    std::vector<double> rate, disk, step, dist;
    for (size_t t = i0; t < i1; ++t) {
        const auto& a = sun_km[t];
        const double ra = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if (ra > sun_radius_km) { disk.push_back(2.0 * std::asin(sun_radius_km / ra)); dist.push_back(ra); }
        if (t + 1 < i1) {
            const auto& b = sun_km[t + 1];
            const double rb = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
            const double dt = times_s[t + 1] - times_s[t];
            if (ra > 0 && rb > 0 && dt > 0) {
                const double cx = a[1] * b[2] - a[2] * b[1], cy = a[2] * b[0] - a[0] * b[2],
                             cz = a[0] * b[1] - a[1] * b[0];
                const double ang = std::atan2(std::sqrt(cx * cx + cy * cy + cz * cz),
                                              a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
                rate.push_back(ang / dt);
                step.push_back(dt);
            }
        }
    }
    auto median = [](std::vector<double> v) {
        if (v.empty()) return 0.0;
        std::nth_element(v.begin(), v.begin() + static_cast<long>(v.size() / 2), v.end());
        return v[v.size() / 2];
    };
    g.theta = median(disk);
    g.r_au = median(dist) / 1.495978707e8;
    if (!rate.empty()) {
        g.omega = median(rate);
        g.dt_out = median(step);
    } else if (day_s > 0) {
        g.omega = 2.0 * M_PI / day_s;
        g.omega_from_day = true;
    }
    g.ok = g.omega > 0 && g.theta > 0;
    return g;
}

// Frames-lit count that the indicator uses: 1 for dt >= theta/omega, 1..2 below it.
inline int indicatorMaxLit(const SunGeometry& g) {
    return (g.ok && g.dt_out > 0 && g.dt_out < g.theta / g.omega) ? 2 : 1;
}

// Area fraction (of the area in mask, or of everything) with 1 <= n_lit <= max_lit.
inline double litFrameFraction(const std::vector<int>& n_lit, const std::vector<double>& area, int max_lit,
                               const std::vector<char>* mask = nullptr, double* denom_out = nullptr) {
    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < n_lit.size() && i < area.size(); ++i) {
        if (mask && !(*mask)[i]) continue;
        den += area[i];
        if (n_lit[i] >= 1 && n_lit[i] <= max_lit) num += area[i];
    }
    if (denom_out) *denom_out = den;
    return den > 0 ? num / den : 0.0;
}

// R_ct(thr): indicator area over ALL faces / area with T_max < thr. Returns -1 when there is no such area.
inline double coldTrapOverstatementBound(const std::vector<int>& n_lit, const std::vector<double>& area,
                                         const std::vector<float>& tmax, double thr, int max_lit,
                                         double* ct_area_frac = nullptr) {
    double num = 0.0, ct = 0.0, tot = 0.0;
    for (size_t i = 0; i < n_lit.size() && i < area.size() && i < tmax.size(); ++i) {
        tot += area[i];
        if (n_lit[i] >= 1 && n_lit[i] <= max_lit) num += area[i];
        if (tmax[i] < thr) ct += area[i];
    }
    if (ct_area_frac) *ct_area_frac = tot > 0 ? ct / tot : 0.0;
    return ct > 0 ? num / ct : -1.0;
}

// The most a missed full-sun flash of duration dt can raise the surface [K]: semi-infinite conductor
// (2 F / Gamma) sqrt(dt / pi), capped by the linearized radiative balance F / (4 eps sigma T^3).
inline double missedFlashDeltaT(double F_abs, double Gamma, double eps, double T, double dt) {
    if (!(F_abs > 0) || !(Gamma > 0) || !(dt > 0)) return 0.0;
    const double cond = 2.0 * F_abs / Gamma * std::sqrt(dt / M_PI);
    const double rad = F_abs / (4.0 * eps * 5.670374419e-8 * T * T * T);
    return std::min(cond, rad);
}

struct CadenceReport {
    SunGeometry sun;
    bool   probed = false;            // occluder probe ran (CPU scene available)
    double frac_resolved = 0.0;       // area fraction with an occluder at some sampled sun position
    double t_rec = 0.0;               // s, recommended output step = 0.4 theta/omega
    double t_peak5 = 0.0;             // s, 5% peak-flux limit
    int    max_lit = 1;               // indicator: frames lit in [1, max_lit]
    double indicator_frac = -1.0;     // mesh-area fraction with that frame count (-1 = not counted)
    double tau110 = 0.0, tau54 = 0.0; // s, surface thermal response time
    double F_abs = 0.0;               // W/m^2, full-sun absorbed flux at normal incidence
    double dT_run = 0.0, dT_rec = 0.0;  // K, dT_miss at the run's step and at the recommended step (110 K)
};

// delta: per-face probe result (rad; +inf = no occluder; empty = probe did not run).
inline CadenceReport evaluateCadence(const SunGeometry& sun, const std::vector<float>& delta,
                                     const std::vector<double>& area) {
    CadenceReport r;
    r.sun = sun;
    if (!sun.ok) return r;
    r.t_rec = kRecC * sun.theta / sun.omega;
    r.t_peak5 = kPeak5W * sun.theta / sun.omega;
    r.max_lit = indicatorMaxLit(sun);
    r.probed = !delta.empty() && delta.size() == area.size();
    if (r.probed) {
        double tot = 0.0, res = 0.0;
        for (size_t i = 0; i < area.size(); ++i) {
            tot += area[i];
            if (std::isfinite(delta[i])) res += area[i];
        }
        r.frac_resolved = tot > 0 ? res / tot : 0.0;
    }
    return r;
}

inline std::string fmtDuration(double s) {
    char b[48];
    if (!std::isfinite(s)) std::snprintf(b, sizeof b, "n/a");
    else if (s < 60.0) std::snprintf(b, sizeof b, "%.3g s", s);
    else if (s < 3600.0) std::snprintf(b, sizeof b, "%.3g min", s / 60.0);
    else if (s < 2.0 * 86400.0) std::snprintf(b, sizeof b, "%.3g h", s / 3600.0);
    else std::snprintf(b, sizeof b, "%.3g d", s / 86400.0);
    return b;
}

inline std::string fmtKelvin(double k) {
    char b[48];
    if (k > 100.0) std::snprintf(b, sizeof b, "more than 100 K");
    else std::snprintf(b, sizeof b, "%.2g K", k);
    return b;
}

inline std::string verdictName(const CadenceReport& r) {
    if (!r.sun.ok) return "UNMEASURED";
    if (r.probed && r.frac_resolved == 0.0) return "not needed";
    return (r.sun.dt_out > 0 && r.sun.dt_out <= kVerdictTolerance * r.t_rec) ? "OK" : "COARSE";
}

inline std::string indicatorPhrase(int max_lit) {
    return max_lit == 1 ? "exactly one output frame" : "one or two output frames";
}

inline std::vector<std::string> formatCadenceAdvice(const CadenceReport& r, size_t n_probe_samples) {
    const std::string P = "[cadence-advisor] ";
    std::vector<std::string> L;
    char b[1000];
    const double deg = 180.0 / M_PI;
    if (!r.sun.ok) {
        L.push_back(P + "WARNING: could not measure the Sun's angular rate or disk from this run's epochs, so "
                        "the output-cadence check did NOT run.");
        return L;
    }
    const std::string v = verdictName(r);
    if (v == "not needed")
        std::snprintf(b, sizeof b, "not needed: no face on this mesh has an occluder between it and the Sun at "
                                   "the sampled sun positions, so there are no short sunlit flashes to miss.");
    else if (v == "OK")
        std::snprintf(b, sizeof b, "OK: the output step (%s) is no longer than the recommended %.1f theta/omega = "
                                   "%s. Very short grazing flashes are still missed, a small residual (see the "
                                   "post-run estimate).",
                      fmtDuration(r.sun.dt_out).c_str(), kRecC, fmtDuration(r.t_rec).c_str());
    else
        std::snprintf(b, sizeof b, "COARSE: the output step (%s) is longer than the recommended %.1f theta/omega = "
                                   "%s. Frame-sampled T_max can miss sunlit flashes (the Sun crossing gaps between "
                                   "occluders, and grazing partial-disk events, which can be arbitrarily short), so "
                                   "cold-trap areas can be overstated; the post-run line estimates by how much.",
                      fmtDuration(r.sun.dt_out).c_str(), kRecC, fmtDuration(r.t_rec).c_str());
    L.push_back(P + b);
    std::snprintf(b, sizeof b, "  Sun: angular rate omega = %.3g deg/h (%s), disk theta = %.3g deg, theta/omega "
                               "= %s.",
                  r.sun.omega * deg * 3600.0,
                  r.sun.omega_from_day ? "2 pi / solar day: fewer than 2 output steps, declination not included"
                                       : "from the run's own sun directions",
                  r.sun.theta * deg, fmtDuration(r.sun.theta / r.sun.omega).c_str());
    L.push_back(P + b);
    if (r.probed)
        std::snprintf(b, sizeof b, "  Occluder probe: %.3g%% of the area has an occluder between it and the Sun at "
                                   "some of the %zu sampled sun positions. Used only to decide whether flashes can "
                                   "occur; the occluder's angular size does not predict how short they are.",
                      100.0 * r.frac_resolved, n_probe_samples);
    else
        std::snprintf(b, sizeof b, "  Occluder probe: did not run (it needs the CPU ray tracer, which this run did "
                                   "not build); the verdict assumes flashes can occur.");
    L.push_back(P + b);
    std::snprintf(b, sizeof b, "  Peak flux of a flash is captured to within 5%% only for an output step <= "
                               "%.1f theta/omega = %s%s.",
                  kPeak5W, fmtDuration(r.t_peak5).c_str(),
                  (r.sun.dt_out > r.t_peak5) ? "; at this step, frame-sampled peak flux and T_max of sunlit "
                                               "flashes are underestimates" : "");
    L.push_back(P + b);
    if (r.tau110 > 0)
        L.push_back(P + "  Surface thermal response time tau = (Gamma/(4 eps sigma T^3))^2: " + fmtDuration(r.tau110) +
                    " at 110 K, " + fmtDuration(r.tau54) + " at 54 K.");
    if (r.F_abs > 0) {
        std::snprintf(b, sizeof b, "  A missed full-sun flash lasting one output step can raise a 110 K surface by up "
                                   "to dT = min[(2F/Gamma) sqrt(dt/pi), F/(4 eps sigma T^3)] = %s at this step and %s "
                                   "at the recommended step (F = S(1-A)/r^2 = %.3g W/m^2, normal incidence).",
                      fmtKelvin(r.dT_run).c_str(), fmtKelvin(r.dT_rec).c_str(), r.F_abs);
        L.push_back(P + b);
    }
    if (r.indicator_frac >= 0) {
        if (r.max_lit == 1)
            std::snprintf(b, sizeof b, "  Faces lit (direct solar > %g W/m^2) in exactly one output frame: %.3g%% of "
                                       "the mesh area. For steps >= theta/omega this was never below the mesh area "
                                       "whose sunlit events the frames miss entirely, over all phases of the n4990 "
                                       "cadence ladder (ratio min 1.00, median 1.7-3.2).",
                          kLitThresholdWm2, 100.0 * r.indicator_frac);
        else
            std::snprintf(b, sizeof b, "  Faces lit (direct solar > %g W/m^2) in one or two output frames: %.3g%% of "
                                       "the mesh area. This step is shorter than theta/omega, where the one-frame "
                                       "count can under-count the missed area (ratio down to 0.34 on the n4990 "
                                       "cadence ladder); the one-or-two-frame count's minimum ratio there was 0.94.",
                          kLitThresholdWm2, 100.0 * r.indicator_frac);
        L.push_back(P + b);
    }
    L.push_back(P + "  Silence with illumination.cadence_advisor: false. Criterion and guidance: "
                    "docs/thermal-model.md, \"Choosing the output cadence\".");
    return L;
}

// The post-run line. r110 / r54: R_ct (or -1 when there is no area below the threshold); ct110 / ct54:
// the cold-trap area fractions of the mesh. has_tmax false: no per-face T_max for this run.
inline std::string formatPostRun(int max_lit, double indicator_frac, bool has_tmax, double r110, double ct110,
                                 double r54, double ct54) {
    char b[1200];
    auto pct = [](double x) { char c[32]; std::snprintf(c, sizeof c, "%.3g%%", 100.0 * x); return std::string(c); };
    if (!has_tmax) {
        std::snprintf(b, sizeof b, "[cadence-advisor] post-run: faces lit in %s: %s of the mesh area (no per-face "
                                   "T_max for this run, so no cold-trap estimate).",
                      indicatorPhrase(max_lit).c_str(), pct(indicator_frac).c_str());
        return b;
    }
    auto est = [&](double r) { return r < 0 ? std::string("n/a (no area below the threshold)") : pct(r); };
    std::snprintf(b, sizeof b, "[cadence-advisor] post-run: estimated upper bound on the fractional overstatement of "
                               "the cold-trap area = (area of all faces lit in %s, %s of the mesh) / (area with "
                               "frame-sampled T_max below the threshold): %s at 110 K (cold traps %s of the mesh), "
                               "%s at 54 K (%s). On the n4990 cadence ladder (one phase) this estimate was "
                               "4.74 / 1.24 / 1.26%% at 5.9 h / 2 h / 1 h steps (one frame, one frame, one or two "
                               "frames), against true overstatements of 2.5 / 0.71-0.87 / 0.38-0.54%%.",
                  indicatorPhrase(max_lit).c_str(), pct(indicator_frac).c_str(), est(r110).c_str(), pct(ct110).c_str(),
                  est(r54).c_str(), pct(ct54).c_str());
    return b;
}

}  // namespace cadence
}  // namespace ptm
