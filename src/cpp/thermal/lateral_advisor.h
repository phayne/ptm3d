// lateral_advisor.h -- when does lateral heat conduction matter for THIS mesh?
//
// A pre-run advisory, log only. It never enables, disables or alters anything: it compares the
// mesh's facet sizes with the critical length below which lateral conduction is not negligible,
// and says whether thermal.lateral_conduction is recommended.
//
// The criterion is the one in the Moon3D paper (Hayne et al., in revision, Section 2; after
// Hayne et al. 2021, Methods): a feature of size D is conductively isolated from its neighbors
// when its own local conductance to radiative and conductive loss,
//     h = 4 eps sigma T^3 + Gamma sqrt(omega),   Gamma = sqrt(k rho c_p),
// outweighs the lateral conductance across it, i.e. when D exceeds
//     D* = 8 k / (pi h).
// D* is small on warm surfaces, where radiation dominates h (about 0.9 mm on sunlit lunar
// regolith), and grows as the surface cools (about 6.5 mm at 110 K, 2 cm at 70 K, 4 cm at 50 K).
// The paper's numbers use the SURFACE-layer material (k_s, rho_s, the Hayne 2017 k(T) and c_p(T)),
// which is what the engine's own D* (ptm::lateral::conductionLength, used by the lateral
// edge-activation rule) evaluates; test_lateral_advisor pins this header's D* to that function.
//
// Facet size is the MEAN EDGE LENGTH, the same measure the engine's lateral activation rule uses
// (lateral_adjacency.cpp, Adjacency::mean_edge), so "smaller than D*" here means what it means
// there.
//
// Header-only and dependency-free apart from conductivity_model.h (itself header-only), so the
// unit test cannot be satisfied by a stale library.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "conductivity_model.h"

namespace ptm {
namespace lateral {

// Representative temperatures [K] at which the advisory evaluates D*.
constexpr double kAdvisorColdT   = 50.0;   // cold trap (the CO2 regime); D* ~ 4 cm on lunar regolith
constexpr double kAdvisorIceT    = 110.0;  // water-ice stability threshold; D* ~ 6.5 mm
constexpr double kAdvisorSunlitT = 300.0;  // warm, sunlit surface; D* ~ 0.9 mm

// Recommend lateral conduction when facets smaller than D*(kAdvisorColdT) cover more than this
// fraction of the surface area.
constexpr double kAdvisorRecommendAreaFraction = 0.01;
// Borderline: not recommended, but facets below D*(cold) cover more than kAdvisorBorderlineFloor
// of the area, or facets below kAdvisorBorderlineFactor * D*(cold) cover more than
// kAdvisorRecommendAreaFraction. The factor is the engine's default edge-activation multiple
// (thermal.lateral_activation_c = 2): edges the selective-activation rule would couple.
constexpr double kAdvisorBorderlineFloor  = 0.001;
constexpr double kAdvisorBorderlineFactor = 2.0;
// "Not needed" still carries the rock caveat when the median facet is within this factor of
// D*(cold): rock (k ~ 0.1-1 W/m/K) has D* of decimeters to meters.
constexpr double kAdvisorRockCaveatFactor = 100.0;

// The surface node's material, exactly as the engine's batch solver sets it at z = 0
// (thermal_model.cpp initializeBatchGrid: kc = k_s, rho = rho_s at zero depth;
// updateMaterialPropertiesBatchSIMD: k = kc (1 + R350 T^3) or Martinez 2021,
// c_p = cp0 + T (c1 + c2 T + ...)).
struct SurfaceMaterial {
    double ks = 0.0;                    // surface contact conductivity k_s [W/m/K]
    double rhos = 0.0;                  // surface density rho_s [kg/m^3]
    double cp0 = 0.0;                   // c_p constant term [J/kg/K]
    std::vector<double> cpCoeff;        // c_p = cp0 + sum_j cpCoeff[j] T^(j+1)
    double emissivity = 0.0;
    double R350 = 0.0;                  // chi / 350^3
    thermal::ConductivityModel model = thermal::ConductivityModel::Hayne2017;
    double day_s = 0.0;                 // solar day [s]; omega = 2 pi / day_s
    double chi = 0.0;                   // for the log only (R350 is what the law uses)
};

inline double surfaceHeatCapacity(const SurfaceMaterial& m, double T) {
    if (m.cpCoeff.empty()) return m.cp0;
    double r = m.cpCoeff.back();
    for (int j = static_cast<int>(m.cpCoeff.size()) - 2; j >= 0; --j) r = r * T + m.cpCoeff[j];
    return m.cp0 + T * r;
}

inline double surfaceConductivity(const SurfaceMaterial& m, double T) {
    return thermal::conductivityK(m.model, T, m.ks, m.rhos, m.R350);
}

// h = 4 eps sigma T^3 + Gamma sqrt(omega) [W/m^2/K].
inline double surfaceConductance(const SurfaceMaterial& m, double T) {
    const double sigma = 5.670374419e-8;   // the literal ptm::lateral::conductionLength uses
    const double omega = 2.0 * M_PI / m.day_s;
    const double k = surfaceConductivity(m, T);
    return 4.0 * m.emissivity * sigma * T * T * T +
           std::sqrt(k * m.rhos * surfaceHeatCapacity(m, T)) * std::sqrt(omega);
}

// D*(T) = 8 k / (pi h) [m]. Same expression tree as ptm::lateral::conductionLength.
inline double criticalLength(const SurfaceMaterial& m, double T) {
    const double k = surfaceConductivity(m, T);
    return 8.0 * k / (M_PI * surfaceConductance(m, T));
}

// Per-face facet size (mean edge length) and area, both in meters.
struct FacetSizes {
    std::vector<double> size_m;
    std::vector<double> area_m2;
};

// V: indexable vertices with [0..2]; F: indexable faces with integer [0..2].
// to_m: meters per mesh length unit (1 for a mesh in meters, 1000 for kilometers).
// Degenerate or out-of-range faces are skipped (counted by the caller as size_m.size() < n).
template <class VertexVec, class FaceVec>
FacetSizes measureFacetSizes(const VertexVec& vertices, const FaceVec& faces, double to_m) {
    FacetSizes out;
    out.size_m.reserve(faces.size());
    out.area_m2.reserve(faces.size());
    const long long nv = static_cast<long long>(vertices.size());
    for (const auto& f : faces) {
        const long long i0 = static_cast<long long>(f[0]);
        const long long i1 = static_cast<long long>(f[1]);
        const long long i2 = static_cast<long long>(f[2]);
        if (i0 < 0 || i1 < 0 || i2 < 0 || i0 >= nv || i1 >= nv || i2 >= nv) continue;
        const auto& a = vertices[static_cast<size_t>(i0)];
        const auto& b = vertices[static_cast<size_t>(i1)];
        const auto& c = vertices[static_cast<size_t>(i2)];
        const double ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
        const double vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
        const double wx = c[0] - b[0], wy = c[1] - b[1], wz = c[2] - b[2];
        const double e = (std::sqrt(ux * ux + uy * uy + uz * uz) + std::sqrt(vx * vx + vy * vy + vz * vz) +
                          std::sqrt(wx * wx + wy * wy + wz * wz)) / 3.0;
        const double cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx;
        const double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
        if (!(std::isfinite(e) && std::isfinite(area) && area > 0.0)) continue;
        out.size_m.push_back(e * to_m);
        out.area_m2.push_back(area * to_m * to_m);
    }
    return out;
}

enum class AdvisorVerdict { Recommended, Borderline, NotNeeded, Unmeasured };

inline const char* verdictName(AdvisorVerdict v) {
    switch (v) {
        case AdvisorVerdict::Recommended: return "RECOMMENDED";
        case AdvisorVerdict::Borderline:  return "BORDERLINE";
        case AdvisorVerdict::NotNeeded:   return "not needed";
        default:                          return "UNMEASURED";
    }
}

struct AdvisorReport {
    AdvisorVerdict verdict = AdvisorVerdict::Unmeasured;
    double dstar_cold = 0.0, dstar_ice = 0.0, dstar_sunlit = 0.0;   // [m]
    double frac_cold = 0.0, frac_ice = 0.0, frac_sunlit = 0.0;      // area fraction with size < D*(T)
    double frac_cold_x2 = 0.0;                                      // ... < kAdvisorBorderlineFactor D*(cold)
    double size_median = 0.0, size_p05 = 0.0, size_p95 = 0.0;       // [m], by face count
    size_t n_faces = 0;                                             // faces measured
};

namespace detail {
// Nearest-rank percentile of an unsorted copy (q in [0, 1]).
inline double percentile(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    const size_t k = static_cast<size_t>(q * static_cast<double>(v.size() - 1));
    std::nth_element(v.begin(), v.begin() + static_cast<long>(k), v.end());
    return v[k];
}
}  // namespace detail

inline AdvisorReport evaluateLateralAdvice(const FacetSizes& fs, const SurfaceMaterial& m) {
    AdvisorReport r;
    r.dstar_cold   = criticalLength(m, kAdvisorColdT);
    r.dstar_ice    = criticalLength(m, kAdvisorIceT);
    r.dstar_sunlit = criticalLength(m, kAdvisorSunlitT);
    r.n_faces = fs.size_m.size();
    double total = 0.0, a_cold = 0.0, a_ice = 0.0, a_sun = 0.0, a_x2 = 0.0;
    for (size_t i = 0; i < fs.size_m.size(); ++i) {
        const double l = fs.size_m[i], a = fs.area_m2[i];
        total += a;
        if (l < r.dstar_cold) a_cold += a;
        if (l < r.dstar_ice) a_ice += a;
        if (l < r.dstar_sunlit) a_sun += a;
        if (l < kAdvisorBorderlineFactor * r.dstar_cold) a_x2 += a;
    }
    const bool finite_dstar = std::isfinite(r.dstar_cold) && r.dstar_cold > 0.0 &&
                              std::isfinite(r.dstar_ice) && std::isfinite(r.dstar_sunlit);
    if (r.n_faces == 0 || !(total > 0.0) || !finite_dstar) return r;   // Unmeasured
    r.frac_cold = a_cold / total;
    r.frac_ice = a_ice / total;
    r.frac_sunlit = a_sun / total;
    r.frac_cold_x2 = a_x2 / total;
    r.size_median = detail::percentile(fs.size_m, 0.50);
    r.size_p05 = detail::percentile(fs.size_m, 0.05);
    r.size_p95 = detail::percentile(fs.size_m, 0.95);
    if (r.frac_cold > kAdvisorRecommendAreaFraction)
        r.verdict = AdvisorVerdict::Recommended;
    else if (r.frac_cold > kAdvisorBorderlineFloor || r.frac_cold_x2 > kAdvisorRecommendAreaFraction)
        r.verdict = AdvisorVerdict::Borderline;
    else
        r.verdict = AdvisorVerdict::NotNeeded;
    return r;
}

// What the formatter needs to know about the run besides the mesh.
struct AdvisorContext {
    bool lateral_on = false;
    std::string unit_label = "meters";        // how the mesh lengths were read
    std::vector<std::string> notes;           // caller-supplied NOTE lines (unit caveats, maps)
    std::string material_label;               // e.g. "Moon surface layer"
};

inline std::string formatLength(double m) {
    char b[64];
    const double a = std::fabs(m);
    if (!std::isfinite(m)) std::snprintf(b, sizeof b, "%g", m);
    else if (a < 0.01) std::snprintf(b, sizeof b, "%.3g mm", m * 1e3);
    else if (a < 1.0) std::snprintf(b, sizeof b, "%.3g cm", m * 1e2);
    else if (a < 1000.0) std::snprintf(b, sizeof b, "%.3g m", m);
    else std::snprintf(b, sizeof b, "%.3g km", m / 1e3);
    return b;
}

// A size ratio without exponent notation: 2 significant figures below 10, else 3.
inline std::string formatRatio(double x) {
    char b[48];
    if (!(x > 0.0) || !std::isfinite(x)) std::snprintf(b, sizeof b, "%g", x);
    else if (x < 10.0) std::snprintf(b, sizeof b, "%.2g", x);
    else {
        const double scale = std::pow(10.0, std::floor(std::log10(x)) - 2.0);
        std::snprintf(b, sizeof b, "%.0f", std::round(x / scale) * scale);
    }
    return b;
}

inline std::string formatPercent(double f) {
    char b[32];
    const double p = 100.0 * f;
    if (p == 0.0) std::snprintf(b, sizeof b, "0%%");
    else if (p < 0.1) std::snprintf(b, sizeof b, "<0.1%%");
    else if (p < 10.0) std::snprintf(b, sizeof b, "%.1f%%", p);
    else std::snprintf(b, sizeof b, "%.0f%%", p);
    return b;
}

// The advisory as log lines, each already prefixed with "[lateral-advisor] ".
inline std::vector<std::string> formatLateralAdvice(const AdvisorReport& r, const SurfaceMaterial& m,
                                                    const AdvisorContext& ctx) {
    const std::string P = "[lateral-advisor] ";
    std::vector<std::string> L;
    char b[768];
    const std::string dc = formatLength(r.dstar_cold);

    if (r.verdict == AdvisorVerdict::Unmeasured) {
        L.push_back(P + "WARNING: could not measure the facet sizes of this mesh (or D* is not finite for "
                        "this material), so the lateral-conduction check did NOT run.");
        return L;
    }

    // ---- headline -----------------------------------------------------------------------------
    if (r.verdict == AdvisorVerdict::Recommended) {
        if (ctx.lateral_on)
            std::snprintf(b, sizeof b, "thermal.lateral_conduction is ON, and appropriate for this mesh: "
                                       "%s of the surface area is in facets smaller than D*(%g K) = %s.",
                          formatPercent(r.frac_cold).c_str(), kAdvisorColdT, dc.c_str());
        else
            std::snprintf(b, sizeof b,
                          "RECOMMENDED: %s of the surface area is in facets smaller than D*(%g K) = %s, so this mesh resolves "
                          "features small enough to exchange significant heat laterally. Without thermal.lateral_conduction, "
                          "small cold traps on this mesh are likely too cold, and cold-trap areas can be overestimated, most "
                          "strongly at the coldest thresholds.",
                          formatPercent(r.frac_cold).c_str(), kAdvisorColdT, dc.c_str());
    } else if (r.verdict == AdvisorVerdict::Borderline) {
        std::snprintf(b, sizeof b, "%s%s of the surface area is in facets smaller than D*(%g K) = %s, and "
                                   "%s in facets smaller than %g D*. "
                                   "Lateral conduction is likely a small correction for water-ice (110 K) cold traps on this mesh, but "
                                   "can still remove the coldest (CO2-class) cold traps, which on lunar regolith sit at 25-45 K where "
                                   "D* is 2-3 times D*(50 K); enable it if those matter.",
                      ctx.lateral_on ? "thermal.lateral_conduction is ON (borderline for this mesh): "
                                     : "BORDERLINE: ",
                      formatPercent(r.frac_cold).c_str(), kAdvisorColdT, dc.c_str(),
                      formatPercent(r.frac_cold_x2).c_str(), kAdvisorBorderlineFactor);
    } else {
        const double ratio = r.size_median / r.dstar_cold;
        if (ctx.lateral_on)
            std::snprintf(b, sizeof b, "thermal.lateral_conduction is ON but unnecessary for this mesh: the "
                                       "median facet (%s) is %s times D*(%g K) = %s, so facets are "
                                       "conductively isolated. It adds cost for a negligible effect.",
                          formatLength(r.size_median).c_str(), formatRatio(ratio).c_str(), kAdvisorColdT, dc.c_str());
        else
            std::snprintf(b, sizeof b, "not needed: the median facet (%s) is %s times D*(%g K) = %s, so "
                                       "facets are conductively isolated and exchange heat by radiation.",
                          formatLength(r.size_median).c_str(), formatRatio(ratio).c_str(), kAdvisorColdT, dc.c_str());
    }
    L.push_back(P + b);

    // ---- criterion and inputs -----------------------------------------------------------------
    char klaw[96];
    if (m.model == thermal::ConductivityModel::Martinez2021)
        std::snprintf(klaw, sizeof klaw, "k(T) Martinez & Siegler 2021");
    else
        std::snprintf(klaw, sizeof klaw, "k(T) Hayne et al. 2017 with chi = %g", m.chi);
    const std::string label = ctx.material_label.empty() ? "" : " (" + ctx.material_label + ")";
    std::snprintf(b, sizeof b, "  D* = 8k/(pi h), h = 4 eps sigma T^3 + Gamma sqrt(omega), from this run's "
                               "surface material%s: %s, k_s = %.3g W/m/K, rho_s = %.4g kg/m^3, eps = %.3g, "
                               "solar day %.4g d",
                  label.c_str(), klaw, m.ks, m.rhos, m.emissivity, m.day_s / 86400.0);
    L.push_back(P + b);
    std::snprintf(b, sizeof b, "  D* = %s at %g K (sunlit), %s at %g K (water-ice threshold), %s at %g K "
                               "(cold trap)",
                  formatLength(r.dstar_sunlit).c_str(), kAdvisorSunlitT, formatLength(r.dstar_ice).c_str(),
                  kAdvisorIceT, dc.c_str(), kAdvisorColdT);
    L.push_back(P + b);
    std::snprintf(b, sizeof b, "  facet size (mean edge, %zu faces, mesh read as %s): median %s, "
                               "5th-95th percentile %s - %s; area in facets below D*: %s (%g K), "
                               "%s (%g K), %s (%g K)",
                  r.n_faces, ctx.unit_label.c_str(), formatLength(r.size_median).c_str(),
                  formatLength(r.size_p05).c_str(), formatLength(r.size_p95).c_str(),
                  formatPercent(r.frac_sunlit).c_str(), kAdvisorSunlitT, formatPercent(r.frac_ice).c_str(),
                  kAdvisorIceT, formatPercent(r.frac_cold).c_str(), kAdvisorColdT);
    L.push_back(P + b);

    // ---- what to expect, and what it costs ------------------------------------------------------
    if (r.verdict != AdvisorVerdict::NotNeeded) {
        L.push_back(P + "  " +
                    "For reference (Moon3D crater, facets ~2 mm; Hayne et al., in revision): at 110 K, cold traps "
                    "smaller than 1 cm kept about half their area with lateral conduction, those of 1-5 cm about three "
                    "quarters, and the dominant region (equivalent diameter about 91 cm) 99%; facets on its margin "
                    "cooled into the class, leaving the total 110 K area nearly unchanged (+0.5%). At the 54 K CO2 "
                    "threshold every cold trap smaller than 5 cm vanished and the total area fell by 69%.");
        if (!ctx.lateral_on)
            L.push_back(P + "  " +
                        "To enable: thermal.lateral_conduction: true, lateral_scheme: monolithic, lateral_picard: 2, "
                        "lateral_activation_c: -1 (the Moon3D configuration), use_gpu: false. Lateral conduction runs on the "
                        "CPU SIMD solver only and solves a coupled system across the mesh at every substep, so expect a "
                        "substantially longer thermal stage.");
    } else if (r.size_median < kAdvisorRockCaveatFactor * r.dstar_cold) {
        L.push_back(P + "  NOTE: this assumes the configured regolith-like material. In rock (k ~ 0.1-1 "
                        "W/m/K), D* grows to decimeters or meters, comparable to these facets.");
    }
    for (const auto& n : ctx.notes) L.push_back(P + "  NOTE: " + n);
    L.push_back(P + "  Silence with thermal.lateral_advisor: false. Criterion and guidance: "
                    "docs/thermal-model.md, \"When to enable lateral conduction\".");
    return L;
}

}  // namespace lateral
}  // namespace ptm
