// ice_retention.cpp - Ice retention criterion (see ice_retention.h).

#include "ice_retention.h"
#include "ice_table.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ice {

namespace {
constexpr double MU = kWaterMoleculeMassKg;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

bool finitePositive(double v) { return std::isfinite(v) && v > 0.0; }

// Seconds -> Gyr as float32; a lifetime beyond the float range is written as
// FLT_MAX and counted (it is "longer than any age of interest", not a number).
float lifetimeGyr(double t_s, long long& n_beyond) {
    const double g = t_s / kSecondsPerGyr;
    if (!(g <= static_cast<double>(FLT_MAX))) { ++n_beyond; return FLT_MAX; }
    return static_cast<float>(g);
}
} // namespace

const char* iceInventoryName(IceInventory inv) {
    return inv == IceInventory::PoreFilling ? "pore_filling" : "weight_fraction";
}

bool parseIceInventory(const std::string& s, IceInventory& out) {
    if (s == "pore_filling")    { out = IceInventory::PoreFilling;    return true; }
    if (s == "weight_fraction") { out = IceInventory::WeightFraction; return true; }
    return false;
}

std::string validateIceRetentionParams(const IceRetentionParams& p) {
    std::ostringstream e;
    if (!finitePositive(p.age_gyr))
        e << "retention_age_gyr must be positive and finite (got " << p.age_gyr << "). ";
    if (!finitePositive(p.hop_length_m))
        e << "grain_hop_length_um must be positive and finite (got " << p.hop_length_m * 1e6 << "). ";
    if (!finitePositive(p.grain_density_kgm3))
        e << "grain_density_kgm3 must be positive and finite (got " << p.grain_density_kgm3 << "). ";
    if (!finitePositive(p.ice_density_kgm3))
        e << "ice density must be positive and finite. ";
    if (!finitePositive(p.rho_surface_kgm3) || !finitePositive(p.rho_deep_kgm3))
        e << "the regolith density profile (rho_surface " << p.rho_surface_kgm3 << ", rho_deep "
          << p.rho_deep_kgm3 << " kg m^-3) is not set. ";
    if (!finitePositive(p.density_scale_m))
        e << "the density e-folding scale H (" << p.density_scale_m << " m) is not set. ";
    if (finitePositive(p.grain_density_kgm3) &&
        (p.rho_surface_kgm3 >= p.grain_density_kgm3 || p.rho_deep_kgm3 >= p.grain_density_kgm3))
        e << "grain_density_kgm3 " << p.grain_density_kgm3 << " is not above the regolith bulk "
             "density (" << p.rho_surface_kgm3 << " to " << p.rho_deep_kgm3
          << " kg m^-3): the porosity 1 - rho/rho_g would not be positive. ";
    if (p.inventory == IceInventory::WeightFraction &&
        !(std::isfinite(p.weight_fraction) && p.weight_fraction > 0.0 && p.weight_fraction < 1.0))
        e << "ice_weight_fraction must be in (0, 1) (got " << p.weight_fraction << "). ";
    // The density profile runs monotonically from rho_surface to rho_deep, and
    // the ice a weight fraction needs, relative to the pore space, grows with
    // rho; so the two end values bound every depth of every column, whatever
    // its H.
    if (e.str().empty() && p.inventory == IceInventory::WeightFraction) {
        for (double rho : {p.rho_surface_kgm3, p.rho_deep_kgm3}) {
            const double need = p.weight_fraction / (1.0 - p.weight_fraction) * rho;
            const double cap = p.ice_density_kgm3 * (1.0 - rho / p.grain_density_kgm3);
            if (need > cap * (1.0 + 1e-12)) {
                e << "ice_weight_fraction " << p.weight_fraction << " needs " << need
                  << " kg m^-3 of ice where the regolith density is " << rho
                  << " kg m^-3, more than its pore space holds (" << cap << " kg m^-3). ";
                break;
            }
        }
    }
    return e.str();
}

double regolithPorosity(const IceRetentionParams& p, double z, double H) {
    const double rho = p.rho_deep_kgm3 + (p.rho_surface_kgm3 - p.rho_deep_kgm3) * std::exp(-z / H);
    return 1.0 - rho / p.grain_density_kgm3;
}

std::vector<double> retentionQuadratureGrid(const std::vector<double>& z) {
    std::vector<double> zq;
    if (z.size() < 2) return zq;
    const int S = kRetentionSubsteps;
    zq.reserve((z.size() - 1) * S + 1);
    for (size_t i = 0; i + 1 < z.size(); ++i)
        for (int k = 0; k < S; ++k)
            zq.push_back(z[i] + (static_cast<double>(k) / S) * (z[i + 1] - z[i]));
    zq.push_back(z.back());
    return zq;
}

IceColumnProfile makeIceColumnFromFunctions(const std::vector<double>& z_nodes,
                                            const std::function<double(double)>& m_of_z,
                                            const std::function<double(double)>& ell_of_z) {
    IceColumnProfile c;
    c.zq = retentionQuadratureGrid(z_nodes);
    if (c.zq.empty()) return c;
    const size_t NQ = c.zq.size();
    c.m_cell.resize(NQ - 1);
    c.lambda.assign(NQ, 0.0);
    for (size_t k = 0; k + 1 < NQ; ++k) {
        const double zm = 0.5 * (c.zq[k] + c.zq[k + 1]);
        c.m_cell[k] = m_of_z(zm);
        c.lambda[k + 1] = c.lambda[k] + (c.zq[k + 1] - c.zq[k]) / ell_of_z(zm);
    }
    return c;
}

IceColumnProfile makeIceColumn(const IceRetentionParams& p,
                               const std::vector<double>& z_nodes,
                               double H_face, std::string* err) {
    // Initial inventory m(z), from the porosity of the thermal model's own
    // density profile:
    //   pore filling     m = rho_ice phi(z),         phi = 1 - rho(z)/rho_g
    //   weight fraction  m = w / (1 - w) rho(z)      (must fit in the pores)
    // Transport length: constant ell (the overburden upgrade replaces this).
    const double H = H_face;
    auto rho = [&](double z) {
        return p.rho_deep_kgm3 + (p.rho_surface_kgm3 - p.rho_deep_kgm3) * std::exp(-z / H);
    };
    auto phi = [&](double z) { return 1.0 - rho(z) / p.grain_density_kgm3; };
    std::function<double(double)> m_of_z;
    if (p.inventory == IceInventory::PoreFilling)
        m_of_z = [&](double z) { return p.ice_density_kgm3 * phi(z); };
    else
        m_of_z = [&](double z) { return p.weight_fraction / (1.0 - p.weight_fraction) * rho(z); };
    const double ell = p.hop_length_m;
    IceColumnProfile c = makeIceColumnFromFunctions(z_nodes, m_of_z, [ell](double) { return ell; });
    // Physical admissibility, at every sub-step: positive porosity, and a
    // weight fraction that fits in the pore space.
    for (size_t k = 0; k < c.m_cell.size(); ++k) {
        const double zm = 0.5 * (c.zq[k] + c.zq[k + 1]);
        const double ph = phi(zm);
        const double cap = p.ice_density_kgm3 * ph;
        if (!(ph > 0.0 && ph < 1.0) || !(c.m_cell[k] > 0.0) || c.m_cell[k] > cap * (1.0 + 1e-12)) {
            if (err) {
                std::ostringstream e;
                e << "the ice inventory is not physical at z = " << zm << " m (H = " << H
                  << " m): porosity " << ph << ", ice " << c.m_cell[k]
                  << " kg m^-3 against a pore capacity of " << cap << " kg m^-3";
                if (p.inventory == IceInventory::WeightFraction)
                    e << " (ice_weight_fraction " << p.weight_fraction
                      << " needs more ice than the pores hold)";
                *err = e.str();
            }
            return IceColumnProfile{};
        }
    }
    return c;
}

namespace {
// phi1(u) = int_0^1 e^{-u s} ds = (1 - e^{-u}) / u and
// phi2(u) = int_0^1 s e^{-u s} ds = (1 - (1 + u) e^{-u}) / u^2, with their
// Taylor series where the closed forms cancel (|u| < 0.1: ten terms, error
// below 1e-16 relative).
inline void expIntegralWeights(double u, double& phi1, double& phi2) {
    if (std::fabs(u) < 0.1) {
        double term = 1.0;            // (-u)^n / n!
        phi1 = 0.0; phi2 = 0.0;
        for (int n = 0; n < 12; ++n) {
            phi1 += term / (n + 1);
            phi2 += term / (n + 2);
            term *= -u / (n + 1);
        }
        return;
    }
    const double em = -std::expm1(-u);    // 1 - e^{-u}
    phi1 = em / u;
    phi2 = (em - u * std::exp(-u)) / (u * u);
}
} // namespace

IceRetentionFace solveIceRetentionColumn(const std::vector<double>& z,
                                         const double* E,
                                         const IceColumnProfile& col,
                                         const IceRetentionParams& p) {
    IceRetentionFace r;
    const int n = static_cast<int>(z.size());
    const int S = kRetentionSubsteps;
    for (int j = 0; j < n; ++j) {
        if (!std::isfinite(E[j]) || E[j] < 0.0) {     // no value: never a class
            r.z_tau_m = kNaN; r.layer = -1; r.cls = kNotEvaluated;
            r.t_top_s = r.t_column_s = r.exposed_loss_m = kNaN;
            return r;
        }
    }
    const double tau_s = p.age_gyr * kSecondsPerGyr;
    r.exposed_loss_m = MU * E[0] * tau_s / p.ice_density_kgm3;

    // On sub-step k (x = z - zq[k] in [0, h]) the integrand of (4) is
    //   g(x) = (A + B x) exp(-(l_k + c x)),
    // with A = m_k (1 + 2 Lambda_k) / mu, B = 2 m_k (Lambda_{k+1} - Lambda_k) / (h mu)
    // (Lambda is linear on the sub-step: ell is constant there), l = ln <E>
    // linear between nodes with slope c. Its integral is exact:
    //   I(x) = exp(-l_k) [A x phi1(c x) + B x^2 phi2(c x)].
    // (The trapezoid rule, used here at first, overestimates the integral of
    // an exponential by about (c h)^2 / 12: 1 to 4% where <E> falls by orders
    // of magnitude between the first two nodes.)
    const size_t k_top = static_cast<size_t>(S / 2);   // zq index of the bottom of layer 0
    double t = 0.0;
    bool found = false;
    size_t k_found = 0;
    for (int i = 0; i + 1 < n; ++i) {
        const double l0 = std::log(std::max(E[i], kRetentionRateFloor));
        const double l1 = std::log(std::max(E[i + 1], kRetentionRateFloor));
        const double c  = (l1 - l0) / (z[i + 1] - z[i]);
        for (int s = 0; s < S; ++s) {
            const size_t k = static_cast<size_t>(i) * S + s;
            const double h  = col.zq[k + 1] - col.zq[k];
            const double lk = l0 + c * (col.zq[k] - z[i]);
            const double ea = std::exp(-lk);
            const double A  = col.m_cell[k] * (1.0 + 2.0 * col.lambda[k]) / MU;
            const double B  = col.m_cell[k] * 2.0 * (col.lambda[k + 1] - col.lambda[k]) / (h * MU);
            auto I = [&](double x) {
                double p1, p2;
                expIntegralWeights(c * x, p1, p2);
                return ea * (A * x * p1 + B * x * x * p2);
            };
            const double dt = I(h);
            if (!found && t + dt >= tau_s) {
                // Root of t + I(x) = tau on [0, h]. I increases strictly
                // (g > 0): Newton, kept inside a shrinking bracket.
                const double d = tau_s - t;
                double lo = 0.0, hi = h;
                const double g0 = ea * A;
                double x = (g0 > 0.0) ? std::min(h, d / g0) : 0.5 * h;
                for (int it = 0; it < 200; ++it) {
                    const double f = I(x) - d;
                    if (f > 0.0) hi = x; else lo = x;
                    if (f == 0.0 || hi - lo <= 1e-15 * h) break;
                    const double g = std::exp(-(lk + c * x)) * (A + B * x);
                    double xn = (g > 0.0 && std::isfinite(g)) ? x - f / g : 0.5 * (lo + hi);
                    if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
                    const bool done = std::fabs(xn - x) <= 1e-16 * h;
                    x = xn;
                    if (done) break;
                }
                x = std::min(std::max(x, 0.0), h);
                r.z_tau_m = col.zq[k] + x;
                found = true;
                k_found = k;
            }
            t += dt;
            if (k + 1 == k_top) r.t_top_s = t;
        }
    }
    r.t_column_s = t;
    if (!found) {
        r.z_tau_m = z.back();
        r.layer = -1;
        r.cls = kNoneWithinColumn;
        return r;
    }
    const int i = static_cast<int>(k_found / S);
    const double boundary = col.zq[static_cast<size_t>(i) * S + S / 2];   // (z_i + z_{i+1}) / 2
    r.layer = (r.z_tau_m < boundary) ? i : i + 1;
    if (r.layer == 0) r.cls = (E[0] == 0.0) ? kRetainedUnderflow : kRetainedTopLayer;
    else              r.cls = kRetainedAtDepth;
    return r;
}

namespace {

std::string checkGrid(const std::vector<double>& z, int n_layers) {
    if (static_cast<int>(z.size()) != n_layers || n_layers < 2)
        return "the depth grid has " + std::to_string(z.size()) + " nodes for " +
               std::to_string(n_layers) + " layers of sums";
    if (z[0] != 0.0) return "the depth grid does not start at the surface (z[0] != 0)";
    for (int j = 1; j < n_layers; ++j)
        if (!(z[j] > z[j - 1])) return "the depth grid is not strictly increasing";
    return "";
}

template <class T>
IceRetentionResult computeImpl(const T* sum_E, int F, int L, long long n_samples,
                               const std::vector<double>& z, const IceRetentionParams& p,
                               const std::vector<double>* Hpf, bool use_openmp) {
    IceRetentionResult r;
    std::string e = validateIceRetentionParams(p);
    if (e.empty()) e = checkGrid(z, L);
    if (e.empty() && (!sum_E || F <= 0)) e = "no sums";
    if (e.empty() && n_samples <= 0) e = "the accumulator holds no samples";
    if (e.empty() && Hpf && static_cast<int>(Hpf->size()) != F)
        e = "the per-face H map has " + std::to_string(Hpf->size()) + " entries for " +
            std::to_string(F) + " faces";
    IceColumnProfile shared;
    if (e.empty()) {
        shared = makeIceColumn(p, z, p.density_scale_m, &e);
        if (shared.zq.empty() && e.empty()) e = "empty column profile";
    }
    if (!e.empty()) { r.error = e; return r; }

    const size_t NF = static_cast<size_t>(F);
    r.depth_m.assign(NF, 0.f); r.layer.assign(NF, 0.f); r.cls.assign(NF, 0.f);
    r.toplayer_gyr.assign(NF, 0.f); r.column_gyr.assign(NF, 0.f); r.exposed_loss_m.assign(NF, 0.f);
    const double inv_n = 1.0 / static_cast<double>(n_samples);
    long long n_top = 0, n_depth = 0, n_none = 0, n_uf = 0, n_ne = 0, n_big = 0, n_deep_uf = 0;
    int profile_error = 0;
    std::string profile_error_text;

    auto body = [&](int f, std::vector<double>& E, IceColumnProfile& own, double& own_H,
                    long long& c_top, long long& c_depth, long long& c_none, long long& c_uf,
                    long long& c_ne, long long& c_big, long long& c_deep_uf) {
        for (int j = 0; j < L; ++j)
            E[j] = static_cast<double>(sum_E[static_cast<size_t>(j) * NF + f]) * inv_n;
        const IceColumnProfile* col = &shared;
        if (Hpf) {
            const double Hf = (*Hpf)[f] > 0.0 ? (*Hpf)[f] : p.density_scale_m;
            if (Hf != p.density_scale_m) {
                if (Hf != own_H) {
                    std::string pe;
                    own = makeIceColumn(p, z, Hf, &pe);
                    own_H = own.zq.empty() ? -1.0 : Hf;
                    if (own.zq.empty()) {
                        #ifdef _OPENMP
                        #pragma omp critical(ice_retention_profile_error)
                        #endif
                        { if (!profile_error++) profile_error_text = pe; }
                        return;
                    }
                }
                col = &own;
            }
        }
        const IceRetentionFace fr = solveIceRetentionColumn(z, E.data(), *col, p);
        r.depth_m[f] = static_cast<float>(fr.z_tau_m);
        r.layer[f]   = fr.cls == kNotEvaluated ? std::numeric_limits<float>::quiet_NaN()
                                               : static_cast<float>(fr.layer);
        r.cls[f]     = static_cast<float>(fr.cls);
        r.exposed_loss_m[f] = static_cast<float>(fr.exposed_loss_m);
        if (fr.cls == kNotEvaluated) {
            r.toplayer_gyr[f] = r.column_gyr[f] = std::numeric_limits<float>::quiet_NaN();
            ++c_ne;
            return;
        }
        r.toplayer_gyr[f] = lifetimeGyr(fr.t_top_s, c_big);
        r.column_gyr[f]   = lifetimeGyr(fr.t_column_s, c_big);
        switch (fr.cls) {
            case kRetainedTopLayer:  ++c_top;   break;
            case kRetainedAtDepth:   ++c_depth; break;
            case kNoneWithinColumn:  ++c_none;  break;
            case kRetainedUnderflow: ++c_uf;    break;
            default: break;
        }
        for (int j = 1; j < L; ++j) if (E[j] == 0.0) { ++c_deep_uf; break; }
    };

#ifdef _OPENMP
    if (use_openmp) {
        #pragma omp parallel reduction(+:n_top,n_depth,n_none,n_uf,n_ne,n_big,n_deep_uf)
        {
            std::vector<double> E(static_cast<size_t>(L));
            IceColumnProfile own; double own_H = -1.0;
            #pragma omp for schedule(static)
            for (int f = 0; f < F; ++f)
                body(f, E, own, own_H, n_top, n_depth, n_none, n_uf, n_ne, n_big, n_deep_uf);
        }
    } else
#endif
    {
        std::vector<double> E(static_cast<size_t>(L));
        IceColumnProfile own; double own_H = -1.0;
        for (int f = 0; f < F; ++f)
            body(f, E, own, own_H, n_top, n_depth, n_none, n_uf, n_ne, n_big, n_deep_uf);
    }
    (void)use_openmp;
    if (profile_error) {
        r = IceRetentionResult{};
        r.error = "per-face density profile on " + std::to_string(profile_error) +
                  " face(s): " + profile_error_text;
        return r;
    }
    r.n_top = n_top; r.n_depth = n_depth; r.n_none = n_none; r.n_underflow = n_uf;
    r.n_not_evaluated = n_ne; r.n_lifetime_beyond_float = n_big; r.n_deep_node_underflow = n_deep_uf;
    return r;
}

std::string fmt(double v, const char* f = "%.6g") {
    char b[64]; std::snprintf(b, sizeof b, f, v); return b;
}

} // namespace

IceRetentionResult computeIceRetention(const float* s, int F, int L, long long n,
                                       const std::vector<double>& z, const IceRetentionParams& p,
                                       const std::vector<double>* H, bool omp) {
    return computeImpl(s, F, L, n, z, p, H, omp);
}
IceRetentionResult computeIceRetention(const double* s, int F, int L, long long n,
                                       const std::vector<double>& z, const IceRetentionParams& p,
                                       const std::vector<double>* H, bool omp) {
    return computeImpl(s, F, L, n, z, p, H, omp);
}

std::string iceRetentionProvenance(const IceRetentionParams& p,
                                   const std::vector<double>& z,
                                   long long n_samples,
                                   const IceRetentionResult& r,
                                   const std::vector<double>* Hpf) {
    std::ostringstream o;
    const double H = p.density_scale_m;
    o << "criterion: ice retention after an age (" << kIceRetentionVersion << ")\n";
    o << "age_gyr: " << fmt(p.age_gyr) << "\n";
    if (p.inventory == IceInventory::PoreFilling) {
        o << "inventory: pore_filling (ice fills the pore space; grain density "
          << fmt(p.grain_density_kgm3) << " kg m^-3, ice density " << fmt(p.ice_density_kgm3)
          << " kg m^-3; porosity " << fmt(regolithPorosity(p, 0.0, H), "%.3f") << " at the surface to "
          << fmt(1.0 - p.rho_deep_kgm3 / p.grain_density_kgm3, "%.3f") << " at depth)\n";
    } else {
        o << "inventory: weight_fraction (ice_weight_fraction " << fmt(p.weight_fraction)
          << " of ice per mass of ice plus dry regolith; grain density "
          << fmt(p.grain_density_kgm3) << " kg m^-3 for the pore-capacity check)\n";
    }
    o << "density_profile: rho(z) = " << fmt(p.rho_deep_kgm3) << " + (" << fmt(p.rho_surface_kgm3)
      << " - " << fmt(p.rho_deep_kgm3) << ") exp(-z/H) kg m^-3, H = " << fmt(H) << " m";
    if (Hpf && !Hpf->empty()) {
        double lo = 1e300, hi = -1e300;
        for (double v : *Hpf) { const double h = v > 0.0 ? v : H; lo = std::min(lo, h); hi = std::max(hi, h); }
        o << " (per-face H map: " << fmt(lo) << " to " << fmt(hi) << " m)";
    }
    o << "\n";
    o << "transport_length: constant ell = " << fmt(p.hop_length_m * 1e6) << " um (D_K = v ell / 8)\n";
    if (z.size() >= 2)
        o << "grid: " << z.size() << " nodes, z_1 = " << fmt(z[1] * 1e3) << " mm, z_bot = "
          << fmt(z.back()) << " m; top layer 0 to " << fmt(0.5 * z[1] * 1e3) << " mm\n";
    o << "n_samples: " << n_samples << "\n";
    o << "classes: " << r.n_top << " top layer, " << r.n_depth << " at depth, " << r.n_none
      << " none within the column (censored), " << r.n_underflow
      << " top layer with surface rate below float range, " << r.n_not_evaluated << " not evaluated\n";
    if (r.n_lifetime_beyond_float)
        o << "lifetimes_beyond_float: " << r.n_lifetime_beyond_float
          << " lifetime values exceed the float32 range and are written as 3.4028235e+38\n";
    return o.str();
}

std::uint32_t crc32(const void* data, std::size_t n, std::uint32_t crc) {
    static std::uint32_t table[256];
    static const bool init = [] {
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1u) ? 0xEDB88320u ^ (c >> 1) : (c >> 1);
            table[i] = c;
        }
        return true;
    }();
    (void)init;
    const auto* p = static_cast<const unsigned char*>(data);
    crc = ~crc;
    for (std::size_t i = 0; i < n; ++i) crc = table[(crc ^ p[i]) & 0xFFu] ^ (crc >> 8);
    return ~crc;
}

bool writeIceRetentionProducts(const std::string& dir, const IceRetentionResult& r,
                               const std::string& prov, std::string* err) {
    std::vector<std::pair<std::string, std::uint32_t>> sums;
    auto put = [&](const char* name, const std::vector<float>& v) {
        std::ofstream os(dir + "/" + name, std::ios::binary);
        const std::uint32_t N = static_cast<std::uint32_t>(v.size());
        os.write(reinterpret_cast<const char*>(&N), sizeof N);
        os.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(float));
        if (!os.good()) { if (err) *err = std::string("write failed: ") + dir + "/" + name; return false; }
        // The checksum of the bytes on disk: the header, then the values.
        sums.emplace_back(name, crc32(v.data(), v.size() * sizeof(float), crc32(&N, sizeof N)));
        return true;
    };
    if (!r.error.empty() || r.depth_m.empty()) { if (err) *err = "no retention product to write"; return false; }
    if (!put("ice_depth.bin", r.depth_m) || !put("ice_table_layer.bin", r.layer) ||
        !put("ice_retention_class.bin", r.cls) || !put("ice_toplayer_lifetime.bin", r.toplayer_gyr) ||
        !put("ice_column_lifetime.bin", r.column_gyr) || !put("ice_exposed_loss_tau.bin", r.exposed_loss_m))
        return false;
    std::ofstream os(dir + "/ice_retention.txt");
    os << prov
       << "files: ice_depth.bin (m; ice-table depth at the age, z_bot where none), "
          "ice_table_layer.bin (layer index, -1 none), ice_retention_class.bin "
          "(0 top layer, 1 at depth, 2 none within the column, 3 top layer with surface rate "
          "below float range, 255 not evaluated), ice_toplayer_lifetime.bin (Gyr), "
          "ice_column_lifetime.bin (Gyr), ice_exposed_loss_tau.bin (m of exposed ice lost in the age); "
          "all uint32 N + N float32\n";
    for (const auto& s : sums) {
        char hex[16];
        std::snprintf(hex, sizeof hex, "%08x", static_cast<unsigned>(s.second));
        os << "file_crc32: " << s.first << " " << hex << "\n";
    }
    if (!os.good()) { if (err) *err = "write failed: " + dir + "/ice_retention.txt"; return false; }
    return true;
}

const std::vector<std::string>& iceProductFileNames() {
    static const std::vector<std::string> names = {
        "ice_depth.bin", "ice_table_layer.bin", "ice_retention_class.bin",
        "ice_toplayer_lifetime.bin", "ice_column_lifetime.bin", "ice_exposed_loss_tau.bin",
        "ice_retention.txt", "ice_depth_flux.bin", "ice_depth_flux_threshold.txt",
        "ice_depth_threshold.txt", "ice_depth_avg.bin"};
    return names;
}

const std::vector<std::string>& iceRawFileNames() {
    static const std::vector<std::string> names = {
        "ice_depth_sum_E.bin", "ice_depth_n_samples.txt", "ice_depth_params.txt",
        "ice_column_H.bin", "depth_grid.bin"};
    return names;
}

std::vector<std::string> listIceProducts(const std::string& dir, bool include_raw) {
    std::vector<std::string> out;
    auto present = [](const std::string& p) { std::ifstream is(p); return static_cast<bool>(is); };
    for (const auto& n : iceProductFileNames())
        if (present(dir + "/" + n)) out.push_back(n);
    if (include_raw)
        for (const auto& n : iceRawFileNames())
            if (present(dir + "/_raw/" + n)) out.push_back("_raw/" + n);
    return out;
}

std::vector<std::string> removeIceProducts(const std::string& dir, bool include_raw,
                                           std::string* err) {
    std::vector<std::string> removed;
    for (const auto& n : listIceProducts(dir, include_raw)) {
        const std::string p = dir + "/" + n;
        if (std::remove(p.c_str()) == 0) removed.push_back(n);
        else if (err) *err += (err->empty() ? "" : ", ") + p;
    }
    return removed;
}

bool writeIceProductRecord(const std::string& path, const IceProductRecord& rec) {
    std::ofstream os(path);
    if (!os) return false;
    os.precision(17);
    // The two legacy lines first: a reader that stops at the first
    // non-numeric value still gets them.
    os << "flux_threshold_kgm2s " << rec.flux_threshold_kgm2s << "\n"
       << "grain_hop_length_m " << rec.hop_length_m << "\n"
       << "flux_product " << (rec.flux_requested ? 1 : 0) << "\n";
    if (rec.has_retention) {
        const auto& p = rec.retention;
        os << "retention_age_gyr " << p.age_gyr << "\n"
           << "ice_weight_fraction " << p.weight_fraction << "\n"
           << "grain_density_kgm3 " << p.grain_density_kgm3 << "\n"
           << "ice_density_kgm3 " << p.ice_density_kgm3 << "\n"
           << "rho_surface_kgm3 " << p.rho_surface_kgm3 << "\n"
           << "rho_deep_kgm3 " << p.rho_deep_kgm3 << "\n"
           << "density_scale_m " << p.density_scale_m << "\n"
           << "per_face_H " << (rec.per_face_H ? 1 : 0) << "\n"
           << "ice_inventory " << iceInventoryName(p.inventory) << "\n";
    }
    os << "criterion " << rec.criterion << "\n";
    return os.good();
}

bool readIceProductRecord(const std::string& path, IceProductRecord& rec) {
    std::ifstream is(path);
    if (!is) return false;
    std::map<std::string, std::string> kv;
    std::string line;
    while (std::getline(is, line)) {
        std::istringstream ls(line);
        std::string k, v;
        if (ls >> k >> v) kv[k] = v;
    }
    auto num = [&](const char* k, double& out) {
        auto it = kv.find(k);
        if (it == kv.end()) return false;
        try { out = std::stod(it->second); } catch (...) { return false; }
        return true;
    };
    rec = IceProductRecord{};
    const bool have_legacy = num("flux_threshold_kgm2s", rec.flux_threshold_kgm2s) &&
                             num("grain_hop_length_m", rec.hop_length_m);
    double fp = 0.0;
    rec.flux_requested = num("flux_product", fp) && fp != 0.0;
    auto& p = rec.retention;
    p.hop_length_m = rec.hop_length_m;
    double pfh = 0.0;
    rec.has_retention = num("retention_age_gyr", p.age_gyr) &&
                        num("ice_weight_fraction", p.weight_fraction) &&
                        num("grain_density_kgm3", p.grain_density_kgm3) &&
                        num("ice_density_kgm3", p.ice_density_kgm3) &&
                        num("rho_surface_kgm3", p.rho_surface_kgm3) &&
                        num("rho_deep_kgm3", p.rho_deep_kgm3) &&
                        num("density_scale_m", p.density_scale_m) &&
                        num("per_face_H", pfh) &&
                        kv.count("ice_inventory") && parseIceInventory(kv["ice_inventory"], p.inventory);
    rec.per_face_H = rec.has_retention && pfh != 0.0;
    // A record written before the retention criterion existed describes the
    // flux product only.
    rec.criterion_recorded = kv.count("criterion") > 0;
    rec.criterion = rec.criterion_recorded ? kv["criterion"] : std::string("flux");
    if (!rec.criterion_recorded) rec.flux_requested = true;
    return have_legacy || rec.has_retention;
}

bool writeFaceH(const std::string& path, const std::vector<double>& H) {
    std::ofstream os(path, std::ios::binary);
    const std::uint32_t N = static_cast<std::uint32_t>(H.size());
    os.write(reinterpret_cast<const char*>(&N), sizeof N);
    os.write(reinterpret_cast<const char*>(H.data()), H.size() * sizeof(double));
    return os.good();
}

bool readFaceH(const std::string& path, std::vector<double>& H) {
    std::ifstream is(path, std::ios::binary);
    std::uint32_t N = 0;
    if (!is.read(reinterpret_cast<char*>(&N), sizeof N)) return false;
    H.resize(N);
    return static_cast<bool>(is.read(reinterpret_cast<char*>(H.data()), N * sizeof(double)));
}

} // namespace ice
