// ice_table.h - Ice stability depth solver (C++ implementation)
//
// Provides a high-throughput Newton solver for ice table depth using
// mean (Tavg) and maximum (Tmax) surface temperatures, mirroring the
// Python vectorized implementation in ptm_ice/ice_table_vectorized.py.
//
// Public API (namespace ice):
//   std::vector<double> computeIceDepths(const std::vector<double>& Tavg,
//                                        const std::vector<double>& Tmax,
//                                        double J0,
//                                        double delta = 1e-4,
//                                        double z0 = 0.0,
//                                        double l = 100e-6,
//                                        double zs = 0.1,
//                                        double Q = 0.02,
//                                        double k = 1e-3,
//                                        int    max_iter = 200,
//                                        double z_min = 0.0,
//                                        double z_max = 5.0,
//                                        std::vector<int>* iterations = nullptr,
//                                        bool   use_openmp = true);
//
// The solver enforces physical guard rails, sets NaN where Tavg > Tmax or
// no convergence within max_iter, and supports optional iteration counts.
//
// Notes:
//  - Temperatures are assumed in Kelvin.
//  - Returned depths are in meters.
//  - Input vectors must be equal length; no reallocation of inputs.
//  - Implementation minimizes temporaries; per-element Newton iteration.
//  - Optional OpenMP parallelism; guarded by compilation flags.
//
// Future integration steps (not implemented here):
//  - Call after thermal simulation to append ice depth to statistics.
//  - Stream processing for very large face counts (>1e7) if needed.

#ifndef ICE_TABLE_H
#define ICE_TABLE_H

#include <vector>
#include <cstddef>
#include <string>

namespace ice {

// ---------------------------------------------------------------------------
// Schorghofer & Taylor (2007) depth-resolved ice-stability model
//
// The thermal solver fills a per-face × per-layer accumulator of E(T(z,t))
// during the integration window. E is the Hertz–Knudsen sublimation NUMBER
// flux of bare ice (molecules m^-2 s^-1); μE is the mass flux (kg m^-2 s^-1).
// Post-solve, computeIceStabilityDepths() classifies each face:
//
//   * stable at the SURFACE (z* = 0) when free sublimation of exposed ice,
//     time-averaged over the window, is at most J₀:  μ⟨E(0)⟩ ≤ J₀;
//   * otherwise stable at the smallest depth z* with
//         J(z) = μ ℓ ⟨E(z)⟩ / (2 z) ≤ J₀        (S&T 2007, Eq. 11)
//     i.e. the vapor loss through the dry overburden falls to J₀;
//   * z* = z_max (the bottom node) when neither holds within the column
//     (right-censored: deeper than the column, or absent).
//
// ⟨E(z)⟩ is the time-averaged rate at depth z over the integration window.
// The S&T flux μℓ⟨E⟩/(2z) is valid for z ≫ ℓ and diverges as z → 0, so it
// cannot decide the surface; before 2026-11 the surface class was assigned
// only when ⟨E(0)⟩ was exactly 0 (float underflow, T < ~56 K throughout) and
// every other face meeting the surface criterion got a spurious depth of
// z₁/512 (a few μm).
//
// J₀ in physical units (ice density 917 kg m^-3): J₀ = 3×10⁻¹⁴ kg m⁻² s⁻¹ is
// 947 kg m⁻² Gyr⁻¹, a loss of about 1 m of ice per Gyr, and equals the free
// sublimation of bare ice at 114.4 K. The literature criterion of S&T (2007)
// and Paige et al. (2010) is 1 kg m⁻² Gyr⁻¹ = 3.17×10⁻¹⁷ kg m⁻² s⁻¹, about
// 1 mm of ice per Gyr, bare ice at 101.3 K. See iceThresholdSummary().
//
// The accumulator buffer is layer-major (sum_E[layer * n_faces + face])
// to match the GPU OpenCL kernel's natural coalesced-write layout.
// ---------------------------------------------------------------------------
struct IceStabilityAccumulator {
    // sum_E[layer * n_faces + face] — running sum of E(T_layer(t)) across
    // all integration timesteps. Allocated only when ice-stability is
    // enabled. float32 to match GPU buffer; ~1 GB at 4.8M faces × 50 layers.
    std::vector<float>  sum_E;
    int                 n_faces  = 0;
    int                 n_layers = 0;
    // n_samples is set by the host once the integration loop is done.
    // (Each substep adds one E sample per layer per face on the GPU path,
    // or one sample per output step on the CPU path — host knows which.)
    int                 n_samples = 0;
    // Native depth grid (m), copied from the thermal Profile after the solve.
    std::vector<double> z;

    void allocate(int faces, int layers) {
        n_faces  = faces;
        n_layers = layers;
        sum_E.assign(static_cast<std::size_t>(faces) * layers, 0.0f);
        z.assign(layers, 0.0);
        n_samples = 0;
    }
    bool ready() const {
        return !sum_E.empty() && n_samples > 0
            && static_cast<std::size_t>(n_faces) * n_layers == sum_E.size()
            && static_cast<int>(z.size()) == n_layers;
    }
};

// Compute per-face ice-stability depths from a finalized accumulator.
// Returns a vector of length acc.n_faces with depths in meters: 0 = stable at
// the surface, z_cap = not stable within the column. Returns an EMPTY vector
// when the accumulator is not ready or J0 / l_m is not a positive finite
// number; callers must treat that as fatal.
//   J0_kgm2s : loss threshold, kg m^-2 s^-1 (default 3e-14: about 1 m of ice
//              per Gyr; see the block comment above)
//   l_m      : grain hop length / pore diameter (default 75 μm)
//   z_max_m  : cap stability depth at this value (default = max depth in acc.z)
std::vector<float>
computeIceStabilityDepths(const IceStabilityAccumulator& acc,
                          double J0_kgm2s = 3.0e-14,
                          double l_m      = 75.0e-6,
                          double z_max_m  = -1.0,
                          bool   use_openmp = true);

// The engine's free-sublimation MASS flux of bare water ice at temperature T
// (kg m^-2 s^-1): μ·E(T), with the same constants as the solvers' accumulators
// (thermal_kernels.cl "sublimation-rate accumulator"), in double precision.
double sublimationMassFlux(double T_K);

// Temperature (K) at which bare ice sublimates at J0_kgm2s, i.e. the
// equivalent surface temperature of the threshold. NaN for J0 <= 0.
double surfaceStabilityTemperature(double J0_kgm2s);

// One line stating J0 in physical units, for logs and provenance:
// "J0 = 3e-14 kg m^-2 s^-1 = 947 kg m^-2 Gyr^-1 (1.03 m of ice per Gyr at
//  917 kg m^-3; bare ice at 114.4 K)".
std::string iceThresholdSummary(double J0_kgm2s);

constexpr double kWaterMoleculeMassKg = 2.99e-26;   // μ, as in the accumulators' E -> μE
constexpr double kSecondsPerGyr = 1.0e9 * 365.25 * 86400.0;
constexpr double kIceDensity    = 917.0;   // kg m^-3, for "ice thickness per Gyr"

// Compute ice stability depths given mean and max temperatures.
// See header comment for parameter definitions.
std::vector<double> computeIceDepths(const std::vector<double>& Tavg,
                                     const std::vector<double>& Tmax,
                                     double J0,
                                     double delta = 1e-4,
                                     double z0 = 0.0,
                                     double l = 100e-6,
                                     double zs = 0.1,
                                     double Q = 0.02,
                                     double k = 1e-3,
                                     int    max_iter = 200,
                                     double z_min = 0.0,
                                     double z_max = 5.0,
                                     std::vector<int>* iterations = nullptr,
                                     bool   use_openmp = true);

// Pointer/stride variant for direct buffer operation (in-place output).
// Returns 0 on success, non-zero on error.
int computeIceDepthsRaw(const double* Tavg,
                        const double* Tmax,
                        std::size_t N,
                        double J0,
                        double* depths_out,
                        int* iterations_out = nullptr,
                        double delta = 1e-4,
                        double z0 = 0.0,
                        double l = 100e-6,
                        double zs = 0.1,
                        double Q = 0.02,
                        double k = 1e-3,
                        int    max_iter = 200,
                        double z_min = 0.0,
                        double z_max = 5.0,
                        bool   use_openmp = true);

} // namespace ice

#endif // ICE_TABLE_H
