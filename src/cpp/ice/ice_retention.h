// ice_retention.h - Ice retention criterion: the ice table after a stated age.
//
// Ice fills (part of) the pore space of every layer of the column at the start
// of an age tau. It is lost by diffusion of vapor through the dry regolith that
// forms above it, at the time-mean rate the thermal model gives for each depth,
// so the ice table retreats downward. The product is the depth of the ice table
// at tau and the model layer that holds it, or "none within the column" when
// the column is exhausted before tau (right-censored: deeper than the column,
// or absent).
//
// Physics (Schorghofer & Taylor 2007, Sec. 4.1, Eqs. 11-12, in series with the
// kinetic resistance of the ice surface). With E the Hertz-Knudsen number flux
// of bare ice (molecules m^-2 s^-1), mu the mass of a water molecule, ell(z)
// the vapor transport length (D_K = v ell / 8, the normalization the engine
// has always used) and Lambda(z) = int_0^z dz'/ell(z'):
//
//     <J(z)> = mu <E(z)> / (1 + 2 Lambda(z))                            (2)
//     m(z) dz/dt = <J(z)>               (moving boundary, quasi-steady)  (3)
//     t(z) = int_0^z m(z') [1 + 2 Lambda(z')] / (mu <E(z')>) dz'         (4)
//     t(z_tau) = tau                                                     (5)
//
// m(z) is the initial ice mass per unit bulk volume. <E(z)> is the time mean of
// the rate at the temperature of depth z over the integration window (the
// accumulator of ice_table.h; never E(<T>)).
//
// Discretization. Layer j is the control volume of node j: layer 0 is
// [0, z_1/2], layer j is [(z_{j-1}+z_j)/2, (z_j+z_{j+1})/2], the last layer
// ends at z_bot. Between nodes, log <E> is linear in z. Each node interval is
// split into kRetentionSubsteps equal sub-steps; m and ell are evaluated at the
// sub-step midpoints (so a material discontinuity on a sub-step boundary is
// exact), Lambda is accumulated exactly for piecewise-constant ell, and the
// integrand of (4), (A + B x) exp(-c x) on a sub-step, is integrated in closed
// form; the root of (5) is found within the sub-step by a bracketed Newton
// iteration on that closed form. The integral is exact for a rate whose
// logarithm is linear between nodes, whatever the rate drop across an
// interval; the remaining error is the log-linear interpolation of <E>
// between nodes, second order in the node spacing.
//
// Material assumptions live in ONE place, makeIceColumn(): the inventory m(z)
// and the transport length ell(z), both from the regolith density profile of
// the thermal model. A later porosity/permeability-versus-overburden model
// replaces only that function; the retreat solver, outputs and tests do not
// change.

#ifndef ICE_RETENTION_H
#define ICE_RETENTION_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace ice {

// Sub-steps per node interval (even, so that every layer boundary, the
// midpoint between two nodes, is a sub-step boundary).
constexpr int kRetentionSubsteps = 16;

// <E> below this (molecules m^-2 s^-1) is taken as this value inside the
// logarithm. A time-mean sum of exactly 0 means the float rate underflowed
// (T below about 56 K throughout the window); a rate this small has a
// lifetime of order 1e200 Gyr, so the floor never decides a class. Faces whose
// SURFACE rate underflowed are reported as their own class.
constexpr double kRetentionRateFloor = 1.0e-200;

// Version tag of this criterion, written to every provenance record. v2: the
// retreat integral in closed form between sub-steps, and a record that carries
// a checksum of every product file.
constexpr const char* kIceRetentionVersion = "retention-v2 (2026-11-18)";

enum class IceInventory { PoreFilling, WeightFraction };

const char* iceInventoryName(IceInventory inv);           // "pore_filling" | "weight_fraction"
bool parseIceInventory(const std::string& s, IceInventory& out);

struct IceRetentionParams {
    double       age_gyr            = 1.0;     // tau
    IceInventory inventory          = IceInventory::PoreFilling;
    double       weight_fraction    = 0.056;   // ice / (ice + dry regolith), WeightFraction only
    double       grain_density_kgm3 = 3100.0;  // rho_g; porosity phi = 1 - rho/rho_g
    double       ice_density_kgm3   = 917.0;   // rho_ice (pore filling; exposed-ice thickness)
    double       hop_length_m       = 75.0e-6; // ell, constant (the engine's grain_hop_length_um)
    // Regolith bulk density of the thermal model:
    //   rho(z) = rho_deep + (rho_surface - rho_deep) exp(-z / H)
    // H may be overridden per face (h_parameter_map / thermal_inertia_map).
    double       rho_surface_kgm3   = -1.0;
    double       rho_deep_kgm3      = -1.0;
    double       density_scale_m    = -1.0;    // H (scalar; the per-face map, if any, is separate)
};

// Empty string when the parameters are usable; otherwise the reason.
std::string validateIceRetentionParams(const IceRetentionParams& p);

// Material of one column on the solver's quadrature grid.
struct IceColumnProfile {
    std::vector<double> zq;      // quadrature points: the nodes plus kRetentionSubsteps-1 per interval
    std::vector<double> m_cell;  // initial ice mass per bulk volume on each sub-step (kg m^-3); zq.size()-1
    std::vector<double> lambda;  // Lambda(z) = int_0^z dz/ell at each zq (dimensionless); zq.size()
};

// The quadrature grid for a node grid (z_nodes[0] must be 0, strictly increasing).
std::vector<double> retentionQuadratureGrid(const std::vector<double>& z_nodes);

// THE profile builder: inventory and transport length from the regolith
// density profile with scale H_face. Empty profile (and *err set) when the
// inventory is not physical (porosity outside (0, 1), or a weight fraction that
// needs more ice than the pore space holds).
IceColumnProfile makeIceColumn(const IceRetentionParams& p,
                               const std::vector<double>& z_nodes,
                               double H_face, std::string* err = nullptr);

// Generic builder: m(z) and ell(z) evaluated at the sub-step midpoints. Used by
// the analytic tests, and the entry point for a different material model.
IceColumnProfile makeIceColumnFromFunctions(const std::vector<double>& z_nodes,
                                            const std::function<double(double)>& m_of_z,
                                            const std::function<double(double)>& ell_of_z);

// Porosity of the thermal model's regolith at depth z.
double regolithPorosity(const IceRetentionParams& p, double z, double H);

enum RetentionClass : std::uint8_t {
    kRetainedTopLayer     = 0,   // the ice table at tau is in layer 0 (the surface layer)
    kRetainedAtDepth      = 1,   // in a deeper layer of the column
    kNoneWithinColumn     = 2,   // column exhausted before tau (right-censored)
    kRetainedUnderflow    = 3,   // surface <E> underflowed to 0 (T < ~56 K throughout): retained in layer 0
    kNotEvaluated         = 255  // a non-finite or negative <E>: no value
};

struct IceRetentionFace {
    double       z_tau_m      = 0.0;  // interpolated ice-table depth at tau; z_bot when none
    int          layer        = -1;   // layer holding the ice table; -1 none
    std::uint8_t cls          = kNotEvaluated;
    double       t_top_s      = 0.0;  // time to exhaust layer 0
    double       t_column_s   = 0.0;  // time to exhaust the whole column
    double       exposed_loss_m = 0.0;// thickness of exposed ice lost in tau: mu <E(0)> tau / rho_ice
};

// Solve one column. E_mean has one time-mean rate per node.
IceRetentionFace solveIceRetentionColumn(const std::vector<double>& z_nodes,
                                         const double* E_mean,
                                         const IceColumnProfile& col,
                                         const IceRetentionParams& p);

// Per-face products, all float32 in the statistics/ convention (uint32 N +
// N float32), so every reader of statistics/*.bin handles them unchanged.
struct IceRetentionResult {
    std::vector<float> depth_m;          // ice_depth.bin
    std::vector<float> layer;            // ice_table_layer.bin (-1 none, NaN no value)
    std::vector<float> cls;              // ice_retention_class.bin
    std::vector<float> toplayer_gyr;     // ice_toplayer_lifetime.bin
    std::vector<float> column_gyr;       // ice_column_lifetime.bin
    std::vector<float> exposed_loss_m;   // ice_exposed_loss_tau.bin
    // census
    long long n_top = 0, n_depth = 0, n_none = 0, n_underflow = 0, n_not_evaluated = 0;
    long long n_lifetime_beyond_float = 0;   // lifetimes written as FLT_MAX
    long long n_deep_node_underflow = 0;     // faces with an underflowed rate below the surface
    std::string error;                       // non-empty = no product
};

// Batch over the accumulator layout: sum_E[layer * n_faces + face], divided by
// n_samples. H_per_face (optional, n_faces entries, <= 0 = use the scalar H)
// gives each face its own density profile. One implementation for the engine
// (float sums) and the chunk merge (double sums).
IceRetentionResult computeIceRetention(const float* sum_E, int n_faces, int n_layers,
                                       long long n_samples, const std::vector<double>& z_nodes,
                                       const IceRetentionParams& p,
                                       const std::vector<double>* H_per_face = nullptr,
                                       bool use_openmp = true);
IceRetentionResult computeIceRetention(const double* sum_E, int n_faces, int n_layers,
                                       long long n_samples, const std::vector<double>& z_nodes,
                                       const IceRetentionParams& p,
                                       const std::vector<double>* H_per_face = nullptr,
                                       bool use_openmp = true);

// Human-readable provenance: criterion, age, inventory model and numbers,
// density profile, transport length, grid and top-layer thickness, samples,
// class census. One line per item, "key: value".
std::string iceRetentionProvenance(const IceRetentionParams& p,
                                   const std::vector<double>& z_nodes,
                                   long long n_samples,
                                   const IceRetentionResult& r,
                                   const std::vector<double>* H_per_face = nullptr);

// Write the per-face products into stats_dir (ice_depth.bin and the five
// companions) and stats_dir/ice_retention.txt (the provenance above, the file
// list, and one "file_crc32: <name> <8 hex digits>" line per product file).
// The checksums tie the record to the exact bytes it describes: a reader that
// finds ice_depth.bin with no record, or with a record whose checksum does not
// match, cannot take it for the retention product (a stale file from an
// earlier run or criterion; see removeIceProducts).
bool writeIceRetentionProducts(const std::string& stats_dir,
                               const IceRetentionResult& r,
                               const std::string& provenance_text,
                               std::string* err = nullptr);

// CRC-32 (IEEE 802.3, the zlib/PNG polynomial; Python's zlib.crc32).
std::uint32_t crc32(const void* data, std::size_t n, std::uint32_t crc = 0);

// ---------------------------------------------------------------------------
// Stale products. A run or merge into a directory that already holds ice
// products from an earlier run would leave them next to its own, and a reader
// would present them as this run's: an ice_depth.bin of the flux criterion, or
// a class file of another integration. Every writer of ice products therefore
// removes all of them first and says which.
// ---------------------------------------------------------------------------
// Product files in statistics/ (both criteria and the earlier names).
const std::vector<std::string>& iceProductFileNames();
// The run's saved ice sums in statistics/_raw/.
const std::vector<std::string>& iceRawFileNames();
// The product files (and, with include_raw, the _raw files as "_raw/<name>")
// present in stats_dir.
std::vector<std::string> listIceProducts(const std::string& stats_dir, bool include_raw);
// Remove them; returns the names removed. *err is set when one could not be
// removed (the caller must then refuse to write next to it).
std::vector<std::string> removeIceProducts(const std::string& stats_dir, bool include_raw,
                                           std::string* err);

// ---------------------------------------------------------------------------
// The record a run keeps in statistics/_raw/ice_depth_params.txt, so that the
// chunk merge and an offline re-evaluation use what the run was configured
// with. Shared by the engine (writer) and ptm_merge_chunks (reader).
// ---------------------------------------------------------------------------
struct IceProductRecord {
    std::string criterion = "retention";   // retention | flux
    bool   criterion_recorded = true;      // false: written before the retention criterion
    bool   flux_requested = false;         // legacy flux product also written
    double flux_threshold_kgm2s = 3.0e-14; // legacy J0
    double hop_length_m = 75.0e-6;
    bool   has_retention = false;          // the retention keys are present
    IceRetentionParams retention;
    bool   per_face_H = false;             // _raw/ice_column_H.bin carries a per-face H
};
bool writeIceProductRecord(const std::string& path, const IceProductRecord& rec);
// Returns false when the file is missing or has neither the legacy pair nor
// the retention keys. A file written before the retention criterion existed
// yields has_retention = false, criterion = "flux".
bool readIceProductRecord(const std::string& path, IceProductRecord& rec);

// Per-face density scale H (uint32 N + N float64), written beside the sums
// only when a per-face map set it.
bool writeFaceH(const std::string& path, const std::vector<double>& H);
bool readFaceH(const std::string& path, std::vector<double>& H);

} // namespace ice

#endif // ICE_RETENTION_H
