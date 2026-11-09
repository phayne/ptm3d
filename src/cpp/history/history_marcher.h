// history_marcher.h -- the streamed history marcher (Fermi build spec 2026-11-09; M1 part 2).
//
// Replaces equilibration orbits before a window with ~2 eclipse years of REAL forcing marched
// explicitly, then hands one column state to the window run (which spins up only N_fine orbits).
//
// The recommended configuration (Fermi REC, numpy-validated; coordinator decisions 2026-11-09):
//   seed      480 samples = every 3rd history epoch of the first 12 lunations; forcing F + Qv
//             plus ONE infrared sweep from instantaneous radiative-equilibrium walls,
//             IR = eps M (F + Qv + Qb); a periodic replay of that block at P/120 x3 with an
//             Aitken extrapolation (deepest node ratio) after the second pass; start state
//             T_0 = ((<F_b> + Qb)/(eps sigma))^(1/4), T(z) = T_0 + Qb int dz/kc.
//   march     every history epoch k: infrared from this epoch's surface, IR_k = eps M (eps sigma T0^4),
//             held over the step; visible TWO bounces, Qv = (1 - A) (b1 + b2), b1 = M s, b2 = M (A b1),
//             s = A_h Q_sun the reflected direct beam; forcing linear between epochs at substep
//             midpoints; substeps ceil(dt / target), target P/120 on old history, P/768 on the last
//             10 lunations. ONE canonical width-3 product per epoch: [eps sigma T0^4, s_{k+2}, A b1_{k+1}]
//             -> IR_k, b1_{k+2}, b2_{k+1} (visible runs one epoch ahead, so both ends of the step have
//             their visible term; visible does not depend on temperature).
//   column    the production CPU SIMD column (ThermalBatchSoA::advanceTimestep), loaded once per
//             history step (MEASURED 3x cheaper than reloading per substep at AGT scale).
//   handoff   a PTMT v2 state at the window start minus N_fine lunations, for PTM_SEED_STATE_FILE /
//             a warm start.
// Epoch grid: the history MUST lie on the site's union epoch grid (every window of a site shares one
// grid, so one .hzc envelope covers them with zero fall-through). A forcing whose epochs are not a
// contiguous run of that grid is refused unless explicitly allowed -- and then said loudly.
// Pause / resume: history_checkpoint.h (byte-identical resume; gate test_history_checkpoint).
//
// Forcing input (PTMHF001; the test path -- in production the illumination stage supplies the same
// two arrays per epoch in memory): little-endian char[8] "PTMHF001"; u64 n_epochs; u64 n_faces;
// f64 et[n_epochs]; then per epoch f64 absorbed_direct[n_faces], f64 reflected_direct[n_faces]
// (reflected = A_h(mu0) Q_sun, the visible SOURCE; absorbed = (1 - A_h) Q_sun + Earthshine etc.).
#pragma once

#include "history_checkpoint.h"
#include "radiation/radiation_exchange.h"
#include "thermal/thermal_model.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace ptm {
namespace history {

class ForcingFile {
public:
    std::string open(const std::string& path);     // "" or an error
    ~ForcingFile();
    uint64_t numEpochs() const { return ne_; }
    uint64_t numFaces() const { return nf_; }
    const std::vector<double>& et() const { return et_; }
    uint64_t digest() const { return digest_; }
    // Reads epoch k (absorbed, reflected), n_faces each. false on a short read.
    bool read(uint64_t k, double* absorbed, double* reflected) const;
private:
    int fd_ = -1; uint64_t ne_ = 0, nf_ = 0, data_off_ = 0, digest_ = 0;
    std::vector<double> et_;
};

// Referee CSR (int32 n, f64 geomhash, i64 indptr[n+1], i64 nnz, i32 cols, f32 vals), the engine's
// invariants (row sums <= 1 + 1e-6, no negative entries) enforced. digest = word hash of the file.
std::unique_ptr<radiation::SparseVFOperator> loadCsr(const std::string& path, uint64_t n_expected,
                                                     uint64_t& digest, std::string& err);

struct MarcherConfig {
    // physics
    double period_s = 29.530588 * 86400.0;   // synodic day (the target substeps are fractions of it)
    double emissivity = 0.95;                // IR emissivity (surface BC and the IR operator)
    double a_bond = 0.12;                    // visible operator albedo (scatter A, absorb 1 - A)
    double qb = 0.018;                       // geothermal flux, W/m^2
    // epochs (indices into the forcing grid)
    uint64_t i0 = 0, ti = 0, fine_from = 0;  // history start, handoff, first epoch at the fine target
    double old_target_div = 120.0, fine_target_div = 768.0;
    // seed (REC)
    int seed_lunations = 12, seed_stride = 3, seed_ir_sweeps = 1, seed_reps = 3;
    double seed_target_div = 120.0;
    bool seed_aitken = true;
    int samples_per_lunation = 120;
    // job
    uint64_t checkpoint_every = 120;
    std::string job_dir;
};

// The column: grid + statics from the engine's own Profile.
struct Column {
    thermal::Planet planet{"Moon"};
    thermal::Config config;
    std::unique_ptr<thermal::Profile> ref;
    int nl = 0;
    std::vector<double> z, dz, kc;
};

class Marcher {
public:
    Marcher(const MarcherConfig& c, const Column& col, const radiation::SparseVFOperator& M, const ForcingFile& F,
            const Identity& id);
    // Runs (or resumes) to the handoff. Returns 0 done, kExitPaused paused, other = error (message printed).
    int run();
    // The handoff state, face-major (n_faces x n_layers), valid after run() returned 0.
    std::vector<double> stateFaceMajor() const;
    double handoffEt() const { return F_.et()[c_.ti]; }
private:
    void seed();
    void step(uint64_t k);
    void columnStep(const std::vector<double>& fa, const std::vector<double>& fb, double dt, double target);
    void visBounces(const std::vector<double>& src, std::vector<double>& b1, std::vector<double>& vis);
    Checkpoint save(uint64_t next, const char* phase) const;
    bool restore(const Checkpoint& c);
    double tgeo() const;

    MarcherConfig c_; const Column& col_; const radiation::SparseVFOperator& M_; const ForcingFile& F_; Identity id_;
    size_t nf_ = 0, nb_ = 0;
    std::vector<double> T_;            // per batch [bi][layer][BS]
    std::vector<double> b1_next_;      // first visible bounce (incident) of epoch k+1
    std::vector<double> vis_k_;        // absorbed visible indirect of epoch k
    uint64_t k_ = 0;
};

}  // namespace history
}  // namespace ptm
