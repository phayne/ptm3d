// bench_thermal_batch -- wall per whole-mesh thermal substep of the engine's CPU SIMD batch column
// (thermal::ThermalBatchSoA::advanceTimestep, the production legacy step), at a chosen face count.
//
// History marcher M1 (2026-11-09), risk 3a: the history's thermal term is n_steps x n_sub
// whole-mesh substeps; its per-substep wall on the Mac Pro CPU was PREDICTED (0.0146-0.029 s at
// AGT scale), never measured. This tool measures it with the engine's own column code, on the
// engine's own grid (thermal::Planet + Config defaults + Profile, as dump_thermal_grid), and with
// the engine's per-batch state layout (timeOuterLoad/Save in ptm_thermal_sim.cpp: T, k, cp, kappa,
// rho, kc per layer, then emissivity, flux, Qb, T_eq_init).
//
// Two drivers, both timed:
//   engine  -- load + save the batch state around EVERY substep (the time-outer CPU driver today);
//   marcher -- load once per history step, n_sub substeps in cache, save once (what a streamed
//              history marcher does: the flux of a history step is fixed except for the linear
//              interpolation between its two epochs, which is computed in-batch).
// No mesh, no illumination, no GPU, no files written. Forcing is synthetic: half the faces dark,
// the rest 0-1300 W/m2 on a diurnal-like pattern, so cp(T) and k(T) span the real range.
//
// usage: bench_thermal_batch [n_faces=1797542] [n_steps=6] [n_sub=7] [threads=t1,t2,...]
#include "../thermal/thermal_model.h"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {
constexpr int BS = thermal::THERMAL_BATCH_SIZE;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// Identical to ptm_thermal_sim.cpp timeOuterSave / timeOuterLoad (static there).
void save(const thermal::ThermalBatchSoA& b, int nl, double* d) {
    const size_t L = static_cast<size_t>(nl) * BS;
    std::memcpy(d, b.T, L * 8); std::memcpy(d + L, b.k, L * 8); std::memcpy(d + 2 * L, b.cp, L * 8);
    std::memcpy(d + 3 * L, b.kappa, L * 8); std::memcpy(d + 4 * L, b.rho, L * 8); std::memcpy(d + 5 * L, b.kc, L * 8);
    double* e = d + 6 * L;
    std::memcpy(e, b.emissivity, BS * 8); std::memcpy(e + BS, b.flux, BS * 8);
    std::memcpy(e + 2 * BS, b.Qb, BS * 8); std::memcpy(e + 3 * BS, b.T_eq_init, BS * 8);
}
void load(thermal::ThermalBatchSoA& b, int nl, const double* d) {
    const size_t L = static_cast<size_t>(nl) * BS;
    std::memcpy(b.T, d, L * 8); std::memcpy(b.k, d + L, L * 8); std::memcpy(b.cp, d + 2 * L, L * 8);
    std::memcpy(b.kappa, d + 3 * L, L * 8); std::memcpy(b.rho, d + 4 * L, L * 8); std::memcpy(b.kc, d + 5 * L, L * 8);
    const double* e = d + 6 * L;
    std::memcpy(b.emissivity, e, BS * 8); std::memcpy(b.flux, e + BS, BS * 8);
    std::memcpy(b.Qb, e + 2 * BS, BS * 8); std::memcpy(b.T_eq_init, e + 3 * BS, BS * 8);
}
}  // namespace

int main(int argc, char** argv) {
    const size_t nf = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 1797542;
    const int nsteps = argc > 2 ? std::atoi(argv[2]) : 6;
    const int nsub = argc > 3 ? std::atoi(argv[3]) : 7;
    std::vector<int> tl;
    if (argc > 4) { std::string a = argv[4]; for (size_t p = 0; p <= a.size();) { size_t q = std::min(a.find(',', p), a.size());
        tl.push_back(std::atoi(a.substr(p, q - p).c_str())); p = q + 1; } }
    else tl.push_back(omp_get_max_threads());
    if (nf < 8 || nsteps < 1 || nsub < 1) { std::fprintf(stderr, "FATAL: bad arguments\n"); return 1; }

    thermal::Planet planet("Moon");
    thermal::Config config;   // engine defaults = the AGT configs (no thermal_m/n/b set there)
    thermal::Profile ref(planet, 0.0, config, false);
    const int nl = ref.nlayers;
    const double P = planet.day;   // solar day, s
    const double dt_step = P / 120.0, dt_sub = dt_step / nsub;
    const size_t nb = (nf + BS - 1) / BS, per = static_cast<size_t>(nl) * BS * 6 + 4 * BS;
    std::printf("Moon grid: %d layers (z_max %.4f m), day %.1f s; faces %zu (%zu batches), state %.2f GB; "
                "history step P/120 = %.1f s, %d substeps of %.1f s (P/%.0f)\n",
                nl, ref.z().back(), P, nf, nb, nb * per * 8 / 1e9, dt_step, nsub, dt_sub, P / dt_sub);
    if (nl != 23) std::printf("WARNING: %d layers, not the 23 of the AGT near-final runs\n", nl);

    std::vector<double> state(nb * per);
    std::vector<double> fa(nf), fb(nf);
    auto forcing = [&](int k, std::vector<double>& f) {
        for (size_t i = 0; i < nf; ++i) {
            const double u = std::fmod(0.6180339887498949 * (double)(i + 1), 1.0);
            const double ph = 2 * M_PI * (u + k / 120.0);
            f[i] = (u < 0.5) ? 1.5 : std::max(0.0, 1300.0 * std::sin(ph) * u) + 1.5;   // +IR floor
        }
    };
    {   // initial state: the engine's grid and properties, T set per face in 40-390 K
        #pragma omp parallel
        {
            thermal::ThermalBatchSoA* b = new thermal::ThermalBatchSoA();
            b->initialize(nl, BS, planet, config); b->copyGridFromProfile(ref);
            #pragma omp for schedule(static)
            for (size_t bi = 0; bi < nb; ++bi) {
                for (int f = 0; f < BS; ++f) {
                    const double u = std::fmod(0.6180339887498949 * (double)(bi * BS + f + 1), 1.0);
                    for (int l = 0; l < nl; ++l) b->T[l * BS + f] = 40.0 + 350.0 * u * std::exp(-l / 6.0) + 10.0;
                    b->emissivity[f] = planet.emissivity; b->Qb[f] = planet.Qb;
                }
                b->n_faces = (int)std::min<size_t>(BS, nf - bi * BS);
                thermal::updateMaterialPropertiesBatchSIMD(*b);
                save(*b, nl, &state[bi * per]);
            }
            delete b;
        }
    }
    const std::vector<double> state0 = state;
    for (int nt : tl) {
        omp_set_num_threads(nt);
        for (int mode = 0; mode < 2; ++mode) {
            state = state0;
            double t_force = 0, t_th = 0;
            for (int k = 0; k < nsteps; ++k) {
                const double a = now(); forcing(k, fa); forcing(k + 1, fb); t_force += now() - a;
                const double t0 = now();
                if (mode == 0) {   // engine: load/save around every substep
                    for (int s = 0; s < nsub; ++s) {
                        const double al = (s + 1.0) / nsub;
                        #pragma omp parallel
                        {
                            thread_local thermal::ThermalBatchSoA* b = nullptr;
                            if (!b) { b = new thermal::ThermalBatchSoA(); b->initialize(nl, BS, planet, config); b->copyGridFromProfile(ref); }
                            #pragma omp for schedule(dynamic, 64)
                            for (size_t bi = 0; bi < nb; ++bi) {
                                const int nfb = (int)std::min<size_t>(BS, nf - bi * BS);
                                b->n_faces = nfb; load(*b, nl, &state[bi * per]);
                                alignas(64) double fl[BS] = {0};
                                for (int f = 0; f < nfb; ++f) fl[f] = (1 - al) * fa[bi * BS + f] + al * fb[bi * BS + f];
                                b->setFluxContiguous(fl, nfb); b->advanceTimestep(dt_sub);
                                save(*b, nl, &state[bi * per]);
                            }
                        }
                    }
                } else {           // marcher: load once per history step
                    #pragma omp parallel
                    {
                        thread_local thermal::ThermalBatchSoA* b = nullptr;
                        if (!b) { b = new thermal::ThermalBatchSoA(); b->initialize(nl, BS, planet, config); b->copyGridFromProfile(ref); }
                        #pragma omp for schedule(dynamic, 64)
                        for (size_t bi = 0; bi < nb; ++bi) {
                            const int nfb = (int)std::min<size_t>(BS, nf - bi * BS);
                            b->n_faces = nfb; load(*b, nl, &state[bi * per]);
                            for (int s = 0; s < nsub; ++s) {
                                const double al = (s + 1.0) / nsub;
                                alignas(64) double fl[BS] = {0};
                                for (int f = 0; f < nfb; ++f) fl[f] = (1 - al) * fa[bi * BS + f] + al * fb[bi * BS + f];
                                b->setFluxContiguous(fl, nfb); b->advanceTimestep(dt_sub);
                            }
                            save(*b, nl, &state[bi * per]);
                        }
                    }
                }
                t_th += now() - t0;
            }
            double cs = 0; for (size_t bi = 0; bi < nb; ++bi) cs += state[bi * per];   // keep the work observable
            const double per_sub = t_th / (nsteps * (double)nsub);
            std::printf("threads %2d %-7s: %.4f s per whole-mesh substep (%.3f s per history step of %d); "
                        "2817 steps -> %.0f s; checksum %.6e\n",
                        nt, mode ? "marcher" : "engine", per_sub, t_th / nsteps, nsub, 2817 * t_th / nsteps, cs);
            std::fflush(stdout);
            (void)t_force;
        }
    }
    return 0;
}
