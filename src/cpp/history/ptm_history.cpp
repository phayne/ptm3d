// ptm_history -- run (or resume) one streamed spin-up history and write the window's handoff state.
// History marcher, M1 part 2 (2026-11-09). See history_marcher.h for the method.
//
// usage: ptm_history JOB.yaml [--threads N] [--nice N]
//
// JOB.yaml (every path required; nothing is guessed):
//   history:
//     job_dir: <dir>                       checkpoints, PAUSE, the log of record
//     forcing: <file.ptmhf>                per-epoch absorbed + reflected direct (PTMHF001)
//     epoch_grid: <file>                   the SITE'S UNION epoch grid, raw float64 ET (all windows share it)
//     vf_csr: <file.csr>                   referee view factors (read-only)
//     grid_reference_state: <state.bin>    a PTMT v2 state written by the window runs' configuration;
//                                          its depth grid must equal the one built here (else FATAL)
//     handoff_state: <file>                output PTMT v2 state at window_start - fine_orbits lunations
//     window_start_et: <et>                must be an epoch of the forcing
//     fine_orbits: 2                       N_fine (window-run spin-up after the handoff)
//     history_lunations: 22                history length before the handoff
//     fine_last_lunations: 10              last lunations at the fine substep target
//     body: Moon                           built-in body (thermal::Planet)
//     geothermal_heat_flow: 0.022          W/m^2 (the window runs' value)
//     emissivity: 0.95
//     bond_albedo: 0.12
//     checkpoint_every_lunations: 1
//     threads: 0                           0 = OpenMP default; --threads overrides
//     nice: 0                              --nice overrides
//     allow_foreign_epoch_grid: false      true = run anyway on a forcing NOT on the union grid (LOUD)
// Exit: 0 handoff written; 75 paused (resume: remove job_dir/PAUSE if present, rerun the same command);
// anything else = error (message on stderr).
#include "history_marcher.h"
#include "thermal/thermal_state_io.h"
#include "ptm_build_info.h"

#include <yaml-cpp/yaml.h>
#include <omp.h>
#include <sys/resource.h>
#include <sys/stat.h>

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace ptm::history;

static int fail(const std::string& m) { std::fprintf(stderr, "[history] FATAL: %s\n", m.c_str()); return 2; }

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: ptm_history JOB.yaml [--threads N] [--nice N]\n"); return 1; }
    YAML::Node root;
    try { root = YAML::LoadFile(argv[1]); } catch (const std::exception& e) { return fail(std::string("cannot read job file: ") + e.what()); }
    const YAML::Node h = root["history"];
    if (!h) return fail("job file has no `history:` block");
    static const char* known[] = {"job_dir", "forcing", "epoch_grid", "vf_csr", "grid_reference_state", "handoff_state",
                                  "window_start_et", "fine_orbits", "history_lunations", "fine_last_lunations", "body",
                                  "geothermal_heat_flow", "emissivity", "bond_albedo", "checkpoint_every_lunations",
                                  "threads", "nice", "allow_foreign_epoch_grid"};
    for (auto it = h.begin(); it != h.end(); ++it) {
        const std::string k = it->first.as<std::string>(); bool ok = false;
        for (const char* q : known) ok |= (k == q);
        if (!ok) return fail("unknown key history." + k + " (a typo here would otherwise be silently ignored)");
    }
    auto req = [&](const char* k) -> YAML::Node { if (!h[k]) { std::fprintf(stderr, "[history] FATAL: history.%s is required\n", k); std::exit(2); } return h[k]; };
    const std::string job_dir = req("job_dir").as<std::string>(), forcing = req("forcing").as<std::string>(),
                      grid_path = req("epoch_grid").as<std::string>(), csr = req("vf_csr").as<std::string>(),
                      ref_state = req("grid_reference_state").as<std::string>(), handoff = req("handoff_state").as<std::string>();
    const double w_et = req("window_start_et").as<double>();
    const int fine_orbits = h["fine_orbits"] ? h["fine_orbits"].as<int>() : 2;
    const int hist_lun = h["history_lunations"] ? h["history_lunations"].as<int>() : 22;
    const int fine_last = h["fine_last_lunations"] ? h["fine_last_lunations"].as<int>() : 10;
    const std::string body = h["body"] ? h["body"].as<std::string>() : "Moon";
    const bool allow_foreign = h["allow_foreign_epoch_grid"] && h["allow_foreign_epoch_grid"].as<bool>();
    int threads = h["threads"] ? h["threads"].as<int>() : 0, nicev = h["nice"] ? h["nice"].as<int>() : 0;
    for (int a = 2; a + 1 < argc; a += 2) {
        if (!std::strcmp(argv[a], "--threads")) threads = std::atoi(argv[a + 1]);
        else if (!std::strcmp(argv[a], "--nice")) nicev = std::atoi(argv[a + 1]);
        else return fail(std::string("unknown option ") + argv[a]);
    }
    if (threads > 0) omp_set_num_threads(threads);
    if (nicev != 0 && setpriority(PRIO_PROCESS, 0, nicev) != 0) std::fprintf(stderr, "[history] WARNING: setpriority(%d) failed\n", nicev);
    ::mkdir(job_dir.c_str(), 0755);
    std::printf("[history] ptm_history %s (%s), job %s, %d threads, nice %d\n", PTM_BUILD_GIT_SHA, PTM_BUILD_VERSION, argv[1],
                omp_get_max_threads(), nicev);

    // ---- forcing and the union epoch grid -------------------------------------------------
    ForcingFile F; std::string err = F.open(forcing);
    if (!err.empty()) return fail(err);
    std::vector<double> grid;
    { std::ifstream g(grid_path, std::ios::binary | std::ios::ate); if (!g) return fail("cannot open epoch grid " + grid_path);
      const auto n = (size_t)g.tellg(); if (n % 8 || !n) return fail(grid_path + ": not a float64 list"); grid.resize(n / 8);
      g.seekg(0); g.read(reinterpret_cast<char*>(grid.data()), (std::streamsize)n); }
    const auto& et = F.et();
    size_t g0 = 0; while (g0 < grid.size() && std::fabs(grid[g0] - et[0]) > 1e-3) ++g0;
    bool on_grid = g0 < grid.size() && g0 + et.size() <= grid.size();
    for (size_t k = 0; on_grid && k < et.size(); ++k) on_grid = std::fabs(grid[g0 + k] - et[k]) <= 1e-3;
    if (!on_grid) {
        const std::string m = "the forcing's " + std::to_string(et.size()) + " epochs are NOT a contiguous run of the site's union "
                              "epoch grid " + grid_path + " (" + std::to_string(grid.size()) + " epochs). Every window of a site must "
                              "march on the union grid, or the site envelope falls through to the exact chord off-grid.";
        if (!allow_foreign) return fail(m + " Refusing (history.allow_foreign_epoch_grid: true overrides, loudly).");
        std::fprintf(stderr, "[history] WARNING *** %s RUNNING ANYWAY (allow_foreign_epoch_grid). ***\n", m.c_str());
    }
    size_t nW = 0; while (nW < et.size() && std::fabs(et[nW] - w_et) > 1e-3) ++nW;
    if (nW == et.size()) return fail("window_start_et is not an epoch of the forcing");
    const thermal::Planet probe(body);
    const double P = probe.day, spl = 120.0;
    MarcherConfig mc;
    mc.period_s = P; mc.samples_per_lunation = 120;
    if ((size_t)(fine_orbits * 120) > nW) return fail("fine_orbits reaches before the forcing");
    mc.ti = nW - (size_t)fine_orbits * 120;
    if ((size_t)hist_lun * 120 > mc.ti) return fail("history_lunations reaches before the forcing");
    mc.i0 = mc.ti - (size_t)hist_lun * 120;
    mc.fine_from = mc.ti - std::min<uint64_t>(mc.ti - mc.i0, (uint64_t)fine_last * 120);
    if (mc.i0 + (uint64_t)mc.seed_lunations * 120 > mc.ti) return fail("history shorter than the seed block");
    if (mc.ti + 2 >= F.numEpochs()) return fail("the forcing must extend >= 2 epochs past the handoff");
    // The history marches [i0, ti] and reads forcing to ti + 2: that span must be on the P/120 cadence (the
    // forcing may continue on another cadence afterwards, e.g. the window's own).
    for (size_t k = mc.i0 + 1; k <= mc.ti + 2; ++k) if (std::fabs((et[k] - et[k - 1]) / (P / spl) - 1.0) > 1e-6)
        return fail("forcing epochs are not on the P/120 history cadence at epoch " + std::to_string(k));
    mc.qb = h["geothermal_heat_flow"] ? h["geothermal_heat_flow"].as<double>() : probe.Qb;
    mc.emissivity = h["emissivity"] ? h["emissivity"].as<double>() : probe.emissivity;
    mc.a_bond = h["bond_albedo"] ? h["bond_albedo"].as<double>() : 0.12;
    mc.checkpoint_every = 120 * (uint64_t)(h["checkpoint_every_lunations"] ? h["checkpoint_every_lunations"].as<int>() : 1);
    mc.job_dir = job_dir;

    // ---- column: the engine's own grid, verified against the window runs' state -----------
    Column col; col.planet = thermal::Planet(body); col.planet.Qb = mc.qb; col.planet.emissivity = mc.emissivity;
    col.ref = std::make_unique<thermal::Profile>(col.planet, 0.0, col.config, false);
    col.nl = col.ref->nlayers; col.z.assign(col.ref->z().begin(), col.ref->z().end()); col.kc.assign(col.ref->kc().begin(), col.ref->kc().end());
    col.dz.resize(col.nl - 1); for (int l = 0; l + 1 < col.nl; ++l) col.dz[l] = col.z[l + 1] - col.z[l];
    thermal::state_io::Trailer rt; std::string terr;
    if (thermal::state_io::readTrailer(ref_state, rt, terr) != thermal::state_io::TrailerStatus::Ok)
        return fail("grid_reference_state " + ref_state + " has no valid v2 trailer: " + terr);
    if ((int)rt.z.size() != col.nl) return fail("depth grid has " + std::to_string(col.nl) + " layers, the reference state " + std::to_string(rt.z.size()));
    double gz = 0; for (int l = 0; l < col.nl; ++l) gz = std::max(gz, std::fabs(rt.z[l] - col.z[l]));
    if (gz > 1e-12) return fail("depth grid differs from the reference state's by " + std::to_string(gz) + " m");
    std::printf("[history] column: %s, %s, Qb %.4g, eps %.3g, A_bond %.3g; grid == reference state (max |dz| %.1e m)\n", body.c_str(),
                thermal::state_io::gridFingerprint(col.z.data(), col.z.size()).c_str(), mc.qb, mc.emissivity, mc.a_bond, gz);

    // ---- view factors --------------------------------------------------------------------------
    uint64_t csr_digest = 0;
    auto M = loadCsr(csr, F.numFaces(), csr_digest, err);
    if (!M) return fail(err);
    std::printf("[history] forcing %s: %" PRIu64 " epochs x %" PRIu64 " faces (digest %s); CSR %s (digest %s); history epochs %" PRIu64
                "..%" PRIu64 " (fine from %" PRIu64 "), handoff ET %.6f\n", forcing.c_str(), F.numEpochs(), F.numFaces(),
                hex(F.digest()).c_str(), csr.c_str(), hex(csr_digest).c_str(), mc.i0, mc.ti, mc.fine_from, et[mc.ti]);
    std::printf("[history] NOTE: no mesh is loaded by a forcing-file job, so the CSR's geometry hash is not checked against a "
                "mesh; the job identity carries the CSR digest instead.\n");

    // ---- identity (resume must match every field) ----------------------------------------------
    Identity id;
    id.set("code", std::string("history-marcher-v1 rule REC two-bounce ") + PTM_BUILD_GIT_SHA);
    id.set("forcing", hex(F.digest())); id.set("vf_csr", hex(csr_digest));
    id.set("epoch_grid", hex(digest(grid.data(), 8 * grid.size())) + (on_grid ? "" : " FOREIGN"));
    id.set("column", thermal::state_io::gridFingerprint(col.z.data(), col.z.size()) + " body " + body);
    char pb[256]; std::snprintf(pb, sizeof pb, "qb %.17g eps %.17g A %.17g i0 %llu ti %llu fine %llu P %.17g", mc.qb, mc.emissivity,
                                mc.a_bond, (unsigned long long)mc.i0, (unsigned long long)mc.ti, (unsigned long long)mc.fine_from, P);
    id.set("params", pb);

    Marcher mr(mc, col, *M, F, id);
    const int rc = mr.run();
    if (rc != 0) return rc;

    // ---- handoff ------------------------------------------------------------------------------
    const auto S = mr.stateFaceMajor(); const size_t nf = F.numFaces();
    std::vector<float> lm((size_t)col.nl * nf);
    for (size_t i = 0; i < nf; ++i) for (int l = 0; l < col.nl; ++l) lm[(size_t)l * nf + i] = (float)S[i * col.nl + l];
    thermal::state_io::Trailer tr; tr.z = col.z; tr.et_seconds = mr.handoffEt(); tr.output_cadence_s = P / 120.0;
    tr.body_name = body; tr.engine_git_sha = PTM_BUILD_GIT_SHA;
    if (!thermal::state_io::write(handoff, nf, (size_t)col.nl, lm, &tr, err)) return fail("handoff: " + err);
    { std::ofstream o(handoff + ".f64", std::ios::binary); o.write(reinterpret_cast<const char*>(S.data()), (std::streamsize)(8 * S.size()));
      if (!o) return fail("cannot write " + handoff + ".f64"); }
    std::printf("[history] handoff written: %s (PTMT v2, ET %.6f, %zu faces x %d layers) + %s.f64 (face-major float64)\n",
                handoff.c_str(), mr.handoffEt(), nf, col.nl, handoff.c_str());
    if (!on_grid) std::fprintf(stderr, "[history] WARNING *** this history was NOT on the site's union epoch grid ***\n");
    return 0;
}
