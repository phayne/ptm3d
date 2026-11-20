// test_flux_cube_pass.cpp
// ---------------------------------------------------------------------------
// Unit test for the flux-cube post-solve pass, on SYNTHETIC geometry whose
// answers are closed forms. No mesh, no Embree, no SPICE, no illumination
// model -- the pass takes its ray casting and primary-body physics as
// callbacks, so the test supplies deterministic ones and knows the truth.
//
// What each arm can falsify:
//
//  A. SEALED NADIR (one partner face, F = 1, no indirect flux)
//       Q_terrain_ir  == eps*sigma*T^4          exactly
//       Q_terrain_vis == A_src * Q_dir_face     exactly
//     Falsifies the reconstruction wiring, the row dot, and the component
//     assignment. This is the "905 W/m^2 under a lander" anchor from the design.
//
//  B. INDIRECT TERM. Same cell, now with a non-zero absorbed indirect flux:
//       Q_terrain_ir == eps*sigma*T^4 + (1-eps)*(Q_ind_ir/eps)
//     Falsifies dropping or mis-weighting the reflected part -- the term a
//     naive implementation forgets, and the one that carries PSR heat.
//
//  C. SHADOW PROBE. A synthetic occluder blocks the four rim probes but not the
//     center. V_sun must be 1/5 and the penumbral flag must be SET. An
//     implementation that probes only the center reports V_sun = 1 and a clear
//     flag, and this arm fails it.
//
//  D. DIRECT-SOLAR ORIENTATION. With the Sun along +U, the zenith face must see
//     S*V_sun and the nadir face exactly 0. Falsifies a sign flip or a face
//     mis-assignment that arm A cannot see (arm A has no direct beam).
//
//  E. CHUNKED WINDOW (positive control). A chunked run integrates one grid
//     point past its export window (thermal.checkpoint_lookahead_steps: 1), so
//     equilibration + output is one SHORT of total_time_steps. That window MUST
//     run, and it must read the FIRST n_steps epochs after equilibration -- not
//     the look-ahead step, and not an epoch window shifted by one. The fixture
//     ramps T with absolute step so a shifted read gives a different number:
//     10 K per step, against a 1e-9 K check.
//     Before the fix this arm was refused outright, which is why flux cubes
//     could not be produced for any chunked window.
//
//  J. TIME-AXIS REFUSAL (negative control). An output window that runs PAST the
//     end of the integrated series -- equilibration_steps + n_steps >
//     total_time_steps -- MUST still be refused. This is the case the check was
//     written for: accepting it reads past the end of every series and cube
//     rows silently stop lining up with time_steps.txt (the output-epoch drift defect
//     class). Relaxing == to <= must not cost this refusal.
//
//  F. CACHE-KEY REFUSAL (negative control). A cache built for a different mesh
//     MUST be refused, never rebuilt or applied anyway.
//
//  G. GRAY-PATH REFUSAL (negative control). Missing band-separated series MUST
//     be refused rather than emitting zeros or a merged column.
//
//  L. MOVED NAMED LOCATION (negative controls plus positive controls). A
//     location whose configured lat/lon, height_m or edge_m differs from the
//     cache MUST be refused; the same coordinates, a 5 mm difference inside the
//     1 cm tolerance, and a name-only entry MUST run. Before the fix the pass
//     compared names only and accepted every case.
//
//  M. SHORT EPHEMERIS (negative controls). A Sun or primary-body series that
//     ends before the window does MUST be refused. Before the fix the pass
//     read the uncovered epochs as "no sun" and reported zero direct flux.
//
// Exit code = number of failures.
// ---------------------------------------------------------------------------

#include "flux_cube_format.h"
#include "flux_cube_pass.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <cstddef>

using namespace ptm::fluxcubes;

namespace {

constexpr double SIGMA = 5.670374419e-8;

// Write a one-cube cache whose NADIR patch sees exactly one partner face with
// F = 1 (a sealed cell) and whose other five faces are pure sky.
std::string writeCache(const std::string& dir, int32_t partner_face,
                       double geomhash, uint64_t n_faces_core,
                       double height, double edge,
                       const std::string& name = std::string()) {
    std::filesystem::create_directories(dir);
    const std::string path = dir + "/synthetic.fcg";

    FcgHeader h{};
    std::memcpy(h.magic, "FCG1", 4);
    h.version = 1;
    h.n_cubes = 1;
    h.n_faces_per_cube = kNumFaces;
    h.total_nnz = 1;
    h.geomhash_core = geomhash;
    h.n_faces_core = n_faces_core;
    h.ring_hash = 0;
    h.center_lat_deg = 0.0; h.center_lon_deg = 0.0;
    h.spacing_m = 25.0; h.radius_m = 0.0;
    h.height_m = height; h.edge_m = edge;
    h.n_rays = 16384; h.sampler_id = kSamplerShirleyChiuCP; h.seed = 0;
    h.truncate_rel = 0.0;
    h.height_factor = 0.1; h.edge_factor = 0.1; h.local_median_edge_m = height / 0.1;
    // Site frame: U = +x, N = +z, E = N x U = +y (right-handed, E x N = U).
    h.frame_origin[0] = 1.0e6; h.frame_origin[1] = 0; h.frame_origin[2] = 0;
    h.frame_up[0]    = 1; h.frame_up[1]    = 0; h.frame_up[2]    = 0;
    h.frame_north[0] = 0; h.frame_north[1] = 0; h.frame_north[2] = 1;
    h.frame_east[0]  = 0; h.frame_east[1]  = 1; h.frame_east[2]  = 0;
    h.frame_max_tilt_deg = 0.0;
    const std::string json = "{\"schema\":1,\"tool\":\"test_flux_cube_pass\"}";
    h.json_len = static_cast<uint32_t>(json.size());
    h.has_names = name.empty() ? 0u : 1u;

    CubeRecord cr{};
    cr.grid_i = 0; cr.grid_j = 0; cr.x_m = 0; cr.y_m = 0;
    cr.center[0] = 1.0e6 + height; cr.center[1] = 0; cr.center[2] = 0;
    cr.face_id = partner_face;
    cr.slope_deg = 0; cr.aspect_deg = 0; cr.area_m2 = 1.0f; cr.elev_m = 0;
    cr.height_achieved_m = height;

    std::vector<PatchRecord> pr(kNumFaces);
    for (int f = 0; f < kNumFaces; ++f) {
        pr[f] = PatchRecord{};
        pr[f].row_offset = 0;
        pr[f].row_len = 0;
        pr[f].f_sky = 1.0f;
    }
    // Nadir: one partner, F = 1, no sky.
    pr[kFaceNadir].row_offset = 0;
    pr[kFaceNadir].row_len = 1;
    pr[kFaceNadir].f_sky = 0.0f;
    pr[kFaceNadir].f_max = 1.0f;
    pr[kFaceNadir].n_faces_50pct = 1;

    std::ofstream o(path, std::ios::binary);
    o.write(reinterpret_cast<const char*>(&h), sizeof(h));
    o.write(reinterpret_cast<const char*>(&cr), sizeof(cr));
    o.write(reinterpret_cast<const char*>(pr.data()),
            static_cast<std::streamsize>(pr.size() * sizeof(PatchRecord)));
    const int32_t col = partner_face;
    const float val = 1.0f;
    o.write(reinterpret_cast<const char*>(&col), sizeof(col));
    o.write(reinterpret_cast<const char*>(&val), sizeof(val));
    o.write(json.data(), static_cast<std::streamsize>(json.size()));
    if (!name.empty()) {
        char nb[kNameLen] = {0};
        std::strncpy(nb, name.c_str(), kNameLen - 1);
        o.write(nb, kNameLen);
    }
    o.close();
    return path;
}

struct Fcb {
    FcbHeader h{};
    std::vector<double> flux;
    std::vector<double> t_face;
    std::vector<uint8_t> penumbral;
};

bool readFcb(const std::string& path, Fcb& f) {
    std::ifstream i(path, std::ios::binary);
    if (!i) return false;
    i.read(reinterpret_cast<char*>(&f.h), sizeof(f.h));
    if (std::memcmp(f.h.magic, "FCB1", 4) != 0) return false;
    const size_t n = static_cast<size_t>(f.h.n_cubes) * f.h.n_faces_per_cube *
                     f.h.n_components * f.h.n_steps;
    f.flux.resize(n);
    i.read(reinterpret_cast<char*>(f.flux.data()),
           static_cast<std::streamsize>(n * sizeof(double)));
    f.t_face.resize(static_cast<size_t>(f.h.n_cubes) * f.h.n_steps);
    i.read(reinterpret_cast<char*>(f.t_face.data()),
           static_cast<std::streamsize>(f.t_face.size() * sizeof(double)));
    f.penumbral.resize(static_cast<size_t>(f.h.n_cubes) * f.h.n_steps);
    i.read(reinterpret_cast<char*>(f.penumbral.data()),
           static_cast<std::streamsize>(f.penumbral.size()));
    return static_cast<bool>(i);
}

double at(const Fcb& f, int cube, int face, int comp, int step) {
    return f.flux[((static_cast<size_t>(cube) * f.h.n_faces_per_cube + face) *
                   f.h.n_components + comp) * f.h.n_steps + step];
}

// ---- the run fixture -------------------------------------------------------
struct Fixture {
    static constexpr size_t NF = 4;      // synthetic mesh faces
    static constexpr int EQ = 3;         // equilibration steps
    static constexpr int NS = 2;         // output steps
    static constexpr int TT = EQ + NS;   // integrated steps with NO look-ahead
    // A chunked window integrates `lookahead` grid points past its export
    // window so the checkpoint it hands the next chunk is the state at that
    // chunk's FIRST epoch (thermal.checkpoint_lookahead_steps). Those points
    // are integrated and never exported, so the series are sized to tt() while
    // `times` stays NS long. Default 0 keeps every pre-existing arm unchanged.
    int lookahead = 0;
    int tt() const { return TT + lookahead; }
    int32_t partner = 2;
    double eps = 0.95, albedo = 0.20, T = 360.0;
    // Temperature as a function of ABSOLUTE step: T + dT_per_step * t. Default
    // 0 makes every step identical, which is what arms A-I assume. A non-zero
    // value is what lets the chunked arm tell the right epoch window from a
    // shifted one.
    double dT_per_step = 0.0;
    double tempAt(int t) const { return T + dT_per_step * static_cast<double>(t); }
    double q_dir_face = 0.0, q_ind_ir = 0.0, q_ind_vis = 0.0;
    double S = 1361.0;
    std::array<double, 3> sun{{0, 0, 0}};
    bool block_rim_only = false;
    bool block_all = false;

    std::vector<double> times{0.0, 21600.0};
    std::vector<std::string> tstr{"2028-08-03T22:21:00.000", "2028-08-04T04:21:00.000"};
    std::vector<double> direct;                      // flat [f*TT + t]
    std::vector<std::vector<double>> temps, isol, iir;
    std::vector<std::array<double, 3>> sunpos, pripos;

    void build() {
        const int n = tt();
        direct.assign(NF * static_cast<size_t>(n), 0.0);
        temps.assign(NF, std::vector<double>(n, T));
        isol.assign(NF, std::vector<double>(n, q_ind_vis));
        iir.assign(NF, std::vector<double>(n, q_ind_ir));
        for (size_t f = 0; f < NF; ++f) {
            for (int t = 0; t < n; ++t) {
                direct[f * static_cast<size_t>(n) + t] = q_dir_face;
                temps[f][t] = tempAt(t);
            }
        }
        sunpos.assign(n, sun);
        pripos.assign(n, {{0, 0, 0}});
    }

    ptm::fluxcubes::PassInputs inputs(const std::string& cache, const std::string& out,
                                      double geomhash) {
        PassInputs in;
        in.cache_file = cache;
        in.geomhash_core = geomhash;
        in.n_faces_core = NF;
        in.times = &times;
        in.time_strings = &tstr;
        in.equilibration_steps = EQ;
        in.total_time_steps = tt();
        in.direct_solar = &direct;
        in.surface_temps = &temps;
        in.indirect_solar_series = &isol;
        in.indirect_ir_series = &iir;
        in.emissivity.assign(NF, eps);
        in.bond_albedo.assign(NF, albedo);
        in.sun_positions = &sunpos;
        in.primary_positions = &pripos;
        in.primary_radius_km = 0.0;
        const double Sv = S;
        in.solar_const = [Sv](int) { return Sv; };
        // No incidence-albedo function: A_src falls back to the Bond value, and
        // the pass records that in the header rather than assuming it silently.
        const bool rim = block_rim_only, all = block_all;
        // Origin x == the cube center x (1e6 + h) is the CENTER probe; the four
        // rim probes are displaced in y/z. Blocking exactly the displaced ones
        // is what makes arm C discriminate a 1-ray probe from a 5-ray one.
        in.occluded = [rim, all](const std::array<double, 3>& o,
                                 const std::array<double, 3>&) {
            if (all) return true;
            if (!rim) return false;
            return (std::fabs(o[1]) > 1e-12 || std::fabs(o[2]) > 1e-12);
        };
        in.output_dir = out;
        in.write_binary = true;
        in.write_csv = false;
        in.shadow_probe_rays = 5;
        return in;
    }
};

int check(const char* what, double got, double want, double tol, int& fails) {
    const double d = std::fabs(got - want);
    const bool ok = d <= tol;
    std::printf("  %-46s got %14.6f  want %14.6f  |d| %9.2e  %s\n",
                what, got, want, d, ok ? "PASS" : "FAIL");
    if (!ok) ++fails;
    return ok;
}

}  // namespace

int main() {
    std::printf("\n=== flux-cube post-solve pass: unit test on synthetic geometry ===\n");
    int fails = 0;
    const std::string base = "output/flux_cubes/unit";
    const double GH = 12345.678;

    // ---------------------------------------------------------------- ARM A
    {
        std::printf("\n--- ARM A: sealed nadir, no indirect flux ---\n");
        Fixture fx;
        fx.q_dir_face = 500.0;   // the FACE under the cube is sunlit
        fx.build();
        const std::string dir = base + "/armA";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        std::string err;
        if (!runFluxCubePass(in, err)) { std::printf("  FAIL: %s\n", err.c_str()); ++fails; }
        else {
            Fcb f;
            if (!readFcb(dir + "/flux_cubes/flux_cubes.fcb", f)) { std::printf("  FAIL: no .fcb\n"); ++fails; }
            else {
                const double want_ir = fx.eps * SIGMA * std::pow(fx.T, 4);
                const double want_vis = fx.albedo * fx.q_dir_face;   // A_src falls back to Bond
                check("nadir terrain_ir  = eps*sigma*T^4", at(f, 0, kFaceNadir, kCompTerrainIR, 0),
                      want_ir, 1e-9, fails);
                check("nadir terrain_vis = A*Q_dir_face", at(f, 0, kFaceNadir, kCompTerrainVis, 0),
                      want_vis, 1e-9, fails);
                check("zenith terrain_ir = 0 (pure sky)", at(f, 0, kFaceZenith, kCompTerrainIR, 0),
                      0.0, 1e-12, fails);
                check("T_face reported", f.t_face[0], fx.T, 1e-12, fails);
                std::printf("  (the design's lander-underside anchor: %.1f W/m^2 at T = %.1f K)\n",
                            want_ir, fx.T);
            }
        }
    }

    // ---------------------------------------------------------------- ARM B
    {
        std::printf("\n--- ARM B: the reflected-IR term (1-eps)*Q_ind_ir/eps ---\n");
        Fixture fx;
        fx.q_ind_ir = 40.0;
        fx.q_ind_vis = 7.0;
        fx.q_dir_face = 500.0;
        fx.build();
        const std::string dir = base + "/armB";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        std::string err;
        if (!runFluxCubePass(in, err)) { std::printf("  FAIL: %s\n", err.c_str()); ++fails; }
        else {
            Fcb f;
            readFcb(dir + "/flux_cubes/flux_cubes.fcb", f);
            const double want_ir = fx.eps * SIGMA * std::pow(fx.T, 4)
                                 + (1.0 - fx.eps) * (fx.q_ind_ir / fx.eps);
            const double want_vis = fx.albedo * fx.q_dir_face
                                  + fx.albedo * (fx.q_ind_vis / (1.0 - fx.albedo));
            const double naive_ir = fx.eps * SIGMA * std::pow(fx.T, 4);
            check("nadir terrain_ir  with reflected term", at(f, 0, kFaceNadir, kCompTerrainIR, 0),
                  want_ir, 1e-9, fails);
            check("nadir terrain_vis with scattered term", at(f, 0, kFaceNadir, kCompTerrainVis, 0),
                  want_vis, 1e-9, fails);
            std::printf("  (dropping the reflected term would give %.6f, off by %.4f W/m^2 -- "
                        "so this arm discriminates)\n", naive_ir, want_ir - naive_ir);
            if (std::fabs(want_ir - naive_ir) < 1e-6) {
                std::printf("  FAIL: arm B is not discriminating (the two spellings agree)\n");
                ++fails;
            }
        }
    }

    // ---------------------------------------------------------------- ARM C+D
    {
        std::printf("\n--- ARM C/D: 5-ray shadow probe, penumbral flag, and orientation ---\n");
        Fixture fx;
        fx.sun = {{1, 0, 0}};        // Sun straight up along +U
        fx.block_rim_only = true;    // only the four rim probes are occluded
        fx.build();
        const std::string dir = base + "/armC";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        std::string err;
        if (!runFluxCubePass(in, err)) { std::printf("  FAIL: %s\n", err.c_str()); ++fails; }
        else {
            Fcb f;
            readFcb(dir + "/flux_cubes/flux_cubes.fcb", f);
            // V_sun = 1 clear of 5 => zenith sees S/5; a center-only probe would
            // report the full S and no penumbral flag.
            check("zenith direct = V_sun*S, V_sun = 1/5", at(f, 0, kFaceZenith, kCompDirect, 0),
                  fx.S / 5.0, 1e-9, fails);
            check("nadir direct = 0 (Sun overhead)", at(f, 0, kFaceNadir, kCompDirect, 0),
                  0.0, 1e-12, fails);
            check("north direct = 0 (normal perp. to Sun)", at(f, 0, kFaceNorth, kCompDirect, 0),
                  0.0, 1e-12, fails);
            const bool pen = f.penumbral[0] != 0;
            std::printf("  %-46s got %14s  want %14s  %s\n", "penumbral flag SET",
                        pen ? "set" : "clear", "set", pen ? "PASS" : "FAIL");
            if (!pen) ++fails;
        }
    }

    // ---------------------------------------------------------------- ARM E
    {
        std::printf("\n--- ARM E: chunked window (checkpoint look-ahead), MUST run ---\n");
        Fixture fx;
        fx.lookahead = 1;         // the chunk integrates one unexported step
        fx.dT_per_step = 10.0;    // T varies with ABSOLUTE step, so a shift shows
        fx.build();
        const std::string dir = base + "/armE";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        std::printf("  time axis: eq %d + output %d = %d of %d integrated "
                    "(look-ahead %d)\n",
                    Fixture::EQ, Fixture::NS, Fixture::EQ + Fixture::NS,
                    fx.tt(), fx.lookahead);
        std::string err;
        if (!runFluxCubePass(in, err)) {
            std::printf("  FAIL: a chunked window was REFUSED: %s\n", err.c_str());
            ++fails;
        } else {
            Fcb f;
            if (!readFcb(dir + "/flux_cubes/flux_cubes.fcb", f)) {
                std::printf("  FAIL: no .fcb\n"); ++fails;
            } else if (static_cast<int>(f.h.n_steps) != Fixture::NS) {
                std::printf("  FAIL: n_steps %u, want %d (the look-ahead step must "
                            "not be exported)\n", f.h.n_steps, Fixture::NS);
                ++fails;
            } else {
                // The exported epochs are the FIRST NS after equilibration:
                // absolute steps EQ and EQ+1, never EQ+1 and EQ+2 (a shift into
                // the look-ahead) and never EQ-1 and EQ (a shift back into
                // equilibration).
                for (int k = 0; k < Fixture::NS; ++k) {
                    const int t_abs = Fixture::EQ + k;
                    char label[80];
                    std::snprintf(label, sizeof(label),
                                  "T_face[k=%d] = T(absolute step %d)", k, t_abs);
                    check(label, f.t_face[k], fx.tempAt(t_abs), 1e-9, fails);
                    std::snprintf(label, sizeof(label),
                                  "nadir terrain_ir[k=%d] = eps*sigma*T(%d)^4", k, t_abs);
                    check(label, at(f, 0, kFaceNadir, kCompTerrainIR, k),
                          fx.eps * SIGMA * std::pow(fx.tempAt(t_abs), 4), 1e-9, fails);
                }
                // Discrimination, stated rather than assumed: a one-step shift
                // would move T_face[0] by dT_per_step, which is 1e10 times the
                // 1e-9 tolerance above.
                const double shifted = fx.tempAt(Fixture::EQ + 1);
                std::printf("  (a one-step shift would report T_face[0] = %.3f K instead "
                            "of %.3f K -- so this arm discriminates)\n",
                            shifted, fx.tempAt(Fixture::EQ));
                if (std::fabs(shifted - fx.tempAt(Fixture::EQ)) < 1e-6) {
                    std::printf("  FAIL: arm E is not discriminating (dT_per_step is 0)\n");
                    ++fails;
                }
            }
        }
    }

    // ------------------------------------------------------- NEGATIVE CONTROLS
    std::printf("\n--- NEGATIVE CONTROLS (each MUST be refused) ---\n");
    {
        // The output window runs PAST the end of the integrated series. The
        // fixture's series are sized to tt() = TT = 5, so accepting this would
        // read absolute step 5 out of a stride-5 array -- the refusal is what
        // stands between this input and an out-of-bounds read.
        Fixture fx; fx.build();
        const std::string dir = base + "/armJ";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        std::vector<double> long_times{0.0, 21600.0, 43200.0};
        std::vector<std::string> long_tstr{"2028-08-03T22:21:00.000",
                                           "2028-08-04T04:21:00.000",
                                           "2028-08-04T10:21:00.000"};
        in.times = &long_times;             // eq 3 + output 3 = 6 > total 5
        in.time_strings = &long_tstr;
        std::string err;
        const bool ran = runFluxCubePass(in, err);
        std::printf("  %-46s %s\n", "J: output window past the end of the series",
                    ran ? "FAIL (accepted!)" : "PASS (refused)");
        if (ran) ++fails; else std::printf("      reason: %.140s\n", err.c_str());
    }
    {
        Fixture fx; fx.build();
        const std::string dir = base + "/armF";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH * 2.0);   // wrong mesh identity
        std::string err;
        const bool ran = runFluxCubePass(in, err);
        std::printf("  %-46s %s\n", "F: cache built for a different mesh",
                    ran ? "FAIL (accepted!)" : "PASS (refused)");
        if (ran) ++fails; else std::printf("      reason: %.140s\n", err.c_str());
    }
    {
        Fixture fx; fx.build();
        const std::string dir = base + "/armG";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        in.indirect_ir_series = nullptr;             // the gray path
        std::string err;
        const bool ran = runFluxCubePass(in, err);
        std::printf("  %-46s %s\n", "G: gray path (no band-separated series)",
                    ran ? "FAIL (accepted!)" : "PASS (refused)");
        if (ran) ++fails; else std::printf("      reason: %.140s\n", err.c_str());
    }
    {
        Fixture fx; fx.build();
        const std::string dir = base + "/armH";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        in.occluded = nullptr;                       // no ray tracer
        std::string err;
        const bool ran = runFluxCubePass(in, err);
        std::printf("  %-46s %s\n", "H: no ray-occlusion function",
                    ran ? "FAIL (accepted!)" : "PASS (refused)");
        if (ran) ++fails; else std::printf("      reason: %.140s\n", err.c_str());
    }
    {
        Fixture fx; fx.eps = 1e-6; fx.build();       // eps below the guard
        const std::string dir = base + "/armI";
        const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
        auto in = fx.inputs(cache, dir, GH);
        std::string err;
        const bool ran = runFluxCubePass(in, err);
        std::printf("  %-46s %s\n", "I: unreconstructible face (eps -> 0)",
                    ran ? "FAIL (accepted!)" : "PASS (refused)");
        if (ran) ++fails; else std::printf("      reason: %.140s\n", err.c_str());
    }

    // ---------------------------------------------------------------- ARM M
    // Ephemeris shorter than the window (PTM-600). The Sun is overhead, so a
    // pass that ran would report zero direct flux at the uncovered epoch and a
    // nonzero value at the covered one; it must refuse instead.
    {
        std::printf("\n--- ARM M: Sun / primary ephemeris shorter than the window ---\n");
        struct Case { const char* tag; int sun_len; bool sun_null; int pri_len; };
        const int need = Fixture::EQ + Fixture::NS;
        const Case cases[] = {
            {"M1 Sun ephemeris one epoch short", need - 1, false, need},
            {"M2 no Sun ephemeris at all", 0, true, need},
            {"M3 primary ephemeris one epoch short", need, false, need - 1},
        };
        int idx = 0;
        for (const Case& cs : cases) {
            Fixture fx;
            fx.sun = {{1, 0, 0}};
            fx.build();
            fx.sunpos.resize(static_cast<size_t>(cs.sun_len));
            fx.pripos.resize(static_cast<size_t>(cs.pri_len));
            const std::string dir = base + "/armM" + std::to_string(idx++);
            const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
            auto in = fx.inputs(cache, dir, GH);
            if (cs.sun_null) in.sun_positions = nullptr;
            std::string err;
            const bool ran = runFluxCubePass(in, err);
            const bool ok = !ran && err.find("ephemeris") != std::string::npos;
            std::printf("  %-46s %s\n", cs.tag,
                        ok ? "PASS (refused: ephemeris)"
                           : (ran ? "FAIL (accepted!)" : "FAIL (refused for another reason)"));
            if (!ran) std::printf("      reason: %.160s\n", err.c_str());
            if (!ok) ++fails;
        }
    }

    // ---------------------------------------------------------------- ARM L
    // Named location MOVED in the config without rebuilding the cache. The
    // cache holds 'pad_A' at the site origin (lat 0, lon 0 in this frame, at
    // 1e6 m from the body center). The name check alone accepts every arm
    // below; only the position check separates them.
    //   L1  same coordinates, height and edge         -> must run
    //   L2  lat moved 1e-3 deg (17.45 m)              -> must be refused
    //   L3  moved 5 mm, inside the 1 cm tolerance     -> must run
    //   L4  height_m 2.0 vs the cache's 0.108         -> must be refused
    //   L5  edge_m 0.5 vs the cache's 0.108           -> must be refused
    //   L6  name only, no lat/lon                     -> must run (name check only)
    //   L7  height_m 2.0 in a multi-height run        -> must run (sizes not compared)
    //   L8  lat moved, multi-height run               -> must be refused (position is)
    {
        std::printf("\n--- ARM L: configured named-location position vs the cache ---\n");
        const double deg_per_m = 180.0 / (M_PI * 1.0e6);
        struct Case { const char* tag; bool has_ll; double lat; double h; double e; bool must_run;
                      bool multi_cache = false; };
        const Case cases[] = {
            {"L1 same position, height, edge", true, 0.0, 0.108, 0.108, true},
            {"L2 lat moved 17.45 m", true, 1.0e-3, -1.0, -1.0, false},
            {"L3 lat moved 5 mm (inside 1 cm)", true, 0.005 * deg_per_m, -1.0, -1.0, true},
            {"L4 height_m 2.0 vs cache 0.108", true, 0.0, 2.0, -1.0, false},
            {"L5 edge_m 0.5 vs cache 0.108", true, 0.0, -1.0, 0.5, false},
            {"L6 name only (no lat/lon)", false, 0.0, -1.0, -1.0, true},
            {"L7 height_m 2.0, multi-height run: not compared", true, 0.0, 2.0, -1.0, true, true},
            {"L8 moved 17.45 m, multi-height run: refused", true, 1.0e-3, 2.0, -1.0, false, true},
        };
        int idx = 0;
        for (const Case& cs : cases) {
            Fixture fx; fx.build();
            const std::string dir = base + "/armL" + std::to_string(idx++);
            std::filesystem::remove_all(dir);
            const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108,
                                                 "pad_A");
            auto in = fx.inputs(cache, dir, GH);
            in.expected_location_names = {"pad_A"};
            PassInputs::ExpectedLocation el;
            el.name = "pad_A";
            el.has_lat_lon = cs.has_ll;
            el.lat_deg = cs.lat;
            el.lon_deg = 0.0;
            el.height_m = cs.h;
            el.edge_m = cs.e;
            in.expected_locations = {el};
            in.compare_location_size = !cs.multi_cache;
            std::string err;
            const bool ran = runFluxCubePass(in, err);
            const bool ok = (ran == cs.must_run) &&
                            (ran || err.find("'pad_A'") != std::string::npos);
            std::printf("  %-46s %s\n", cs.tag,
                        ok ? (ran ? "PASS (ran)" : "PASS (refused, names pad_A)")
                           : (ran ? "FAIL (accepted!)" : "FAIL (refused)"));
            if (!ran) std::printf("      reason: %.220s\n", err.c_str());
            if (!ok) ++fails;
        }
    }

    // ---------------------------------------------------------------- ARM K
    // PTM-823: an open patch has nothing below its terrain. Here the mesh blocks
    // nothing at all; the body sphere R sits 1000 m below the cube, so its cap
    // edge is acos(R/r) = 2.56 deg below the horizontal. A Sun 5 deg below the
    // horizontal (toward +E) is behind the body: with body_blocks the east face
    // must get 0 (without it, S cos 5 deg: the leak). A Sun 1 deg below the
    // horizontal clears the sphere and must be bit-identical with and without.
    {
        std::printf("\n--- ARM K: open-patch body sphere (PTM-823) ---\n");
        const double R = 1.0e6 - 1000.0;
        auto sphere = [R](const std::array<double, 3>& o, const std::array<double, 3>& d) {
            const double b = o[0]*d[0] + o[1]*d[1] + o[2]*d[2];
            if (b >= 0.0) return false;
            const double c = o[0]*o[0] + o[1]*o[1] + o[2]*o[2] - R * R;
            const double disc = b * b - c;
            if (disc < 0.0) return false;
            const double sq = std::sqrt(disc);
            return (-b - sq) > 0.0 || (-b + sq) > 0.0;
        };
        auto runArm = [&](double el_deg, bool with_sphere, const std::string& tag, Fcb& f,
                          std::string& prov) {
            Fixture fx;
            const double el = el_deg * M_PI / 180.0;
            fx.sun = {{std::sin(el), std::cos(el), 0.0}};   // U = +x, E = +y
            fx.build();
            const std::string dir = base + "/armK_" + tag;
            std::filesystem::remove_all(dir);
            const std::string cache = writeCache(dir, fx.partner, GH, Fixture::NF, 0.108, 0.108);
            auto in = fx.inputs(cache, dir, GH);
            if (with_sphere) { in.body_blocks = sphere; in.body_sphere_radius = R; }
            std::string err;
            if (!runFluxCubePass(in, err)) { std::printf("  FAIL: %s\n", err.c_str()); ++fails; return false; }
            readFcb(dir + "/flux_cubes/flux_cubes.fcb", f);
            std::ifstream pv(dir + "/provenance.txt");
            prov.assign(std::istreambuf_iterator<char>(pv), std::istreambuf_iterator<char>());
            return true;
        };
        Fcb f_off, f_on, g_off, g_on; std::string p_off, p_on, q_off, q_on;
        if (runArm(-5.0, false, "m5_off", f_off, p_off) && runArm(-5.0, true, "m5_on", f_on, p_on)) {
            check("Sun -5 deg, no sphere: east direct = S cos 5 (the leak)",
                  at(f_off, 0, kFaceEast, kCompDirect, 0), 1361.0 * std::cos(5.0 * M_PI / 180.0), 1e-9, fails);
            check("Sun -5 deg, sphere: east direct = 0", at(f_on, 0, kFaceEast, kCompDirect, 0), 0.0, 0.0, fails);
            const bool rec = p_on.find("flux_cube_body_cull[flux_cubes]=on") != std::string::npos
                && p_on.find("cube ray(s) blocked by the body sphere") != std::string::npos
                && p_on.find("flux_cube_body_cull_effect[flux_cubes]=0 ") == std::string::npos
                && p_on.find("LEAK BAND") != std::string::npos;
            std::printf("  %-46s %s\n", "provenance: cull on, rays blocked, leak band", rec ? "PASS" : "FAIL");
            if (!rec) { ++fails; std::printf("%s\n", p_on.c_str()); }
            const bool rec_off = p_off.find("flux_cube_body_cull[flux_cubes]=not used") != std::string::npos;
            std::printf("  %-46s %s\n", "provenance without a sphere: 'not used'", rec_off ? "PASS" : "FAIL");
            if (!rec_off) ++fails;
        }
        if (runArm(-1.0, false, "m1_off", g_off, q_off) && runArm(-1.0, true, "m1_on", g_on, q_on)) {
            const bool same = g_off.flux.size() == g_on.flux.size()
                && std::memcmp(g_off.flux.data(), g_on.flux.data(), g_off.flux.size() * sizeof(double)) == 0
                && g_off.penumbral == g_on.penumbral;
            std::printf("  %-46s %s\n", "Sun -1 deg (clears the sphere): bit-identical", same ? "PASS" : "FAIL");
            if (!same) ++fails;
            check("Sun -1 deg: east direct = S cos 1", at(g_on, 0, kFaceEast, kCompDirect, 0),
                  1361.0 * std::cos(1.0 * M_PI / 180.0), 1e-9, fails);
        }
    }

    std::printf("\n==== flux-cube pass unit test: %d failure(s) ====\n\n", fails);
    return fails;
}
