// flux_cube_pass.cpp -- see flux_cube_pass.h for the design rationale.
#include "flux_cube_pass.h"
#include <cstdint>

#include "flux_cube_format.h"
#include "../radiation/radiosity_reconstruction.h"
#include "../illumination/primary_disk.h"
#include "../illumination/solar_disk_visibility.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ptm {
namespace fluxcubes {
namespace {

constexpr double kPi = 3.14159265358979323846;

struct V3 {
    double x = 0, y = 0, z = 0;
};
inline V3 mk(const double* p) { return {p[0], p[1], p[2]}; }
inline V3 add(const V3& a, const V3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 sub(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 mul(const V3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot3(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double len3(const V3& a) { return std::sqrt(dot3(a, a)); }
inline V3 unit3(const V3& a) { const double n = len3(a); return n > 0 ? mul(a, 1.0 / n) : a; }
inline std::array<double, 3> arr(const V3& a) { return {a.x, a.y, a.z}; }

struct Cache {
    FcgHeader h{};
    std::vector<CubeRecord> cubes;
    std::vector<PatchRecord> patches;
    std::vector<int32_t> cols;
    std::vector<float> vals;
    std::string json;
    std::vector<std::string> names;   // empty, or one per cube
};

bool loadCache(const std::string& path, Cache& c, std::string& err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { err = "cannot open flux-cube cache '" + path + "'"; return false; }
    f.read(reinterpret_cast<char*>(&c.h), sizeof(c.h));
    if (!f || std::memcmp(c.h.magic, "FCG1", 4) != 0) {
        err = "'" + path + "' is not an FCG1 flux-cube cache";
        return false;
    }
    if (c.h.version != 1) {
        err = "flux-cube cache version " + std::to_string(c.h.version) + " != 1";
        return false;
    }
    if (c.h.n_faces_per_cube != kNumFaces) {
        err = "flux-cube cache has " + std::to_string(c.h.n_faces_per_cube) +
              " faces per cube, expected " + std::to_string(int(kNumFaces));
        return false;
    }
    c.cubes.resize(c.h.n_cubes);
    f.read(reinterpret_cast<char*>(c.cubes.data()),
           static_cast<std::streamsize>(c.cubes.size() * sizeof(CubeRecord)));
    c.patches.resize(static_cast<size_t>(c.h.n_cubes) * kNumFaces);
    f.read(reinterpret_cast<char*>(c.patches.data()),
           static_cast<std::streamsize>(c.patches.size() * sizeof(PatchRecord)));
    c.cols.resize(c.h.total_nnz);
    f.read(reinterpret_cast<char*>(c.cols.data()),
           static_cast<std::streamsize>(c.cols.size() * sizeof(int32_t)));
    c.vals.resize(c.h.total_nnz);
    f.read(reinterpret_cast<char*>(c.vals.data()),
           static_cast<std::streamsize>(c.vals.size() * sizeof(float)));
    c.json.resize(c.h.json_len);
    if (c.h.json_len) f.read(&c.json[0], c.h.json_len);
    if (!f) { err = "flux-cube cache '" + path + "' is truncated"; return false; }
    // Names, if present, sit AFTER the JSON -- see flux_cube_format.h.
    if (c.h.has_names) {
        std::vector<char> nb(static_cast<size_t>(c.h.n_cubes) * kNameLen, 0);
        f.read(nb.data(), static_cast<std::streamsize>(nb.size()));
        if (!f) {
            err = "flux-cube cache '" + path + "' declares named locations but its names "
                  "block is truncated";
            return false;
        }
        c.names.resize(c.h.n_cubes);
        for (uint32_t k = 0; k < c.h.n_cubes; ++k) {
            const char* p0 = &nb[static_cast<size_t>(k) * kNameLen];
            c.names[k].assign(p0, ::strnlen(p0, kNameLen));
        }
    }
    return true;
}

// Sample points on the primary body's disk as seen from `origin`. Returns unit
// directions. With radius_km <= 0 the primary is a POINT source and this
// returns the single center direction -- recorded in the header, not hidden.
// EVERY argument here is in KILOMETERS -- origin, position and radius alike.
// The returned directions are unit vectors and therefore frame-independent, so
// they may be cast in the mesh frame without conversion.
void sampleDisk(const V3& origin_km, const V3& primary_pos_km, double radius_km,
                int n_samples, std::vector<V3>& out) {
    out.clear();
    const V3 d = sub(primary_pos_km, origin_km);
    const double dist = len3(d);
    if (dist <= 0.0) { out.push_back(V3{0, 0, 1}); return; }
    const V3 c = mul(d, 1.0 / dist);
    if (radius_km <= 0.0 || n_samples <= 1 || dist <= radius_km) {
        out.push_back(c);
        return;
    }
    // Orthonormal basis spanning the disk.
    const V3 a = (std::fabs(c.z) < 0.9) ? V3{0, 0, 1} : V3{1, 0, 0};
    const V3 t = unit3({a.y * c.z - a.z * c.y, a.z * c.x - a.x * c.z, a.x * c.y - a.y * c.x});
    const V3 b = {c.y * t.z - c.z * t.y, c.z * t.x - c.x * t.z, c.x * t.y - c.y * t.x};
    // Sunflower (Vogel) spiral: equal-area, no clumping, deterministic.
    const double golden = kPi * (3.0 - std::sqrt(5.0));
    out.reserve(n_samples);
    for (int k = 0; k < n_samples; ++k) {
        const double r = radius_km * std::sqrt((k + 0.5) / n_samples);
        const double th = k * golden;
        const V3 p = add(add(mul(t, r * std::cos(th)), mul(b, r * std::sin(th))),
                         mul(c, dist));
        out.push_back(unit3(p));
    }
}

// Unoccluded fraction of the SOLAR DISK seen from `o`.
//
// The ladder used to live here. It now lives in
// illumination/solar_disk_visibility.h, because the shadow-depth pass asks the
// identical question from a point above the surface and two copies of it would
// be the drift class of the primary/eclipse fallbacks. This is a thin V3 adapter over that one
// spelling; `{0,0,0}` selects the historical disk basis, so the arithmetic is
// unchanged.
double solarDiskVisibility(const V3& o, const V3& s_hat, double alpha,
                           int samples,
                           const std::function<bool(const std::array<double,3>&,
                                                    const std::array<double,3>&)>& occluded,
                           bool& split) {
    return ptm::disk::solarDiskVisibility(arr(o), arr(s_hat), alpha, samples,
                                          ptm::disk::Vec3{0.0, 0.0, 0.0},
                                          occluded, split);
}

}  // namespace

bool runFluxCubePass(const PassInputs& in, std::string& err) {
    // ---------------------------------------------------------- preconditions
    if (!in.occluded) {
        err = "flux cubes: no ray-occlusion function was supplied. Refusing to run: "
              "assuming every cube is unshadowed would silently fabricate direct solar "
              "flux at every permanently shadowed site in the grid.";
        return false;
    }
    if (!in.solar_const) { err = "flux cubes: no solar-constant function supplied"; return false; }
    if (!in.times || !in.time_strings) { err = "flux cubes: no time axis supplied"; return false; }
    if (!in.direct_solar || !in.surface_temps ||
        !in.indirect_solar_series || !in.indirect_ir_series) {
        err = "flux cubes: a required field array is missing. The band-separated series "
              "(indirect_solar_series / indirect_ir_series) are produced ONLY by the "
              "two-band radiosity solve; with band_separated:false there is a single mixed "
              "radiosity and no visible/IR split, so a component-resolved product cannot be "
              "made. Set radiation.band_separated:true or disable flux_cubes -- this pass "
              "will not emit zeros or a merged column in place of the split.";
        return false;
    }

    Cache c;
    if (!loadCache(in.cache_file, c, err)) return false;

    // ------------------------------------------------------------- cache key
    // Refuse on mismatch and print the exact regeneration command. NEVER
    // rebuild inline: a 7,542-row hemisphere sweep starting unannounced inside
    // a production window is the kind of surprise that ends campaigns.
    const double gh_tol = 1e-3 * std::fabs(in.geomhash_core) + 1e-9;
    if (c.h.n_faces_core != in.n_faces_core ||
        std::fabs(c.h.geomhash_core - in.geomhash_core) > gh_tol) {
        std::ostringstream m;
        m << "flux cubes: the geometry cache was built for a DIFFERENT mesh (cache "
          << c.h.n_faces_core << " faces / geomhash " << c.h.geomhash_core
          << " vs run " << in.n_faces_core << " / " << in.geomhash_core
          << "). Never apply a foreign cache. Regenerate with:\n"
             "    ptm_flux_cube_geom --core <core.ply> --ring <ring.ply> \\\n"
             "        --center-lat " << c.h.center_lat_deg
          << " --center-lon " << c.h.center_lon_deg
          << " --spacing " << c.h.spacing_m << " --radius " << c.h.radius_m
          << " --n-rays " << c.h.n_rays << " --seed " << c.h.seed
          << " --out " << in.cache_file;
        err = m.str();
        return false;
    }

    // --- the configured names must actually be in the cache -----------------
    if (!in.expected_location_names.empty()) {
        std::vector<std::string> missing;
        for (const std::string& want : in.expected_location_names) {
            bool found = false;
            for (const std::string& got : c.names) if (got == want) { found = true; break; }
            if (!found) missing.push_back(want);
        }
        if (!missing.empty()) {
            std::ostringstream m;
            m << "flux cubes: the geometry cache does not contain " << missing.size()
              << " location(s) named in flux_cubes.locations:";
            for (const std::string& s2 : missing) m << " '" << s2 << "'";
            m << ". The cache is built by ptm_flux_cube_geom --locations, so editing the YAML "
                 "does not move a cube. Rebuild the cache with the current location list, or "
                 "remove the entries. Refusing rather than returning a product that silently "
                 "describes the PREVIOUS set of locations.";
            err = m.str();
            return false;
        }
        size_t n_named = 0;
        for (const std::string& got : c.names) if (!got.empty()) ++n_named;
        if (n_named != in.expected_location_names.size())
            std::cout << "  [flux-cubes] NOTE: the cache carries " << n_named
                      << " named location(s) but the config lists "
                      << in.expected_location_names.size()
                      << "; all configured names are present, the extras come from the cache.\n";
    }

    // --- ...and they must be where the config says they are -----------------
    // A name match does not prove the location: moving a point in the YAML
    // without rebuilding the cache keeps the name and silently returns the
    // product for the OLD position. The site-frame offset is recomputed from
    // the configured lat/lon with the geometry tool's own formula
    // (flux_cube_geom.cpp, named-location placement): d = unit(lat/lon),
    // rel = (d - U) * |frame_origin|, x = rel.E, y = rel.N.
    if (!in.expected_locations.empty()) {
        const V3 U0 = mk(c.h.frame_up), N0 = mk(c.h.frame_north), E0 = mk(c.h.frame_east);
        const double r_mean = len3(mk(c.h.frame_origin));
        std::ostringstream bad;
        size_t n_bad = 0, n_unpositioned = 0, n_size_skipped = 0;
        for (const PassInputs::ExpectedLocation& L : in.expected_locations) {
            size_t k = c.names.size();
            for (size_t q = 0; q < c.names.size(); ++q)
                if (c.names[q] == L.name) { k = q; break; }
            if (k == c.names.size()) continue;            // reported by the name check above
            const CubeRecord& cr = c.cubes[k];
            std::ostringstream why;
            if (L.has_lat_lon) {
                const double la = L.lat_deg * kPi / 180.0, lo = L.lon_deg * kPi / 180.0;
                const V3 d = unit3(V3{std::cos(la) * std::cos(lo), std::cos(la) * std::sin(lo),
                                      std::sin(la)});
                const V3 rel = mul(sub(d, U0), r_mean);
                const double x = dot3(rel, E0), y = dot3(rel, N0);
                const double off = std::hypot(x - cr.x_m, y - cr.y_m);
                if (!(off <= kLocationTolM))
                    why << " position: config lat " << std::setprecision(10) << L.lat_deg
                        << " lon " << L.lon_deg << " is site offset (E " << std::setprecision(6)
                        << x << ", N " << y << ") m, the cache has (E " << cr.x_m << ", N "
                        << cr.y_m << ") m, " << off << " m apart (tolerance "
                        << kLocationTolM << " m);";
            } else {
                ++n_unpositioned;
            }
            // A run over caches at SEVERAL heights cannot match one configured
            // height_m in every pass, so the sizes are compared only when the
            // run reads a single cache (compare_location_size); otherwise a
            // configured value is reported as not compared.
            if (L.height_m > 0.0 || L.edge_m > 0.0) {
                if (!in.compare_location_size) {
                    ++n_size_skipped;
                } else {
                    if (L.height_m > 0.0 &&
                        !(std::fabs(L.height_m - cr.height_achieved_m) <= kLocationSizeTolM))
                        why << " height_m: config " << L.height_m << ", cache "
                            << cr.height_achieved_m << ";";
                    const double cube_edge = (cr.edge_m > 0.0f)
                        ? static_cast<double>(cr.edge_m) : c.h.edge_m;
                    if (L.edge_m > 0.0 &&
                        !(std::fabs(L.edge_m - cube_edge) <= kLocationSizeTolM))
                        why << " edge_m: config " << L.edge_m << ", cache " << cube_edge << ";";
                }
            }
            if (!why.str().empty()) {
                ++n_bad;
                bad << "\n    '" << L.name << "':" << why.str();
            }
        }
        if (n_bad > 0) {
            std::ostringstream m;
            m << "flux cubes: " << n_bad << " named location(s) in flux_cubes.locations do "
                 "not match the geometry cache '" << in.cache_file << "':" << bad.str()
              << "\n  The cache fixes where each cube is; editing the YAML does not move it. "
                 "Refusing rather than labeling the cached positions with the configured ones. "
                 "Put the current coordinates in the locations CSV and rebuild the cache:\n"
                 "    ptm_flux_cube_geom --core <core.ply> [--ring <ring.ply>] \\\n"
                 "        --center-lat " << std::setprecision(10) << c.h.center_lat_deg
              << " --center-lon " << c.h.center_lon_deg
              << " --spacing " << c.h.spacing_m << " --radius " << c.h.radius_m;
            if (c.h.height_factor > 0.0) m << " --height-factor " << c.h.height_factor;
            else m << " --height-m " << c.h.height_m;
            if (c.h.edge_factor > 0.0) m << " --edge-factor " << c.h.edge_factor;
            else m << " --edge-m " << c.h.edge_m;
            if (c.h.truncate_rel > 0.0) m << " --truncate-rel " << c.h.truncate_rel;
            m << " --n-rays " << c.h.n_rays << " --seed " << c.h.seed
              << " \\\n        --locations <locations.csv> --out " << in.cache_file
              << "\n  or set the YAML entries back to the coordinates the cache was built with.";
            err = m.str();
            return false;
        }
        if (n_unpositioned > 0)
            std::cout << "  [flux-cubes] NOTE: " << n_unpositioned << " configured location(s) "
                         "give no lat_deg/lon_deg, so only their names were checked against "
                         "the cache.\n";
        if (n_size_skipped > 0)
            std::cout << "  [flux-cubes] NOTE: " << n_size_skipped << " configured location(s) "
                         "give height_m or edge_m, which are NOT compared with the cache because "
                         "this run reads caches at more than one height; their positions were "
                         "compared.\n";
    }

    // Minor 6: the ring identity, not just the core's.
    if (c.h.ring_hash != 0 && in.ring_hash != 0 && c.h.ring_hash != in.ring_hash) {
        std::ostringstream m;
        m << "flux cubes: the cache was built against a DIFFERENT occluder ring (cache hash "
          << c.h.ring_hash << " vs run " << in.ring_hash
          << "). The ring sets every sky fraction in the file, so this is as wrong as a "
             "mismatched core mesh. Rebuild the cache against the run's ring.";
        err = m.str();
        return false;
    }
    if (c.h.ring_hash == 0 || in.ring_hash == 0)
        std::cout << "  [flux-cubes] NOTE: ring identity NOT verified (cache hash "
                  << c.h.ring_hash << ", run hash " << in.ring_hash
                  << "); only the core mesh was checked.\n";

    const size_t nf = in.n_faces_core;
    if (in.emissivity.size() != nf || in.bond_albedo.size() != nf) {
        err = "flux cubes: per-face emissivity/albedo arrays do not match the mesh";
        return false;
    }

    // ------------------------------------------------ time-axis assertions
    const int n_steps = static_cast<int>(in.times->size());
    if (static_cast<int>(in.time_strings->size()) != n_steps) {
        err = "flux cubes: times and time_strings disagree in length";
        return false;
    }
    // ⛔ THE RELATION IS `<=`, NOT `==`, AND THAT IS NOT A LOOSENING.
    // What this check exists to protect is the ABSOLUTE-step indexing below:
    // every series read is at t = equilibration_steps + k for k in [0, n_steps),
    // against arrays strided by total_time_steps (direct_solar) or sized to it
    // (surface_temps, the indirect series, the ephemerides). The condition that
    // keeps every one of those reads in range is
    //     equilibration_steps + n_steps <= total_time_steps,
    // and the `==` spelling was a stricter accident of the era before the
    // checkpoint look-ahead existed. Loosening to `<=` cannot admit a read the
    // `==` form forbade: under `==` the largest absolute index touched is
    // total_time_steps - 1, and under `<=` it is equilibration_steps + n_steps
    // - 1, which is no larger.
    //
    // A CHUNKED run sets thermal.checkpoint_lookahead_steps: 1, and
    // ptm_thermal_sim.cpp:18284 then builds
    //     total_time_steps = equilibration_steps + output_steps + look-ahead,
    // while the exported epoch slice (the `times` handed here) is the contiguous
    // range [equilibration_steps, equilibration_steps + time_steps) -- see the
    // three export paths at ptm_thermal_sim.cpp:13719, 23559 and 23816, each of
    // which slices contiguously and clamps at total_time_steps. So the engine
    // GUARANTEES `<=` and produces `==` only when the look-ahead is zero.
    // Requiring equality refused every chunked window: measured 118 + 59 != 178
    // on a 2-chunk crater80N run. The same defect was fixed in
    // shadow_depth_pass.cpp:137 first; this is the sibling it named.
    //
    // What must still hold, and what the absolute-step indexing below assumes,
    // is that the exported epochs are the FIRST n_steps after equilibration.
    if (in.equilibration_steps + n_steps > in.total_time_steps) {
        std::ostringstream m;
        m << "flux cubes: time-axis relation violated -- equilibration_steps ("
          << in.equilibration_steps << ") + output steps (" << n_steps
          << ") EXCEEDS total_time_steps (" << in.total_time_steps
          << "). The per-epoch rows would index past the end of the series and "
             "silently stop lining up with time_steps.txt.";
        err = m.str();
        return false;
    }
    {
        const int extra = in.total_time_steps - in.equilibration_steps - n_steps;
        if (extra > 0) {
            // Say it out loud. One extra step is the documented chunk
            // look-ahead; more than one is a configuration this pass has never
            // been shown correct for, and it must not pass in silence.
            std::cout << "[flux-cubes] time axis: " << in.equilibration_steps
                      << " equilibration + " << n_steps << " output = "
                      << (in.equilibration_steps + n_steps) << " of "
                      << in.total_time_steps << " integrated steps; the extra "
                      << extra << " step(s) are"
                      << (extra == 1 ? " the chunk checkpoint look-ahead"
                                     : " NOT the single-step chunk look-ahead -- "
                                       "confirm what they are")
                      << " and are integrated but not exported.\n";
            if (extra > 1) {
                std::cerr << "  [flux-cubes] WARNING: " << extra << " integrated steps "
                             "are neither equilibration nor output. The exported epochs "
                             "are assumed to be the FIRST " << n_steps
                          << " after equilibration; if that is not true for this run, "
                             "every cube row is mislabeled in time.\n";
            }
        }
    }
    // ------------------------------------------------- ephemeris coverage
    // Every epoch reads the Sun (and the primary, when there is one) at its
    // absolute step. A series shorter than the window used to read as "no sun"
    // past its end: zero direct flux, no message. The shadow-depth pass refuses
    // the same input (shadow_depth_pass.cpp), and so does this one.
    {
        const size_t need = static_cast<size_t>(in.equilibration_steps + n_steps);
        if (!in.sun_positions || in.sun_positions->size() < need) {
            std::ostringstream m;
            m << "flux cubes: the Sun ephemeris holds "
              << (in.sun_positions ? in.sun_positions->size() : 0)
              << " epoch(s) but this window reads absolute steps up to " << (need - 1)
              << " (equilibration " << in.equilibration_steps << " + output " << n_steps
              << "). The direct component would be zero at every uncovered epoch with no "
                 "error. Refusing.";
            err = m.str();
            return false;
        }
        if (in.primary_positions && !in.primary_positions->empty() &&
            in.primary_positions->size() < need) {
            std::ostringstream m;
            m << "flux cubes: the primary-body ephemeris holds "
              << in.primary_positions->size() << " epoch(s) but this window reads absolute "
                 "steps up to " << (need - 1) << ". primary_refl and primary_ir would be "
                 "zero at every uncovered epoch with no error. Refusing.";
            err = m.str();
            return false;
        }
    }
    if (n_steps > 1) {
        const double dt0 = (*in.times)[1] - (*in.times)[0];
        double spread = 0.0;
        for (int k = 1; k < n_steps; ++k)
            spread = std::max(spread, std::fabs(((*in.times)[k] - (*in.times)[k - 1]) - dt0));
        if (spread > 1e-6)
            std::cerr << "  [flux-cubes] WARNING: output epoch spacing is not uniform; "
                         "measured spread " << spread << " s about dt = " << dt0 << " s.\n";
    }

    // --------------------------------------------- reconstructibility guard
    {
        std::string reason;
        if (!radiation::checkReconstructible(in.emissivity, in.bond_albedo, reason)) {
            err = "flux cubes: the terrain radiosity cannot be reconstructed -- " + reason +
                  ". Refusing rather than dividing by a collapsing denominator.";
            return false;
        }
    }

    // ------------------------------- the union of partner faces, gathered once
    // B is reconstructed only for faces some patch actually looks at. On a
    // 2.4 M-face mesh that is a small fraction, and it turns the per-step cost
    // from "touch every face" into "touch what the product needs".
    std::vector<int32_t> partners(c.cols.begin(), c.cols.end());
    std::sort(partners.begin(), partners.end());
    partners.erase(std::unique(partners.begin(), partners.end()), partners.end());
    std::unordered_map<int32_t, int32_t> slot;
    slot.reserve(partners.size() * 2);
    for (size_t i = 0; i < partners.size(); ++i)
        slot.emplace(partners[i], static_cast<int32_t>(i));
    const size_t n_part = partners.size();

    const size_t n_cubes = c.h.n_cubes;
    const size_t n_patch = n_cubes * kNumFaces;
    std::cout << "\n[flux-cubes] " << n_cubes << " cubes, " << n_patch << " patches, "
              << c.h.total_nnz << " view-factor entries over " << n_part
              << " distinct partner faces; " << n_steps << " output steps.\n";
    std::cout << "[flux-cubes] cache " << in.cache_file << "\n"
              << "[flux-cubes] height " << c.h.height_m << " m, edge " << c.h.edge_m
              << " m (local median facet edge " << c.h.local_median_edge_m << " m, "
              << (c.h.height_factor > 0.0 ? "site-relative" : "ABSOLUTE override") << ")\n";
    if (!in.source_albedo)
        std::cout << "[flux-cubes] NOTE: no incidence-albedo function supplied; the visible "
                     "source albedo falls back to the Bond value (legacy "
                     "behavior). Recorded in the file header.\n";
    if (in.primary_radius_km <= 0.0)
        std::cout << "[flux-cubes] NOTE: primary body treated as a POINT source "
                     "(no radius supplied); mean_cos is the single center direction. "
                     "Recorded in the file header.\n";

    // ---------------------------------------------------------- flux storage
    // [cube][face][component][step], float64: the dynamic range at a polar site
    // runs 0.01 -> 1000 W/m^2 and the PSR end is the end that sizes a radiator.
    std::vector<double> flux(n_cubes * kNumFaces * kNumComponents * n_steps, 0.0);
    std::vector<double> t_face(n_cubes * n_steps, 0.0);
    std::vector<uint8_t> penumbral(n_cubes * n_steps, 0);
    auto FI = [&](size_t cu, int fa, int co, int st) {
        return ((cu * kNumFaces + fa) * kNumComponents + co) * n_steps + st;
    };

    const V3 U = mk(c.h.frame_up), N = mk(c.h.frame_north), E = mk(c.h.frame_east);
    const std::array<V3, kNumFaces> fn = {U, mul(U, -1.0), N, mul(N, -1.0), E, mul(E, -1.0)};

    // ---- PTM-823: the open-mesh body sphere -------------------------------
    // Every sun and primary ray first takes the engine's direct-solar R_occ test
    // (in.body_blocks), then the mesh. On an open patch a ray aimed below the
    // limb that misses the terrain is otherwise counted as clear sky. When the
    // sphere never blocks a ray the mesh misses, every result is unchanged.
    const std::string stem_id = in.output_stem.empty() ? std::string("flux_cubes")
                                                       : in.output_stem;
    using OccFn = std::function<bool(const std::array<double, 3>&,
                                     const std::array<double, 3>&)>;
    std::atomic<size_t> n_body_blocked{0};
    const OccFn occ = in.body_blocks
        ? OccFn([&in, &n_body_blocked](const std::array<double, 3>& o,
                                        const std::array<double, 3>& d) {
              if (in.body_blocks(o, d)) {
                  n_body_blocked.fetch_add(1, std::memory_order_relaxed);
                  return true;
              }
              return in.occluded(o, d);
          })
        : in.occluded;
    std::ostringstream prov;   // appended to <output_dir>/provenance.txt
    prov << "flux_cube_body_cull[" << stem_id << "]="
         << (in.body_blocks ? "on (R_occ = " : "not used: no body sphere (closed mesh or no R_occ)");
    if (in.body_blocks) prov << std::setprecision(10) << in.body_sphere_radius << ")";
    prov << "\n";

    // Pre-flight, per site (named location, else "lattice"): can the sphere block
    // anything the mesh does not? From each cube, the sphere covers every
    // direction more than dep = acos(R_occ/|x|) below the local horizontal (body
    // radial). A ray at exactly that depression, at 72 azimuths, is cast against
    // the MESH ONLY: an azimuth the mesh leaves open there is where a leak can
    // happen. For the worst cube of each site the mesh horizon is also bisected
    // per azimuth, so the record shows the margin, not just a flag.
    if (in.body_blocks && in.body_sphere_radius > 0.0) {
        const double R = in.body_sphere_radius;
        constexpr int kAz = 72;
        auto cross = [](const V3& a, const V3& b) {
            return V3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
        };
        auto dirAt = [&](const V3& up, const V3& e1, const V3& e2, double az, double el) {
            const V3 h = add(mul(e1, std::cos(az)), mul(e2, std::sin(az)));
            return unit3(add(mul(h, std::cos(el)), mul(up, std::sin(el))));
        };
        auto frameOf = [&](const V3& xc, V3& up, V3& e1, V3& e2) {
            up = unit3(xc);
            V3 t = sub(E, mul(up, dot3(E, up)));
            if (len3(t) < 1e-6) t = sub(N, mul(up, dot3(N, up)));
            e1 = unit3(t); e2 = cross(up, e1);
        };
        std::vector<int> n_open(n_cubes, 0);
        std::vector<double> dep(n_cubes, 0.0);
        #pragma omp parallel for schedule(dynamic, 16)
        for (long long cu = 0; cu < static_cast<long long>(n_cubes); ++cu) {
            const V3 xc = mk(c.cubes[cu].center);
            const double r = len3(xc);
            if (!(r > R)) { dep[cu] = -1.0; continue; }   // inside the sphere: reported below
            V3 up, e1, e2; frameOf(xc, up, e1, e2);
            const double d = std::acos(R / r);
            dep[cu] = d;
            int open = 0;
            for (int a = 0; a < kAz; ++a)
                if (!in.occluded(arr(xc), arr(dirAt(up, e1, e2, 2.0 * M_PI * a / kAz, -d)))) ++open;
            n_open[cu] = open;
        }
        std::vector<std::string> gname;
        std::unordered_map<std::string, size_t> gidx;
        std::vector<size_t> g_n, g_open, g_inside, g_rep;
        std::vector<double> g_dmin, g_dmax;
        for (size_t cu = 0; cu < n_cubes; ++cu) {
            const std::string nm = (cu < c.names.size() && !c.names[cu].empty()) ? c.names[cu]
                                                                                 : std::string("lattice");
            auto it = gidx.find(nm);
            size_t g;
            if (it == gidx.end()) {
                g = gname.size(); gidx.emplace(nm, g); gname.push_back(nm);
                g_n.push_back(0); g_open.push_back(0); g_inside.push_back(0); g_rep.push_back(cu);
                g_dmin.push_back(1e300); g_dmax.push_back(-1e300);
            } else {
                g = it->second;
            }
            ++g_n[g];
            if (dep[cu] < 0.0) { ++g_inside[g]; continue; }
            if (n_open[cu] > 0) ++g_open[g];
            g_dmin[g] = std::min(g_dmin[g], dep[cu]); g_dmax[g] = std::max(g_dmax[g], dep[cu]);
            if (n_open[cu] > n_open[g_rep[g]] || dep[g_rep[g]] < 0.0) g_rep[g] = cu;
        }
        const double deg = 180.0 / M_PI;
        size_t total_open = 0;
        for (size_t g = 0; g < gname.size(); ++g) {
            total_open += g_open[g];
            std::ostringstream line;
            line << std::fixed << std::setprecision(2);
            line << g_n[g] << " cubes";
            if (g_inside[g]) line << " (" << g_inside[g] << " at or inside the sphere: not scanned)";
            if (g_n[g] > g_inside[g]) {
                const size_t rc = g_rep[g];
                const V3 xc = mk(c.cubes[rc].center);
                V3 up, e1, e2; frameOf(xc, up, e1, e2);
                const double d = dep[rc];
                double hmin = 1e300, hmin_az = 0.0;
                for (int a = 0; a < kAz; ++a) {
                    const double az = 2.0 * M_PI * a / kAz;
                    double h;
                    if (!in.occluded(arr(xc), arr(dirAt(up, e1, e2, az, -d)))) {
                        h = -d;                                 // open at the cap edge: horizon below it
                    } else {
                        double lo = -d, hi = 60.0 / deg;        // lo blocked
                        if (in.occluded(arr(xc), arr(dirAt(up, e1, e2, az, hi)))) { h = hi; }
                        else {
                            for (int it = 0; it < 20; ++it) {
                                const double mid = 0.5 * (lo + hi);
                                if (in.occluded(arr(xc), arr(dirAt(up, e1, e2, az, mid)))) lo = mid; else hi = mid;
                            }
                            h = hi;
                        }
                    }
                    if (h < hmin) { hmin = h; hmin_az = az; }
                }
                line << "; cap depression " << g_dmin[g] * deg << ".." << g_dmax[g] * deg
                     << " deg; cubes with a mesh-open azimuth at the cap edge (" << kAz
                     << " azimuths): " << g_open[g] << "; worst cube min mesh horizon "
                     << std::showpos << hmin * deg << std::noshowpos << " deg (azimuth "
                     << hmin_az * deg << " deg from east) vs its cap edge "
                     << std::showpos << -d * deg << std::noshowpos << " deg";
                line << (g_open[g] == 0 ? " -> no leak band found: the sphere should block nothing the mesh misses"
                                        : " -> LEAK BAND: the sphere culls rays the mesh leaves open");
            }
            std::cout << "[flux-cubes] body sphere (PTM-823) " << stem_id << " site '" << gname[g]
                      << "': " << line.str() << "\n";
            prov << "flux_cube_horizon_vs_cap[" << stem_id << "][" << gname[g] << "]="
                 << line.str() << "\n";
        }
        if (total_open > 0)
            std::cerr << "  [flux-cubes] NOTE (PTM-823): " << total_open << " cube(s) see the "
                         "body-sphere cap through a gap in the mesh horizon; their sun/primary "
                         "rays there are now blocked by the sphere (before PTM-823 they counted "
                         "as clear sky).\n";
    }
    // PER-CUBE edge: a named location may override the cache-level value, so the
    // shadow-probe spread must come from the cube, not the header. Falls back to
    // the header for caches written before the field existed (it was a pad, so
    // an old cache reads 0 there).
    const double header_half_edge = 0.5 * c.h.edge_m;
    // Disk sampling and the distance handed to the engine's primary-body
    // functions are done in KM, because that is the unit those functions and
    // primary_positions use. Ray DIRECTIONS are dimensionless, so they are
    // valid in the mesh frame without conversion.
    const double to_km = in.mesh_to_km;
    size_t n_primary_lit = 0;

    std::vector<double> B_vis(n_part, 0.0), B_ir(n_part, 0.0);

    for (int k = 0; k < n_steps; ++k) {
        const int t = in.equilibration_steps + k;
        const double S_t = in.solar_const(t);

        // ---- reconstruct the terrain radiosity for the partner faces
        #pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(n_part); ++i) {
            const size_t j = static_cast<size_t>(partners[i]);
            const double eps = in.emissivity[j];
            const double A = in.bond_albedo[j];
            const double T = (*in.surface_temps)[j][t];
            const double q_dir = (*in.direct_solar)[j * in.total_time_steps + t];
            const double q_ir = (*in.indirect_ir_series)[j][t];
            const double q_vis = (*in.indirect_solar_series)[j][t];
            const double a_src = in.source_albedo ? in.source_albedo(A, q_dir, S_t) : A;
            B_ir[i] = radiation::reconstructBIR(eps, T, q_ir);
            B_vis[i] = radiation::reconstructBVis(a_src, q_dir, A, q_vis);
        }

        // ---- sun direction for this epoch, same normalization as the sweep
        V3 s_hat{0, 0, 0};
        bool have_sun = false;
        if (in.sun_positions && t < static_cast<int>(in.sun_positions->size())) {
            const auto& sp = (*in.sun_positions)[t];
            const V3 v{sp[0], sp[1], sp[2]};
            if (len3(v) > 0) { s_hat = unit3(v); have_sun = true; }
        }
        // Solar angular radius this epoch, from the actual Sun distance.
        double solar_alpha = 0.0;
        if (have_sun && in.solar_radius_km > 0.0 && in.sun_positions
            && t < static_cast<int>(in.sun_positions->size())) {
            const auto& sp = (*in.sun_positions)[t];
            const double dist_km = std::sqrt(sp[0]*sp[0] + sp[1]*sp[1] + sp[2]*sp[2]);
            if (dist_km > in.solar_radius_km)
                solar_alpha = std::asin(in.solar_radius_km / dist_km);
        }
        V3 p_pos{0, 0, 0};
        bool have_primary = false;
        if (in.primary_positions && t < static_cast<int>(in.primary_positions->size())) {
            const auto& pp = (*in.primary_positions)[t];
            p_pos = V3{pp[0], pp[1], pp[2]};
            if (len3(p_pos) > 0) have_primary = true;
        }

        #pragma omp parallel for schedule(dynamic, 8)
        for (long long cu = 0; cu < static_cast<long long>(n_cubes); ++cu) {
            const CubeRecord& cr = c.cubes[cu];
            const V3 xc = mk(cr.center);
            t_face[static_cast<size_t>(cu) * n_steps + k] =
                (*in.surface_temps)[cr.face_id][t];

            // ---- direct solar: ONE probe per CUBE, not per face.
            // Terrain occlusion is a property of the POINT; the normal enters
            // only through max(0, n.s). A per-face shadow ray could return a
            // rigid cube with one face lit and the opposite face shadowed --
            // an answer with no physical referent.
            double v_sun = 0.0;
            if (have_sun) {
                // ANGULAR (finite solar disk) at the cube center, plus the
                // SPATIAL spread over the cube's own extent. The two smear the
                // shadow edge on comparable scales at this height and neither
                // subsumes the other: the disk is ~0.045 m at a 10 m occluder,
                // the cube is edge_m ~0.11 m. Averaging the disk fraction over
                // the spatial probe points captures both.
                const double he = (cr.edge_m > 0.0f) ? 0.5 * cr.edge_m : header_half_edge;
                std::array<V3, 5> pts = {xc, xc, xc, xc, xc};
                int n_pt = 1;
                if (in.shadow_probe_rays >= 5) {
                    // the physics review's falsifying test (PTM_FLUXCUBE_PROBE_AXIS=U): the
                    // horizontal probe does NOT sample the cube's VERTICAL
                    // extent, and at grazing sun -- winter is entirely
                    // -4.3..+7.3 deg -- the unsampled +-edge/2 of HEIGHT maps to
                    // 8-19x the horizontal extent along the shadow direction.
                    // This arm displaces the same five points along U instead,
                    // so the diff against the shipped product measures exactly
                    // the sensitivity the default probe cannot see.
                    static const bool probe_u = [] {
                        const char* e = std::getenv("PTM_FLUXCUBE_PROBE_AXIS");
                        return e && (e[0] == 'U' || e[0] == 'u');
                    }();
                    if (probe_u) {
                        pts[1] = add(xc, mul(U, he)); pts[2] = sub(xc, mul(U, he));
                        pts[3] = add(xc, mul(E, he)); pts[4] = sub(xc, mul(E, he));
                    } else {
                        pts[1] = add(xc, mul(E, he)); pts[2] = sub(xc, mul(E, he));
                        pts[3] = add(xc, mul(N, he)); pts[4] = sub(xc, mul(N, he));
                    }
                    n_pt = 5;
                }
                double acc = 0.0; bool any_split = false;
                for (int q = 0; q < n_pt; ++q) {
                    bool sp = false;
                    acc += solarDiskVisibility(pts[q], s_hat, solar_alpha,
                                               in.solar_disk_samples, occ, sp);
                    any_split = any_split || sp;
                }
                v_sun = acc / n_pt;
                // Penumbral now means what it says: the disk was partially
                // occluded, or the spatial probes disagreed.
                if (any_split || (v_sun > 0.0 && v_sun < 1.0))
                    penumbral[static_cast<size_t>(cu) * n_steps + k] = 1;
            }

            // ---- primary body: sample the disk ONCE per cube, then every face
            // reuses the clear directions. calculatePrimary*() take mean_cos,
            // not a normal, which is exactly this decomposition.
            std::array<double, kNumFaces> mean_cos{};
            mean_cos.fill(0.0);
            double p_dist_km = 0.0;
            V3 p_hat{0, 0, 0};
            if (have_primary) {
                const V3 xc_km = mul(xc, to_km);          // cube center in KM
                std::vector<V3> dirs;
                sampleDisk(xc_km, p_pos, in.primary_radius_km,
                           in.primary_disk_samples, dirs);
                p_dist_km = len3(sub(p_pos, xc_km));
                p_hat = unit3(sub(p_pos, xc_km));
                const double invK = 1.0 / static_cast<double>(dirs.size());
                for (const V3& d : dirs) {
                    // Directions are dimensionless: cast in the MESH frame.
                    if (occ(arr(xc), arr(d))) continue;
                    for (int fa = 0; fa < kNumFaces; ++fa)
                        mean_cos[fa] += std::max(0.0, dot3(fn[fa], d)) * invK;
                }
                double mc_sum = 0.0;
                for (int fa = 0; fa < kNumFaces; ++fa) mc_sum += mean_cos[fa];
                if (mc_sum > 0.0) {
                    #pragma omp atomic
                    ++n_primary_lit;
                }
            }

            for (int fa = 0; fa < kNumFaces; ++fa) {
                const PatchRecord& pr = c.patches[static_cast<size_t>(cu) * kNumFaces + fa];

                const double q_dir = (have_sun && v_sun > 0.0)
                    ? v_sun * S_t * std::max(0.0, dot3(fn[fa], s_hat)) : 0.0;

                double q_tvis = 0.0, q_tir = 0.0;
                for (uint32_t e = 0; e < pr.row_len; ++e) {
                    const size_t idx = pr.row_offset + e;
                    const auto it = slot.find(c.cols[idx]);
                    const double F = c.vals[idx];
                    q_tvis += F * B_vis[it->second];
                    q_tir += F * B_ir[it->second];
                }

                double q_pref = 0.0, q_pir = 0.0;
                if (have_primary && mean_cos[fa] > 0.0) {
                    if (in.primary_heating)
                        q_pir = in.primary_heating(mean_cos[fa], p_dist_km);
                    if (in.primary_reflected && in.sun_positions &&
                        t < static_cast<int>(in.sun_positions->size()))
                        q_pref = in.primary_reflected(mean_cos[fa], arr(p_hat),
                                                      (*in.sun_positions)[t],
                                                      (*in.primary_positions)[t], p_dist_km);
                }

                flux[FI(cu, fa, kCompDirect, k)] = q_dir;
                flux[FI(cu, fa, kCompTerrainVis, k)] = q_tvis;
                flux[FI(cu, fa, kCompTerrainIR, k)] = q_tir;
                flux[FI(cu, fa, kCompPrimaryRefl, k)] = q_pref;
                flux[FI(cu, fa, kCompPrimaryIR, k)] = q_pir;
                flux[FI(cu, fa, kCompTotal, k)] = q_dir + q_tvis + q_tir + q_pref + q_pir;
            }
        }
        if ((k % 10) == 0 || k == n_steps - 1)
            std::cout << "\r[flux-cubes] step " << (k + 1) << "/" << n_steps << std::flush;
    }
    std::cout << "\n";

    // Make a silent zero impossible. If the primary body is present in the
    // ephemeris but NO cube ever sees it, that is either a real geometric fact
    // (the primary is permanently below every horizon) or a units/frame slip --
    // and the two must not look alike. At a lunar PSR the primary's IR is
    // 12.8-38.8% of the floor flux, so a quiet zero here is a large error.
    if (in.primary_positions && !in.primary_positions->empty() && n_primary_lit == 0) {
        std::cerr << "  [flux-cubes] WARNING: the primary body is in the ephemeris but NO "
                     "cube saw it at ANY epoch, so primary_refl and primary_ir are "
                     "identically zero in this product. Confirm that is geometry and not a "
                     "units/frame mismatch (mesh_to_km = " << to_km << ").\n";
    } else if (in.primary_positions && !in.primary_positions->empty()) {
        std::cout << "[flux-cubes] primary body visible at " << n_primary_lit << "/"
                  << (n_cubes * static_cast<size_t>(n_steps)) << " cube-epochs\n";
    }

    // ----------------------------------------------------------------- write
    const std::string dir = in.output_dir + "/flux_cubes";
    const std::string stem = in.output_stem.empty() ? std::string("flux_cubes")
                                                    : in.output_stem;
    try { std::filesystem::create_directories(dir); }
    catch (const std::exception& e) { err = std::string("flux cubes: cannot create output dir: ") + e.what(); return false; }

    // Site-level f_ring summary: the occluder ring carries NO thermal solution,
    // so terrain IR here is CORE-ONLY and the ring's contribution is absent.
    // Reported, never quietly folded into sky.
    double fring_max = 0.0, fring_sum = 0.0;
    for (const auto& p : c.patches) { fring_max = std::max<double>(fring_max, p.f_ring); fring_sum += p.f_ring; }
    const double fring_mean = n_patch ? fring_sum / n_patch : 0.0;
    const char* kRingNote =
        "terrain_ir and terrain_vis EXCLUDE the occluder-ring contribution: the ring is a "
        "shadow-only body with no thermal solution and therefore no radiosity. The omitted "
        "irradiance on a patch is bounded by f_ring * sigma*T_ring^4 for any assumed ring "
        "temperature; per-patch f_ring is in this file.";
    std::cout << "[flux-cubes] f_ring over all patches: mean " << fring_mean
              << ", max " << fring_max << " -- " << kRingNote << "\n";

    std::ostringstream js;
    js << "{\"ptm3d_flux_cubes_schema\":1,"
       << "\"units\":\"W m-2 incident (no albedo/emissivity applied)\","
       << "\"face_order\":[\"zenith\",\"nadir\",\"north\",\"south\",\"east\",\"west\"],"
       << "\"component_order\":[\"direct\",\"terrain_vis\",\"terrain_ir\",\"primary_refl\","
          "\"primary_ir\",\"total\"],"
       << "\"terrain_ir_scope\":\"core_mesh_only\","
       << "\"ring_note\":\"" << kRingNote << "\","
       << "\"f_ring_mean\":" << fring_mean << ",\"f_ring_max\":" << fring_max << ","
       << "\"height_m\":" << c.h.height_m << ",\"edge_m\":" << c.h.edge_m << ","
       << "\"height_factor\":" << c.h.height_factor
       << ",\"local_median_edge_m\":" << c.h.local_median_edge_m << ","
       << "\"n_rays\":" << c.h.n_rays << ",\"sampler_id\":" << c.h.sampler_id
       << ",\"seed\":" << c.h.seed << ","
       << "\"frame_max_tilt_deg\":" << c.h.frame_max_tilt_deg << ","
       << "\"equilibration_steps\":" << in.equilibration_steps
       << ",\"total_time_steps\":" << in.total_time_steps
       << ",\"n_output_steps\":" << n_steps << ","
       << "\"time_axis_source\":\"time_steps.txt\","
       << "\"source_albedo\":\"" << (in.source_albedo ? "ptm106_incidence" : "bond_fallback")
       << "\",\"primary_model\":\""
       << (in.primary_radius_km > 0.0 ? "finite_disk" : "point_source")
       << "\",\"primary_disk_samples\":"
       << (in.primary_radius_km > 0.0 ? in.primary_disk_samples : 1) << ","
       << "\"shadow_probe_rays\":" << (in.shadow_probe_rays >= 5 ? 5 : 1) << ","
       << "\"solar_model\":\"" << (in.solar_radius_km > 0.0 ? "finite_disk" : "point_sun")
       << "\",\"solar_disk_samples\":"
       << (in.solar_radius_km > 0.0 ? in.solar_disk_samples : 1) << ","
       << "\"geometry_cache\":\"" << in.cache_file << "\","
       << "\"provenance\":" << (in.provenance_json.empty() ? "null" : in.provenance_json)
       << "}";
    const std::string json = js.str();

    if (in.write_binary) {
        FcbHeader h{};
        std::memcpy(h.magic, "FCB1", 4);
        h.version = 1;
        h.n_cubes = static_cast<uint32_t>(n_cubes);
        h.n_faces_per_cube = kNumFaces;
        h.n_components = kNumComponents;
        h.n_steps = static_cast<uint32_t>(n_steps);
        h.equilibration_steps = in.equilibration_steps;
        h.total_time_steps = in.total_time_steps;
        h.height_m = c.h.height_m;
        h.edge_m = c.h.edge_m;
        h.spacing_m = c.h.spacing_m;
        h.radius_m = c.h.radius_m;
        std::memcpy(h.frame_origin, c.h.frame_origin, sizeof(h.frame_origin));
        std::memcpy(h.frame_east, c.h.frame_east, sizeof(h.frame_east));
        std::memcpy(h.frame_north, c.h.frame_north, sizeof(h.frame_north));
        std::memcpy(h.frame_up, c.h.frame_up, sizeof(h.frame_up));
        h.json_len = static_cast<uint32_t>(json.size());
        h.has_names = c.names.empty() ? 0u : 1u;

        const std::string fp = dir + "/" + stem + ".fcb";
        std::ofstream o(fp, std::ios::binary);
        if (!o) { err = "flux cubes: cannot write " + fp; return false; }
        o.write(reinterpret_cast<const char*>(&h), sizeof(h));
        o.write(reinterpret_cast<const char*>(flux.data()),
                static_cast<std::streamsize>(flux.size() * sizeof(double)));
        o.write(reinterpret_cast<const char*>(t_face.data()),
                static_cast<std::streamsize>(t_face.size() * sizeof(double)));
        o.write(reinterpret_cast<const char*>(penumbral.data()),
                static_cast<std::streamsize>(penumbral.size()));
        std::vector<float> fsky(n_patch), fring(n_patch), fmax(n_patch);
        for (size_t p = 0; p < n_patch; ++p) {
            fsky[p] = c.patches[p].f_sky;
            fring[p] = c.patches[p].f_ring;
            fmax[p] = c.patches[p].f_max;
        }
        o.write(reinterpret_cast<const char*>(fsky.data()), static_cast<std::streamsize>(n_patch * sizeof(float)));
        o.write(reinterpret_cast<const char*>(fring.data()), static_cast<std::streamsize>(n_patch * sizeof(float)));
        o.write(reinterpret_cast<const char*>(fmax.data()), static_cast<std::streamsize>(n_patch * sizeof(float)));
        o.write(reinterpret_cast<const char*>(c.cubes.data()),
                static_cast<std::streamsize>(c.cubes.size() * sizeof(CubeRecord)));
        o.write(reinterpret_cast<const char*>(in.times->data()),
                static_cast<std::streamsize>(n_steps * sizeof(double)));
        for (int k = 0; k < n_steps; ++k) {
            char buf[24];
            std::memset(buf, 0, sizeof(buf));
            std::strncpy(buf, (*in.time_strings)[k].c_str(), sizeof(buf) - 1);
            o.write(buf, sizeof(buf));
        }
        o.write(json.data(), static_cast<std::streamsize>(json.size()));
        // Names LAST, after the JSON, exactly as in the .fcg: a reader written
        // against the original .fcb layout parses everything before this at
        // unchanged offsets and ignores the tail.
        if (!c.names.empty()) {
            std::vector<char> nb(static_cast<size_t>(n_cubes) * kNameLen, 0);
            for (size_t k = 0; k < n_cubes && k < c.names.size(); ++k)
                std::strncpy(&nb[k * kNameLen], c.names[k].c_str(), kNameLen - 1);
            o.write(nb.data(), static_cast<std::streamsize>(nb.size()));
        }
        o.close();
        std::cout << "[flux-cubes] wrote " << fp << " ("
                  << (flux.size() * 8 / 1048576.0) << " MB of flux)\n";
    }

    if (in.write_csv) {
        const std::string cdir = dir + ((stem == "flux_cubes") ? std::string("/csv")
                                                              : ("/csv_" + stem));
        try { std::filesystem::create_directories(cdir); }
        catch (const std::exception& e) { err = std::string("flux cubes: cannot create csv dir: ") + e.what(); return false; }
        // A cap on the lattice fan-out must never drop a NAMED location: those
        // were asked for explicitly by name, and silently omitting one would be
        // the same class of failure as silently dropping a cube at placement.
        size_t n_write = (in.write_csv_max_cubes > 0)
            ? std::min<size_t>(n_cubes, static_cast<size_t>(in.write_csv_max_cubes)) : n_cubes;
        for (size_t k = 0; k < n_cubes; ++k)
            if (k < c.names.size() && !c.names[k].empty()) n_write = n_cubes;
        for (size_t cu = 0; cu < n_write; ++cu) {
            char nm[96];
            // Outputs are keyed by NAME for named locations -- that is the point
            // of naming them -- and by index for lattice cubes.
            if (cu < c.names.size() && !c.names[cu].empty())
                std::snprintf(nm, sizeof(nm), "/cube_%s.csv", c.names[cu].c_str());
            else
                std::snprintf(nm, sizeof(nm), "/cube_%05zu.csv", cu);
            std::ofstream o(cdir + nm);
            if (!o) { err = "flux cubes: cannot write " + cdir + nm; return false; }
            const CubeRecord& cr = c.cubes[cu];
            const bool is_named = (cu < c.names.size() && !c.names[cu].empty());
            o << "# ptm3d_flux_cubes v1  cube=" << cu;
            if (is_named) o << "  name=" << c.names[cu];
            o << (is_named ? "  placement=named_location"
                           : "  placement=lattice")
              << "  grid_ij=(" << cr.grid_i << "," << cr.grid_j << ")"
              << "  x_m=" << cr.x_m << " y_m=" << cr.y_m << "\n"
              << "# face_id=" << cr.face_id << " slope_deg=" << cr.slope_deg
              << " aspect_deg=" << cr.aspect_deg << " area_m2=" << cr.area_m2
              << " elev_m=" << cr.elev_m << " height_m=" << cr.height_achieved_m
              // PER-CUBE edge, not the cache-level one: a named location may
              // override it, and printing the header value here reported 2.71 m
              // for a cube that was actually built with 1.00 m.
              << " edge_m=" << ((cr.edge_m > 0.0f) ? static_cast<double>(cr.edge_m) : c.h.edge_m)
              << "\n"
              << "# frame: site-local ENU at " << c.h.center_lat_deg << ", "
              << c.h.center_lon_deg << "; up=+U north=+N east=+E, right-handed\n"
              << "# units: W/m^2 INCIDENT (no albedo/emissivity applied); T in K\n"
              << "# " << kRingNote << "\n";
            for (int fa = 0; fa < kNumFaces; ++fa) {
                const PatchRecord& pr = c.patches[cu * kNumFaces + fa];
                o << "# " << faceName(fa) << ": f_sky=" << pr.f_sky
                  << " f_ring=" << pr.f_ring << " f_max=" << pr.f_max
                  << " n_faces_50pct=" << pr.n_faces_50pct << "\n";
            }
            o << "utc,et_s,T_face_K,penumbral";
            for (int fa = 0; fa < kNumFaces; ++fa)
                for (int co = 0; co < kNumComponents; ++co)
                    o << "," << faceName(fa) << "_" << componentName(co);
            o << "\n";
            o.setf(std::ios::scientific);
            o.precision(9);
            for (int k = 0; k < n_steps; ++k) {
                o << (*in.time_strings)[k] << "," << (*in.times)[k] << ","
                  << t_face[cu * n_steps + k] << ","
                  << static_cast<int>(penumbral[cu * n_steps + k]);
                for (int fa = 0; fa < kNumFaces; ++fa)
                    for (int co = 0; co < kNumComponents; ++co)
                        o << "," << flux[FI(cu, fa, co, k)];
                o << "\n";
            }
        }
        std::cout << "[flux-cubes] wrote " << n_write << " per-cube CSV files to " << cdir
                  << (n_write < n_cubes ? "  (capped by write_csv_max_cubes)" : "") << "\n";
    }

    // PTM-823 record: what the body sphere did in this pass, next to the run's
    // other provenance. Kept OUT of the .fcb so a pass the sphere does not touch
    // writes the same product file as before.
    if (in.body_blocks) {
        prov << "flux_cube_body_cull_effect[" << stem_id << "]=" << n_body_blocked.load()
             << " cube ray(s) blocked by the body sphere\n";
        std::cout << "[flux-cubes] body sphere (PTM-823) " << stem_id << ": "
                  << n_body_blocked.load() << " cube ray(s) blocked by the body sphere\n";
    }
    {
        const std::string pp = in.output_dir + "/provenance.txt";
        std::ofstream pv(pp, std::ios::app);
        if (pv) pv << prov.str();
        else std::cerr << "  [flux-cubes] WARNING: cannot append the PTM-823 body-sphere record to "
                       << pp << "; this run's directory will not say whether the sphere was "
                          "applied to the cube rays.\n";
    }

    return true;
}

}  // namespace fluxcubes
}  // namespace ptm
