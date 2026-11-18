// ptm_merge_chunks.cpp
//
// Merges per-chunk statistics from a chunked PTM3D run (produced by
// scripts/run_chunked_simulation.py) into a single set of final
// statistics, mathematically equivalent to a hypothetical single-shot
// run over the concatenated time window.
//
// Each input chunk dir is expected to contain a statistics/_raw/
// subdirectory with:
//   <datatype>_sum.bin     uint32 N + N float64
//   <datatype>_count.bin   uint32 N + N uint32
//   ice_depth_sum_E.bin    uint32 n_layers + uint32 n_faces + n_layers*n_faces float32
//   ice_depth_n_samples.txt   single integer line
//   depth_grid.bin         uint32 n_layers + n_layers float64
// plus the existing statistics/<datatype>_max.bin and _min.bin files.
//
// Output: a target dir <outdir>/statistics/<datatype>_{max,avg,min}.bin
// + the ice products from the merged accumulator: ice_depth.bin and its
// companions (retention criterion, ice/ice_retention.h, the same code the
// engine runs) and, when requested, ice_depth_flux.bin (flux threshold).
//
// CLI:
//   ptm_merge_chunks --chunks DIR1 DIR2 ... --output OUTDIR [ice options]
//
// The ice criterion is read from each chunk's statistics/_raw/ice_depth_params.txt
// (written by the engine), so the merged products use the settings the run was
// configured with. The ice options override them, with a warning: that is how
// a finished run is re-evaluated at another age, inventory or threshold from
// the same saved sums (--ice-only writes the ice products alone).

#include "../ice/ice_table.h"
#include "../ice/ice_retention.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ---- I/O helpers ----------------------------------------------------------

bool read_uint32(std::ifstream& is, uint32_t& v) {
    is.read(reinterpret_cast<char*>(&v), sizeof(uint32_t));
    return is.good();
}

bool read_float32_array(const fs::path& p, uint32_t& N_out,
                        std::vector<float>& data) {
    std::ifstream is(p, std::ios::binary);
    if (!is) { std::cerr << "ERROR opening " << p << "\n"; return false; }
    if (!read_uint32(is, N_out)) { std::cerr << "ERROR reading N from " << p << "\n"; return false; }
    data.resize(N_out);
    is.read(reinterpret_cast<char*>(data.data()), N_out * sizeof(float));
    return is.good();
}

bool read_double_array(const fs::path& p, uint32_t& N_out,
                       std::vector<double>& data) {
    std::ifstream is(p, std::ios::binary);
    if (!is) { std::cerr << "ERROR opening " << p << "\n"; return false; }
    if (!read_uint32(is, N_out)) { std::cerr << "ERROR reading N from " << p << "\n"; return false; }
    data.resize(N_out);
    is.read(reinterpret_cast<char*>(data.data()), N_out * sizeof(double));
    return is.good();
}

bool read_uint32_array(const fs::path& p, uint32_t& N_out,
                       std::vector<uint32_t>& data) {
    std::ifstream is(p, std::ios::binary);
    if (!is) { std::cerr << "ERROR opening " << p << "\n"; return false; }
    if (!read_uint32(is, N_out)) { std::cerr << "ERROR reading N from " << p << "\n"; return false; }
    data.resize(N_out);
    is.read(reinterpret_cast<char*>(data.data()), N_out * sizeof(uint32_t));
    return is.good();
}

bool write_float32_array(const fs::path& p, const std::vector<float>& data) {
    std::ofstream os(p, std::ios::binary);
    if (!os) return false;
    uint32_t N = static_cast<uint32_t>(data.size());
    os.write(reinterpret_cast<const char*>(&N), sizeof(uint32_t));
    os.write(reinterpret_cast<const char*>(data.data()),
             data.size() * sizeof(float));
    return os.good();
}

// ---- Per-datatype merge ---------------------------------------------------

// Combine max/min/sum/count across N chunks for one datatype.
// Writes <outdir>/statistics/<name>_{max,avg,min}.bin.
//
// `strict` changes what a MISSING component means.
//   false (the five fixed datatypes): all four side files or nothing. A
//         thermal-only run has no illumination statistics at all, and skipping
//         that quietly is the documented behavior.
//   true  (subsurface depth maps): the components present are whatever
//         output.map_depth_stats asked for, so "avg is absent" is normal and
//         "avg is present in chunk 0 and absent in chunk 3" is a defect. Each
//         component is merged independently and any disagreement between chunks
//         is an ERROR, because the alternative is a merged product that quietly
//         drops a map the chunks all carried.
//
// `nan_when_empty`: what a face with NO contribution in any chunk gets. For the
// surface fields that is 0.0 and always has been (a face with no samples had no
// flux). For SHADOW DEPTH it must be NaN, because 0.0 is a MEANING in that field
// -- "the Sun is visible at the surface" -- so scrubbing an empty face to 0.0
// would turn "this face was never resolved" into "this face was always lit", at
// every permanently shadowed face, silently. NaN also matches what the engine
// writes for the same case in a single-chunk run, so a merged window and an
// unchunked one say the same thing.
bool merge_one_datatype(const std::string& name,
                        const std::vector<fs::path>& chunk_dirs,
                        const fs::path& outdir,
                        uint32_t expected_N,
                        bool strict = false,
                        bool nan_when_empty = false) {
    std::vector<float>   merged_max, merged_min;
    std::vector<double>  merged_sum;
    std::vector<uint64_t> merged_count;
    merged_max.assign(expected_N, std::numeric_limits<float>::lowest());
    merged_min.assign(expected_N, std::numeric_limits<float>::max());
    merged_sum.assign(expected_N, 0.0);
    merged_count.assign(expected_N, 0u);

    auto components = [&](const fs::path& cd, bool& has_avg, bool& has_max,
                          bool& has_min) {
        const fs::path stats = cd / "statistics";
        const fs::path raw   = stats / "_raw";
        has_max = fs::exists(stats / (name + "_max.bin"));
        has_min = fs::exists(stats / (name + "_min.bin"));
        has_avg = fs::exists(raw / (name + "_sum.bin")) &&
                  fs::exists(raw / (name + "_count.bin"));
    };

    bool want_avg = false, want_max = false, want_min = false;
    components(chunk_dirs.front(), want_avg, want_max, want_min);
    if (strict && !(want_avg || want_max || want_min)) {
        std::cerr << "[merge] ERROR: " << name << " has no mergeable component in "
                  << chunk_dirs.front() << "\n";
        return false;
    }
    // An `_avg.bin` MAP with no `_raw` sum/count behind it is a chunked run
    // that was not given `output.write_raw_stats: true`. Averages cannot be
    // merged from the per-chunk averages alone (the chunks need not carry equal
    // sample counts), so the honest answer is to refuse -- not to merge max and
    // min, drop the average, and report success.
    if (strict && !want_avg &&
        fs::exists(chunk_dirs.front() / "statistics" / (name + "_avg.bin"))) {
        std::cerr << "[merge] ERROR: " << name << " has a per-chunk _avg.bin but no "
                     "statistics/_raw/" << name << "_{sum,count}.bin to merge it "
                     "from. Re-run the chunks with output.write_raw_stats: true; an "
                     "average cannot be recovered from per-chunk averages alone.\n";
        return false;
    }

    bool any_loaded = false;
    for (const auto& cd : chunk_dirs) {
        const fs::path stats = cd / "statistics";
        const fs::path raw   = stats / "_raw";
        const fs::path p_max   = stats / (name + "_max.bin");
        const fs::path p_min   = stats / (name + "_min.bin");
        const fs::path p_sum   = raw   / (name + "_sum.bin");
        const fs::path p_count = raw   / (name + "_count.bin");

        bool has_avg, has_max, has_min;
        components(cd, has_avg, has_max, has_min);
        if (!strict) {
            // If any of the four files is missing in any chunk, the datatype
            // wasn't computed (e.g. illumination on a thermal-only run);
            // skip this datatype for the whole merge.
            if (!has_max || !has_min || !has_avg) {
                return true;  // not an error — just nothing to merge
            }
        } else if (has_avg != want_avg || has_max != want_max ||
                   has_min != want_min) {
            std::cerr << "[merge] ERROR: " << name << " carries a different set of "
                         "statistics in " << cd << " than in "
                      << chunk_dirs.front() << " (avg/max/min present: "
                      << has_avg << has_max << has_min << " vs "
                      << want_avg << want_max << want_min
                      << "). The chunks were not run with the same "
                         "output.map_depth_stats; refusing to merge.\n";
            return false;
        }

        uint32_t N;
        std::vector<float>   max_arr, min_arr;
        std::vector<double>  sum_arr;
        std::vector<uint32_t> count_arr;
        if (want_max && (!read_float32_array(p_max,  N, max_arr) || N != expected_N)) return false;
        if (want_min && (!read_float32_array(p_min,  N, min_arr) || N != expected_N)) return false;
        if (want_avg && (!read_double_array (p_sum,  N, sum_arr) || N != expected_N)) return false;
        if (want_avg && (!read_uint32_array (p_count,N, count_arr) || N != expected_N)) return false;

        for (uint32_t i = 0; i < expected_N; ++i) {
            if (want_max && max_arr[i] > merged_max[i]) merged_max[i] = max_arr[i];
            if (want_min && min_arr[i] < merged_min[i]) merged_min[i] = min_arr[i];
            if (want_avg) {
                merged_sum[i]   += sum_arr[i];
                merged_count[i] += count_arr[i];
            }
        }
        any_loaded = true;
    }

    if (!any_loaded) return true;

    std::vector<float> merged_avg(expected_N, 0.0f);
    for (uint32_t i = 0; i < expected_N; ++i) {
        if (merged_count[i] > 0) {
            merged_avg[i] = static_cast<float>(merged_sum[i] / static_cast<double>(merged_count[i]));
        } else {
            merged_avg[i] = std::numeric_limits<float>::quiet_NaN();
        }
        // Scrub sentinels for faces with no contribution at all.
        const float empty = nan_when_empty ? std::numeric_limits<float>::quiet_NaN() : 0.0f;
        if (merged_max[i] == std::numeric_limits<float>::lowest())
            merged_max[i] = empty;
        if (merged_min[i] == std::numeric_limits<float>::max())
            merged_min[i] = empty;
    }

    fs::path stats_out = outdir / "statistics";
    fs::create_directories(stats_out);
    bool ok = true;
    if (want_max) ok &= write_float32_array(stats_out / (name + "_max.bin"), merged_max);
    if (want_avg) ok &= write_float32_array(stats_out / (name + "_avg.bin"), merged_avg);
    if (want_min) ok &= write_float32_array(stats_out / (name + "_min.bin"), merged_min);
    if (ok) {
        std::cout << "[merge] " << name << ": merged " << chunk_dirs.size()
                  << " chunks  (" << expected_N << " faces)\n";
    }
    return ok;
}

// ---- Subsurface depth-map discovery ---------------------------------------
//
// The depth field names encode the resolved depth (`temperature_z0p050m`), so
// they are a property of the RUN, not a compile-time list. Discover them from a
// chunk's own statistics directory rather than re-deriving them from a config
// the merge tool does not read.
std::vector<std::string> discover_depth_fields(const fs::path& chunk_dir) {
    std::vector<std::string> out;
    const fs::path stats = chunk_dir / "statistics";
    std::error_code ec;
    if (!fs::is_directory(stats, ec)) return out;
    const std::string prefix = "temperature_z";
    // ⛔ ALL THREE suffixes, de-duplicated -- not just `_max.bin`.
    // `output.map_depth_stats` is a free subset, so a run asked for
    // `[avg, min]` writes no `_max.bin` at all. Anchoring on one suffix made
    // discovery return nothing for that run, the whole depth product vanish
    // from the merged directory, and the tool still print success.
    for (const auto& entry : fs::directory_iterator(stats, ec)) {
        if (ec) break;
        const std::string name = entry.path().filename().string();
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        for (const char* suffix : { "_max.bin", "_min.bin", "_avg.bin" }) {
            const std::string suf(suffix);
            if (name.size() <= prefix.size() + suf.size()) continue;
            if (name.compare(name.size() - suf.size(), suf.size(), suf) != 0)
                continue;
            std::string field = name.substr(0, name.size() - suf.size());
            if (std::find(out.begin(), out.end(), field) == out.end())
                out.push_back(field);
            break;
        }
    }
    return out;
}

// ---- Scalar per-face fraction merge ---------------------------------------
//
// A FRACTION is not a max and not a per-chunk average: merging shadow_lit_frac
// by taking chunk 0's file, or by averaging the chunks' own fractions, is wrong
// whenever the chunks carry different epoch counts -- and chunked windows
// routinely do, because the last chunk is short. Without this the merged
// product would silently ship ONE chunk's fraction under a whole-window name.
//
// So it is merged from sum + count, the same way `avg` is, and written back as
// the scalar `statistics/<name>.bin` the Python side reads (STATISTICS_FIELDS
// entries with 'scalar': True).
bool merge_scalar_fraction(const std::string& name,
                           const std::vector<fs::path>& chunk_dirs,
                           const fs::path& outdir,
                           uint32_t expected_N) {
    std::vector<double>   merged_sum(expected_N, 0.0);
    std::vector<uint64_t> merged_count(expected_N, 0u);

    // Present in chunk 0 or not at all. Present in SOME chunks is a refusal,
    // not a partial merge: a fraction assembled from a subset of the window is
    // a number that is true of no window.
    const fs::path raw0 = chunk_dirs.front() / "statistics" / "_raw";
    const bool want = fs::exists(raw0 / (name + "_sum.bin")) &&
                      fs::exists(raw0 / (name + "_count.bin"));
    if (!want) return true;   // the product was not requested; nothing to merge

    for (const auto& cd : chunk_dirs) {
        const fs::path raw = cd / "statistics" / "_raw";
        const fs::path p_sum   = raw / (name + "_sum.bin");
        const fs::path p_count = raw / (name + "_count.bin");
        if (!fs::exists(p_sum) || !fs::exists(p_count)) {
            std::cerr << "[merge] ERROR: " << name << " is present in chunk 0 but its "
                         "statistics/_raw sum/count is missing in " << cd
                      << ". A fraction merged from a subset of the chunks is true of "
                         "no window; refusing.\n";
            return false;
        }
        uint32_t N = 0;
        std::vector<double>   sum_arr;
        std::vector<uint32_t> count_arr;
        if (!read_double_array(p_sum, N, sum_arr)   || N != expected_N) return false;
        if (!read_uint32_array(p_count, N, count_arr) || N != expected_N) return false;
        for (uint32_t i = 0; i < expected_N; ++i) {
            merged_sum[i]   += sum_arr[i];
            merged_count[i] += count_arr[i];
        }
    }

    std::vector<float> out(expected_N, std::numeric_limits<float>::quiet_NaN());
    for (uint32_t i = 0; i < expected_N; ++i)
        if (merged_count[i] > 0)
            out[i] = static_cast<float>(merged_sum[i] / static_cast<double>(merged_count[i]));

    fs::path stats_out = outdir / "statistics";
    fs::create_directories(stats_out);
    if (!write_float32_array(stats_out / (name + ".bin"), out)) return false;
    std::cout << "[merge] " << name << ": merged " << chunk_dirs.size()
              << " chunks  (" << expected_N << " faces, from _raw sum/count)\n";
    return true;
}

// ---- Ice-stability merge --------------------------------------------------

// Command-line choices for the ice products. Every value given here overrides
// what the chunks recorded, with a warning; that is how an existing run is
// re-evaluated at another age, inventory or threshold without a rerun.
struct IceCli {
    bool ice_only = false;
    std::string criterion;                 // "" = what the chunks recorded
    double J0 = 3.0e-14;     bool J0_given = false;
    double ell_m = 75.0e-6;  bool ell_given = false;
    double age_gyr = 0;      bool age_given = false;
    std::string inventory;                 // "" = recorded
    double wf = 0;           bool wf_given = false;
    double rho_g = 0;        bool rho_g_given = false;
    double rho_s = 0;        bool rho_s_given = false;
    double rho_d = 0;        bool rho_d_given = false;
    double H_m = 0;          bool H_given = false;
};

bool same_retention(const ice::IceRetentionParams& a, const ice::IceRetentionParams& b) {
    return a.age_gyr == b.age_gyr && a.inventory == b.inventory &&
           a.weight_fraction == b.weight_fraction && a.grain_density_kgm3 == b.grain_density_kgm3 &&
           a.ice_density_kgm3 == b.ice_density_kgm3 && a.hop_length_m == b.hop_length_m &&
           a.rho_surface_kgm3 == b.rho_surface_kgm3 && a.rho_deep_kgm3 == b.rho_deep_kgm3 &&
           a.density_scale_m == b.density_scale_m;
}

bool merge_ice_stability(const std::vector<fs::path>& chunk_dirs,
                         const fs::path& outdir,
                         uint32_t expected_N,
                         const IceCli& cli) {
    // Probe chunk 0 for n_layers + depth grid.
    const fs::path raw0 = chunk_dirs.front() / "statistics" / "_raw";
    const fs::path p_grid0    = raw0 / "depth_grid.bin";
    const fs::path p_sumE0    = raw0 / "ice_depth_sum_E.bin";
    const fs::path p_n0       = raw0 / "ice_depth_n_samples.txt";
    if (!fs::exists(p_sumE0) || !fs::exists(p_grid0) || !fs::exists(p_n0)) {
        if (cli.ice_only) {
            std::cerr << "[merge] ERROR: --ice-only, but chunk[0] has no saved ice sums ("
                      << raw0 << "/ice_depth_sum_E.bin, ice_depth_n_samples.txt, depth_grid.bin); "
                         "the run needs thermal.compute_ice_stability and output.write_raw_stats.\n";
            return false;
        }
        std::cout << "[merge] ice_depth: raw accumulator missing in chunk[0]; skipping.\n";
        return true;
    }

    // ---- What the chunks recorded (all must agree) -------------------------
    // Chunks written before the retention criterion carry only J0 and l (or no
    // record at all). A run that spans the upgrade mixes the two formats; that
    // is allowed, and each part of the record is compared among the chunks
    // that carry it.
    int n_with = 0, n_ret = 0, n_crit = 0;
    ice::IceProductRecord rec0, rec_ret, rec_crit;
    for (const auto& cd : chunk_dirs) {
        ice::IceProductRecord r;
        const fs::path pp = cd / "statistics" / "_raw" / "ice_depth_params.txt";
        if (!ice::readIceProductRecord(pp.string(), r)) continue;
        if (n_with > 0 &&
            (r.flux_threshold_kgm2s != rec0.flux_threshold_kgm2s || r.hop_length_m != rec0.hop_length_m)) {
            std::cerr << "[merge] ERROR: chunks were run with different ice criteria ("
                      << rec0.flux_threshold_kgm2s << " kg/m2/s, " << rec0.hop_length_m << " m vs "
                      << r.flux_threshold_kgm2s << ", " << r.hop_length_m
                      << " in " << cd << "); refusing to merge them into one ice product.\n";
            return false;
        }
        if (n_with == 0) rec0 = r;
        ++n_with;
        if (r.criterion_recorded) {
            if (n_crit > 0 && (r.criterion != rec_crit.criterion ||
                               r.flux_requested != rec_crit.flux_requested)) {
                std::cerr << "[merge] ERROR: chunks were run with different ice criteria ("
                          << rec_crit.criterion << " vs " << r.criterion << " in " << cd
                          << "); refusing to merge them into one ice product.\n";
                return false;
            }
            if (n_crit == 0) rec_crit = r;
            ++n_crit;
        }
        if (r.has_retention) {
            if (n_ret > 0 && (!same_retention(r.retention, rec_ret.retention) ||
                              r.per_face_H != rec_ret.per_face_H)) {
                std::cerr << "[merge] ERROR: chunks were run with different ice criteria (the "
                             "retention settings recorded in " << cd << " differ from those of "
                             "the other chunks); refusing to merge them into one ice product.\n";
                return false;
            }
            if (n_ret == 0) rec_ret = r;
            ++n_ret;
        }
    }
    const int n_chunks = static_cast<int>(chunk_dirs.size());
    if (n_crit > 0 && n_crit < n_chunks)
        std::cerr << "[merge] WARNING: " << n_chunks - n_crit << " of " << n_chunks
                  << " chunks were written before the retention criterion; their record is "
                     "the flux threshold only.\n";

    // ---- Legacy flux criterion: J0 and l (unchanged rules) -----------------
    double J0 = cli.J0, ell_m = cli.ell_m;
    if (n_with > 0 && n_with != static_cast<int>(chunk_dirs.size())) {
        if (!(cli.J0_given && cli.ell_given)) {
            std::cerr << "[merge] ERROR: only " << n_with << " of " << chunk_dirs.size()
                      << " chunks record their ice criterion (statistics/_raw/"
                         "ice_depth_params.txt); refusing to guess for the others. Pass "
                         "both --J0 and --ell-um to state it.\n";
            return false;
        }
        std::cerr << "[merge] WARNING: only " << n_with << " of " << chunk_dirs.size()
                  << " chunks record their ice criterion; using the --J0 / --ell-um given.\n";
    }
    if (n_with > 0) {
        if (!cli.J0_given) J0 = rec0.flux_threshold_kgm2s;
        else if (cli.J0 != rec0.flux_threshold_kgm2s)
            std::cerr << "[merge] WARNING: --J0 " << cli.J0 << " OVERRIDES the chunks' "
                         "configured threshold " << rec0.flux_threshold_kgm2s << " kg/m2/s.\n";
        if (!cli.ell_given) ell_m = rec0.hop_length_m;
        else if (cli.ell_m != rec0.hop_length_m)
            std::cerr << "[merge] WARNING: --ell-um " << cli.ell_m * 1e6 << " OVERRIDES the "
                         "chunks' configured hop length " << rec0.hop_length_m * 1e6 << " um.\n";
    } else if (!cli.J0_given || !cli.ell_given) {
        std::cerr << "[merge] WARNING: the chunks do not record their ice criterion "
                     "(written before 2026-11); using J0 = " << J0 << " kg/m2/s and l = "
                  << ell_m * 1e6 << " um" << (cli.J0_given && cli.ell_given ? "" : " (defaults)")
                  << ". Pass --J0 / --ell-um if the run was configured otherwise.\n";
    }

    // ---- Which products ------------------------------------------------------
    // Chunks written before the retention criterion (no "criterion" record)
    // produced the flux product; a re-merge keeps writing it (as
    // ice_depth_flux.bin) and adds the retention product, the default.
    const bool all_new = (n_crit == n_chunks);
    std::string criterion = !cli.criterion.empty() ? cli.criterion
                          : (all_new ? rec_crit.criterion : std::string("retention"));
    if (criterion != "retention" && criterion != "flux") {
        std::cerr << "[merge] ERROR: --criterion must be retention or flux, not '" << criterion << "'.\n";
        return false;
    }
    bool want_retention = (criterion == "retention");
    const bool want_flux = (criterion == "flux") || cli.J0_given ||
                           (all_new ? rec_crit.flux_requested : true);

    // ---- Retention parameters: recorded, then overridden -------------------
    ice::IceRetentionParams rp;
    const bool have_rp = (n_ret == n_chunks);
    if (have_rp) rp = rec_ret.retention;
    const ice::IceProductRecord& rec_r = have_rp ? rec_ret : rec0;
    std::vector<std::string> overrides;
    if (cli.ell_given && (!have_rp || cli.ell_m != rec_ret.retention.hop_length_m))
        overrides.push_back("ell-um " + std::to_string(cli.ell_m * 1e6));
    rp.hop_length_m = ell_m;
    auto over = [&](bool given, double v, double& field, const char* name, double rec_v) {
        if (!given) return;
        if (have_rp && v != rec_v) {
            std::cerr << "[merge] WARNING: --" << name << " " << v << " OVERRIDES the chunks' configured "
                      << rec_v << ".\n";
            overrides.push_back(std::string(name) + " " + std::to_string(v));
        } else if (!have_rp) {
            overrides.push_back(std::string(name) + " " + std::to_string(v) + " (not recorded by the run)");
        }
        field = v;
    };
    over(cli.age_given, cli.age_gyr, rp.age_gyr, "retention-age-gyr", rec_r.retention.age_gyr);
    over(cli.wf_given, cli.wf, rp.weight_fraction, "ice-weight-fraction", rec_r.retention.weight_fraction);
    over(cli.rho_g_given, cli.rho_g, rp.grain_density_kgm3, "grain-density", rec_r.retention.grain_density_kgm3);
    over(cli.rho_s_given, cli.rho_s, rp.rho_surface_kgm3, "rho-s", rec_r.retention.rho_surface_kgm3);
    over(cli.rho_d_given, cli.rho_d, rp.rho_deep_kgm3, "rho-d", rec_r.retention.rho_deep_kgm3);
    over(cli.H_given, cli.H_m, rp.density_scale_m, "H-m", rec_r.retention.density_scale_m);
    if (!cli.inventory.empty()) {
        ice::IceInventory inv;
        if (!ice::parseIceInventory(cli.inventory, inv)) {
            std::cerr << "[merge] ERROR: --ice-inventory must be pore_filling or weight_fraction.\n";
            return false;
        }
        if (have_rp && inv != rec_r.retention.inventory)
            std::cerr << "[merge] WARNING: --ice-inventory " << cli.inventory
                      << " OVERRIDES the chunks' configured " << ice::iceInventoryName(rec_r.retention.inventory) << ".\n";
        if (!have_rp || inv != rec_r.retention.inventory)
            overrides.push_back("ice-inventory " + cli.inventory);
        rp.inventory = inv;
    }
    if (want_retention && !have_rp) {
        // Written before the retention criterion: the density profile is not
        // recorded. It must be stated, never assumed. Without it, a merge the
        // user did not explicitly ask for retention still writes the product
        // these chunks were run for (the flux product), and says so.
        if (!(cli.rho_s_given && cli.rho_d_given && cli.H_given)) {
            const char* how = "Pass --rho-s, --rho-d and --H-m (kg m^-3, kg m^-3, m: the run's "
                              "regolith density profile, from the body defaults or its "
                              "generic_body section, and the '[property-map]' lines of its log)";
            if (cli.criterion == "retention" || cli.ice_only) {
                std::cerr << "[merge] ERROR: the chunks do not record the regolith density profile "
                             "the retention criterion needs (they were written before it "
                             "existed). " << how << ", or --criterion flux for the flux-threshold "
                             "product alone.\n";
                return false;
            }
            std::cerr << "[merge] WARNING: the chunks were written before the retention criterion "
                         "and do not record the regolith density profile it needs, so NO "
                         "ice_depth.bin (retention) is written; the flux-threshold product they "
                         "were run for is written as ice_depth_flux.bin. " << how
                      << " with --ice-only to add the retention product.\n";
            want_retention = false;
        }
    }
    if (want_retention && !have_rp) {
        std::cerr << "[merge] WARNING: the chunks do not record their retention settings; using "
                  << "age " << rp.age_gyr << " Gyr, " << ice::iceInventoryName(rp.inventory)
                  << ", grain density " << rp.grain_density_kgm3 << " kg m^-3, rho "
                  << rp.rho_surface_kgm3 << " -> " << rp.rho_deep_kgm3 << " kg m^-3, H "
                  << rp.density_scale_m << " m, l " << ell_m * 1e6 << " um, as given or by default."
                  << " A per-face H map, if the run used one, is not applied.\n";
    }
    if (!want_retention && (cli.age_given || cli.wf_given || cli.rho_g_given || cli.rho_s_given ||
                            cli.rho_d_given || cli.H_given || !cli.inventory.empty()))
        std::cerr << "[merge] WARNING: retention options were given but no retention product is "
                     "written; they are ignored.\n";
    if (!want_flux && cli.J0_given)
        std::cerr << "[merge] WARNING: --J0 is ignored.\n";
    std::cout << "[merge] ice criterion: " << (want_retention ? "retention" : "flux")
              << (want_retention && want_flux ? " (plus the flux-threshold product)" : "") << "\n";
    if (want_retention) {
        const std::string e = ice::validateIceRetentionParams(rp);
        if (!e.empty()) { std::cerr << "[merge] ERROR: retention criterion: " << e << "\n"; return false; }
    }

    // Per-face H, when the run had a map and no scalar override replaces it.
    std::vector<double> face_H;
    const bool use_face_H = want_retention && have_rp && rec_ret.per_face_H && !cli.H_given;
    if (want_retention && have_rp && rec_ret.per_face_H && cli.H_given)
        std::cerr << "[merge] WARNING: --H-m replaces the run's per-face H map on every face.\n";
    if (use_face_H) {
        for (size_t c = 0; c < chunk_dirs.size(); ++c) {
            std::vector<double> h;
            const fs::path ph = chunk_dirs[c] / "statistics" / "_raw" / "ice_column_H.bin";
            if (!ice::readFaceH(ph.string(), h) || h.size() != expected_N) {
                std::cerr << "[merge] ERROR: the chunks record a per-face H map but " << ph
                          << " is missing or has the wrong length.\n";
                return false;
            }
            if (c == 0) face_H.swap(h);
            else if (h != face_H) {
                std::cerr << "[merge] ERROR: per-face H map in " << chunk_dirs[c]
                          << " differs from chunk[0]'s.\n";
                return false;
            }
        }
    }

    // ---- Sums ---------------------------------------------------------------
    uint32_t n_layers = 0;
    std::vector<double> z_grid;
    {
        uint32_t N;
        if (!read_double_array(p_grid0, N, z_grid)) return false;
        n_layers = N;
    }

    // The flux product keeps its float32 accumulation (byte-for-byte what it
    // was); the retention product sums the chunks in double.
    ice::IceStabilityAccumulator acc;
    if (want_flux) {
        acc.allocate(static_cast<int>(expected_N), static_cast<int>(n_layers));
        for (int i = 0; i < (int)n_layers; ++i) acc.z[i] = z_grid[i];
    }
    std::vector<double> sum_d;
    if (want_retention) sum_d.assign(static_cast<size_t>(expected_N) * n_layers, 0.0);
    long long n_samples_total = 0;

    for (const auto& cd : chunk_dirs) {
        const fs::path raw = cd / "statistics" / "_raw";
        const fs::path p_sumE = raw / "ice_depth_sum_E.bin";
        const fs::path p_n    = raw / "ice_depth_n_samples.txt";
        const fs::path p_grid = raw / "depth_grid.bin";
        if (!fs::exists(p_sumE) || !fs::exists(p_n) || !fs::exists(p_grid)) {
            std::cerr << "[merge] ERROR: missing ice_depth raw files in " << cd << "\n";
            return false;
        }

        // Verify depth grid matches.
        std::vector<double> z_chunk;
        uint32_t Nl;
        if (!read_double_array(p_grid, Nl, z_chunk)) return false;
        if (Nl != n_layers || z_chunk != z_grid) {
            std::cerr << "[merge] ERROR: depth-grid mismatch between chunks (chunk "
                      << cd << ")\n";
            return false;
        }

        // Read sum_E for this chunk.
        std::ifstream is(p_sumE, std::ios::binary);
        if (!is) { std::cerr << "ERROR opening " << p_sumE << "\n"; return false; }
        uint32_t nl_h, nf_h;
        if (!read_uint32(is, nl_h) || !read_uint32(is, nf_h)) return false;
        if (nl_h != n_layers || nf_h != expected_N) {
            std::cerr << "[merge] ERROR: ice_depth sum_E shape mismatch in " << cd
                      << " (got " << nl_h << "x" << nf_h << ", expected "
                      << n_layers << "x" << expected_N << ")\n";
            return false;
        }
        const size_t total = static_cast<size_t>(nl_h) * nf_h;
        std::vector<float> sumE_chunk(total);
        is.read(reinterpret_cast<char*>(sumE_chunk.data()), total * sizeof(float));
        if (!is) { std::cerr << "ERROR reading sum_E from " << p_sumE << "\n"; return false; }

        if (want_flux)
            for (size_t i = 0; i < total; ++i) acc.sum_E[i] += sumE_chunk[i];
        if (want_retention)
            for (size_t i = 0; i < total; ++i) sum_d[i] += static_cast<double>(sumE_chunk[i]);

        // Read n_samples.
        long long n;
        std::ifstream is2(p_n);
        if (!is2 || !(is2 >> n)) { std::cerr << "ERROR reading " << p_n << "\n"; return false; }
        n_samples_total += n;
    }

    if (want_flux) acc.n_samples = static_cast<int>(n_samples_total);
    std::cout << "[merge] ice_depth: total samples = " << n_samples_total
              << " across " << chunk_dirs.size() << " chunks\n";

    fs::path stats_out = outdir / "statistics";
    fs::create_directories(stats_out);

    if (want_retention) {
        std::cout << "[merge] solving the retention criterion (ice table after "
                  << rp.age_gyr << " Gyr, " << ice::iceInventoryName(rp.inventory) << ")...\n";
        const ice::IceRetentionResult rr = ice::computeIceRetention(
            sum_d.data(), static_cast<int>(expected_N), static_cast<int>(n_layers),
            n_samples_total, z_grid, rp, use_face_H ? &face_H : nullptr, true);
        if (!rr.error.empty()) {
            std::cerr << "[merge] ERROR: retention criterion: " << rr.error << "\n";
            return false;
        }
        std::string prov = ice::iceRetentionProvenance(rp, z_grid, n_samples_total, rr,
                                                       use_face_H ? &face_H : nullptr);
        prov += "source: ptm_merge_chunks, " + std::to_string(chunk_dirs.size()) +
                " chunk(s) summed in double precision\n";
        if (!overrides.empty()) {
            prov += "overrides:";
            for (const auto& o : overrides) prov += " --" + o + ";";
            prov += "\n";
        }
        std::string werr;
        if (!ice::writeIceRetentionProducts(stats_out.string(), rr, prov, &werr)) {
            std::cerr << "[merge] ERROR: " << werr << "\n";
            return false;
        }
        std::cout << "[merge] Wrote " << (stats_out / "ice_depth.bin") << " and companions: "
                  << rr.n_top << " top layer, " << rr.n_depth << " at depth, " << rr.n_none
                  << " none within the column, " << rr.n_underflow << " surface-rate underflow, "
                  << rr.n_not_evaluated << " not evaluated\n";
        if (rr.n_not_evaluated)
            std::cerr << "[merge] WARNING: " << rr.n_not_evaluated << " faces have a non-finite or "
                         "negative mean rate and are NOT evaluated (class 255).\n";
    }

    if (want_flux) {
        std::cout << "[merge] flux criterion: " << ice::iceThresholdSummary(J0)
                  << "; l = " << ell_m * 1e6 << " um\n";
        std::vector<float> ice_depth = ice::computeIceStabilityDepths(
            acc, J0, ell_m, /*z_max_m*/ -1.0, /*use_openmp*/ true);
        if (ice_depth.size() != expected_N) {
            std::cerr << "[merge] ERROR: no ice-stability depths could be solved (J0 = " << J0
                      << ", l = " << ell_m << " m; both must be positive and finite).\n";
            return false;
        }
        if (!write_float32_array(stats_out / "ice_depth_flux.bin", ice_depth)) {
            std::cerr << "[merge] ERROR writing ice_depth_flux.bin\n";
            return false;
        }
        std::cout << "[merge] Wrote " << (stats_out / "ice_depth_flux.bin") << "\n";
        std::ofstream os(stats_out / "ice_depth_flux_threshold.txt");
        os << ice::iceThresholdSummary(J0) << "; grain_hop_length_um=" << ell_m * 1e6
           << "; surface rule mu<E(0)> <= J0\n";
        if (!os.good()) {
            std::cerr << "[merge] ERROR writing " << (stats_out / "ice_depth_flux_threshold.txt") << "\n";
            return false;
        }
    }
    return true;
}

void usage(const char* argv0) {
    std::cerr <<
"Usage: " << argv0 << " --chunks DIR1 DIR2 ... --output OUTDIR [options]\n"
"\n"
"Combine per-chunk statistics from a chunked PTM3D run into a single set\n"
"of merged statistics. Reads each chunk's statistics/_raw/ to combine\n"
"max/min/sum/count and the time-mean sublimation sums of the ice products.\n"
"\n"
"Ice products (the settings the chunks record are used unless overridden;\n"
"every override is reported):\n"
"  --ice-only                 write only the ice products (re-evaluate a\n"
"                             finished run; one run directory is one chunk,\n"
"                             with --allow-nonstandard-names)\n"
"  --criterion C              retention (statistics/ice_depth.bin) or flux\n"
"                             (statistics/ice_depth_flux.bin only)\n"
"  --retention-age-gyr A      age of the retention criterion\n"
"  --ice-inventory I          pore_filling or weight_fraction\n"
"  --ice-weight-fraction W    ice mass fraction (weight_fraction)\n"
"  --grain-density RHO        regolith grain density, kg m^-3\n"
"  --rho-s RHO --rho-d RHO --H-m H\n"
"                             regolith density profile (needed only for runs\n"
"                             written before the retention criterion)\n"
"  --ell-um L                 vapor transport length, um (both criteria)\n"
"  --J0 J                     flux threshold, kg m^-2 s^-1; also writes the\n"
"                             flux product (ice_depth_flux.bin)\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<fs::path> chunk_dirs;
    fs::path outdir;
    IceCli ice_cli;

    bool reading_chunks = false;
    bool allow_nonstandard_names = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--chunks") { reading_chunks = true; continue; }
        if (a == "--allow-nonstandard-names") {
            allow_nonstandard_names = true; reading_chunks = false; continue;
        }
        if (a == "--output" && i + 1 < argc) {
            outdir = argv[++i]; reading_chunks = false; continue;
        }
        auto num_arg = [&](const char* flag, double& v, bool& given) {
            if (a != flag || i + 1 >= argc) return false;
            v = std::stod(argv[++i]); given = true; reading_chunks = false; return true;
        };
        if (num_arg("--J0", ice_cli.J0, ice_cli.J0_given)) continue;
        if (num_arg("--ell-um", ice_cli.ell_m, ice_cli.ell_given)) { ice_cli.ell_m *= 1e-6; continue; }
        if (num_arg("--retention-age-gyr", ice_cli.age_gyr, ice_cli.age_given)) continue;
        if (num_arg("--ice-weight-fraction", ice_cli.wf, ice_cli.wf_given)) continue;
        if (num_arg("--grain-density", ice_cli.rho_g, ice_cli.rho_g_given)) continue;
        if (num_arg("--rho-s", ice_cli.rho_s, ice_cli.rho_s_given)) continue;
        if (num_arg("--rho-d", ice_cli.rho_d, ice_cli.rho_d_given)) continue;
        if (num_arg("--H-m", ice_cli.H_m, ice_cli.H_given)) continue;
        if (a == "--ice-only") { ice_cli.ice_only = true; reading_chunks = false; continue; }
        if (a == "--criterion" && i + 1 < argc) {
            ice_cli.criterion = argv[++i]; reading_chunks = false; continue;
        }
        if (a == "--ice-inventory" && i + 1 < argc) {
            ice_cli.inventory = argv[++i]; reading_chunks = false; continue;
        }
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        if (reading_chunks) { chunk_dirs.emplace_back(a); continue; }
        std::cerr << "Unknown argument: " << a << "\n"; usage(argv[0]); return 2;
    }

    if (chunk_dirs.empty() || outdir.empty()) {
        usage(argv[0]); return 2;
    }
    for (const auto& cd : chunk_dirs) {
        if (!fs::is_directory(cd)) {
            std::cerr << "Not a directory: " << cd << "\n"; return 2;
        }
    }

    // Defensive guard (2026-07-21): refuse chunk dirs whose name isn't the
    // canonical chunk_<NNN> the writers emit. This tool takes an EXPLICIT
    // chunk list (no glob), but a shell glob `chunk_*` in a re-merge command
    // would sweep in QUARANTINED dirs such as chunk_001_NOVF_GARBAGE and
    // silently corrupt the time-mean statistics (a max survives a cold garbage
    // chunk, a mean does not — the 2026-07-21 near-miss). Fail loud instead.
    // --allow-nonstandard-names overrides for a deliberate non-canonical merge.
    if (!allow_nonstandard_names) {
        auto is_canonical = [](const std::string& name) {
            if (name.rfind("chunk_", 0) != 0 || name.size() <= 6) return false;
            for (size_t k = 6; k < name.size(); ++k)
                if (!std::isdigit(static_cast<unsigned char>(name[k]))) return false;
            return true;
        };
        bool bad = false;
        for (const auto& cd : chunk_dirs) {
            const std::string base = cd.filename().string();
            if (!is_canonical(base)) {
                std::cerr << "ERROR: chunk dir '" << cd.string()
                          << "' has a non-canonical name (expected chunk_<NNN>). "
                             "This is the shell-glob-swept-a-quarantine-dir footgun "
                             "(e.g. chunk_*_GARBAGE); refusing to merge it into the "
                             "statistics. Pass an explicit canonical chunk list, or "
                             "--allow-nonstandard-names to override.\n";
                bad = true;
            }
        }
        if (bad) return 2;
    }

    fs::create_directories(outdir / "statistics");

    // --ice-only: the ice products alone, face count from the saved sums.
    if (ice_cli.ice_only) {
        uint32_t nf = 0;
        {
            std::ifstream is(chunk_dirs.front() / "statistics" / "_raw" / "ice_depth_sum_E.bin",
                             std::ios::binary);
            uint32_t nl = 0;
            if (!is || !read_uint32(is, nl) || !read_uint32(is, nf) || nf == 0) {
                std::cerr << "ERROR: --ice-only: cannot read the face count from "
                          << chunk_dirs.front() / "statistics" / "_raw" / "ice_depth_sum_E.bin" << "\n";
                return 1;
            }
        }
        std::cout << "[merge] --ice-only: " << nf << " faces, " << chunk_dirs.size() << " chunk(s)\n";
        return merge_ice_stability(chunk_dirs, outdir, nf, ice_cli) ? 0 : 1;
    }

    // Determine N (face count) from chunk 0.
    uint32_t expected_N = 0;
    {
        const fs::path probe = chunk_dirs.front() / "statistics" / "temperature_max.bin";
        if (!fs::exists(probe)) {
            std::cerr << "ERROR: cannot find " << probe
                      << " — is chunk 0 complete? was write_raw_stats=true set?\n";
            return 1;
        }
        std::ifstream is(probe, std::ios::binary);
        if (!is || !read_uint32(is, expected_N) || expected_N == 0) {
            std::cerr << "ERROR probing N from " << probe << "\n";
            return 1;
        }
    }
    std::cout << "[merge] face count: " << expected_N
              << "  chunks: " << chunk_dirs.size() << "\n";

    // Merge each per-face datatype that has the necessary _raw side files.
    static const std::vector<std::string> DATATYPES = {
        "temperature", "illumination", "direct_solar",
        "primary_reflected", "primary_ir"
    };
    bool ok = true;
    for (const auto& dt : DATATYPES) {
        if (!merge_one_datatype(dt, chunk_dirs, outdir, expected_N)) {
            std::cerr << "[merge] ERROR merging " << dt << "\n";
            ok = false;
        }
    }

    // Subsurface temperature maps (output.map_depths_m). The field names are
    // per-run, not a fixed list, so they are DISCOVERED from chunk 0 and then
    // required to be identical in every other chunk. A chunk that asked for a
    // different depth set is a refusal, not a partial merge: averaging one
    // window at 0.05 m with another at 0.06 m would produce a file whose name
    // is true of neither.
    {
        std::vector<std::string> depth_fields = discover_depth_fields(chunk_dirs.front());
        std::sort(depth_fields.begin(), depth_fields.end());
        for (size_t i = 1; i < chunk_dirs.size(); ++i) {
            std::vector<std::string> other = discover_depth_fields(chunk_dirs[i]);
            std::sort(other.begin(), other.end());
            if (other != depth_fields) {
                std::cerr << "[merge] ERROR: chunk " << chunk_dirs[i]
                          << " carries a DIFFERENT set of subsurface depth maps from "
                             "chunk 0 (" << other.size() << " vs "
                          << depth_fields.size() << " fields). The chunks were not "
                             "run with the same output.map_depths_m; refusing to "
                             "merge them into one product.\n";
                ok = false;
                depth_fields.clear();
                break;
            }
        }
        for (const auto& df : depth_fields) {
            if (!merge_one_datatype(df, chunk_dirs, outdir, expected_N,
                                    /*strict=*/true)) {
                std::cerr << "[merge] ERROR merging " << df << "\n";
                ok = false;
            }
        }
        if (!depth_fields.empty()) {
            // Carry the z-manifest through. It describes the depths, weights and
            // error budget, all of which are chunk-invariant; without it the
            // merged directory would hold maps whose depths are only readable
            // from their file names, at three decimals.
            const fs::path src = chunk_dirs.front() / "statistics" / "depth_map_z.json";
            if (fs::exists(src)) {
                std::error_code ec;
                fs::copy_file(src, outdir / "statistics" / "depth_map_z.json",
                              fs::copy_options::overwrite_existing, ec);
                if (ec) {
                    std::cerr << "[merge] WARNING: could not copy depth_map_z.json: "
                              << ec.message() << "\n";
                }
            } else {
                std::cerr << "[merge] WARNING: chunk 0 has subsurface depth maps but "
                             "no statistics/depth_map_z.json; the merged product will "
                             "not say which grid nodes they came from.\n";
            }
            std::cout << "[merge] subsurface depth maps: " << depth_fields.size()
                      << " field(s) merged\n";
        }
    }

    // Shadow depth. STRICT, and only when chunk 0 actually carries it: strict
    // mode is what turns "the chunks were run with write_raw_stats off" into a
    // refusal naming the missing files, instead of a merged directory that
    // silently has no shadow-depth maps at all while every chunk had them.
    if (fs::exists(chunk_dirs.front() / "statistics" / "shadow_depth_max.bin") ||
        fs::exists(chunk_dirs.front() / "statistics" / "shadow_depth_avg.bin") ||
        fs::exists(chunk_dirs.front() / "statistics" / "shadow_depth_min.bin")) {
        if (!merge_one_datatype("shadow_depth", chunk_dirs, outdir, expected_N,
                                /*strict=*/true, /*nan_when_empty=*/true)) {
            std::cerr << "[merge] ERROR merging shadow_depth\n";
            ok = false;
        }
    }

    // Shadow-depth fractions. These are per-face scalars over the whole window,
    // so they merge from sum + count and not by any per-chunk rule.
    for (const char* frac : { "shadow_lit_frac", "shadow_unresolved_frac" }) {
        if (!merge_scalar_fraction(frac, chunk_dirs, outdir, expected_N)) {
            std::cerr << "[merge] ERROR merging " << frac << "\n";
            ok = false;
        }
    }
    // Carry the shadow-depth manifest through. It describes the search ladder,
    // the sentinel encoding and the vertical convention, all chunk-invariant;
    // without it the merged directory holds maps whose sentinel values are not
    // documented anywhere. Written by the pass LAST, so its presence in chunk 0
    // is also evidence chunk 0's maps are complete.
    if (fs::exists(chunk_dirs.front() / "statistics" / "shadow_depth_max.bin")) {
        const fs::path src = chunk_dirs.front() / "statistics" / "shadow_depth.json";
        if (fs::exists(src)) {
            std::error_code ec;
            fs::copy_file(src, outdir / "statistics" / "shadow_depth.json",
                          fs::copy_options::overwrite_existing, ec);
            if (ec) {
                std::cerr << "[merge] WARNING: could not copy shadow_depth.json: "
                          << ec.message() << "\n";
            }
        } else {
            std::cerr << "[merge] WARNING: chunk 0 has shadow-depth maps but no "
                         "statistics/shadow_depth.json; the merged product will not "
                         "say what its +inf entries mean.\n";
        }
    }

    // Merge ice-stability accumulator (optional).
    if (!merge_ice_stability(chunk_dirs, outdir, expected_N, ice_cli)) {
        ok = false;
    }

    return ok ? 0 : 1;
}
