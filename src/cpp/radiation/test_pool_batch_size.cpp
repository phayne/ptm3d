// Unit test of radiation::planPoolBatchSize, the multi-GPU radiosity pool batch
// rule (radiation.gpu_batch_size). Pure function: no GPU, no files, no env.
//
// Covers: host-RAM tiers (ample / tight / below reserve), the corrected host
// bytes per face-step (PTM-833), device-count rounding (1-4 devices), the sweep
// cap and the split-in-two rule, the floor, the explicit override and its
// floor-raise WARNING (PTM-829), the env budget and the unknown-memory fallback,
// plus invariants over a grid of inputs. Exit 0 = all pass.
#include "radiation_exchange.h"

#include <cstdio>
#include <string>

using radiation::PoolBatchInputs;
using radiation::PoolBatchPlan;
using radiation::planPoolBatchSize;

static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; \
    std::printf("FAIL %s:%d: %s -- ", __FILE__, __LINE__, #cond); \
    std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static constexpr size_t GiB = size_t{1} << 30;

static PoolBatchInputs autoIn(int n_steps, int dev, size_t faces, size_t avail_gib,
                              size_t total_gib = 384) {
    PoolBatchInputs in;
    in.configured = 0; in.n_steps = n_steps; in.n_devices = dev; in.n_faces = faces;
    in.host_available_bytes = avail_gib * GiB; in.host_total_bytes = total_gib * GiB;
    return in;
}

static bool hasWarning(const PoolBatchPlan& p, const char* needle) {
    for (const auto& w : p.warnings) if (w.find(needle) != std::string::npos) return true;
    return false;
}

int main() {
    const size_t kAGT = 1797542;   // the 1.8 M-face core+ring site
    const size_t kW6h = 4990;

    // --- PTM-833: the per-step host cost is 9 doubles, not 6 ----------------
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(312, 4, kAGT, 300));
        CHECK(p.bytes_per_step == 72 * kAGT, "bytes_per_step %zu", p.bytes_per_step);
        CHECK(radiation::kPoolHostBytesPerFaceStep >= 63,
              "per-step cost must not under-count the measured 63 B footprint");
    }

    // --- RAM tier: ample. 3-day AGT slice (312 steps): whole sweep fits, so it
    //     is split in two (156 = the measured saturation point). -------------
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(312, 4, kAGT, 300));
        CHECK(p.batch == 156 && p.n_dispatch == 2, "batch %d dispatch %d", p.batch, p.n_dispatch);
        CHECK(p.budget_source == "host RAM", "%s", p.budget_source.c_str());
        CHECK(p.reserve_bytes == static_cast<size_t>(0.10 * 384.0 * GiB), "reserve %zu", p.reserve_bytes);
        CHECK(p.warnings.empty(), "unexpected warning");
    }
    // 30-day hourly window (~975 steps incl. eq span), same host: budget
    // 0.25 x (300 - 38.4) GiB = 65.4 GiB affords 542 steps -> 2 even batches.
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(975, 4, kAGT, 300));
        CHECK(p.affordable == 542, "affordable %lld", p.affordable);
        CHECK(p.batch == 488 && p.n_dispatch == 2, "batch %d dispatch %d", p.batch, p.n_dispatch);
        CHECK(p.batch * p.bytes_per_step <= p.budget_bytes, "over budget");
    }
    // --- RAM tier: tight (60 GiB free of 384): budget 5.4 GiB affords 44 steps;
    //     8 even batches of 39 -> rounded up to 40 (device multiple, <= 44). --
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(312, 4, kAGT, 60));
        CHECK(p.affordable == 44, "affordable %lld", p.affordable);
        CHECK(p.batch == 40 && p.n_dispatch == 8, "batch %d dispatch %d", p.batch, p.n_dispatch);
        CHECK(p.batch % 4 == 0, "not a device multiple");
        CHECK(p.warnings.empty(), "unexpected warning");
    }
    // --- RAM tier: below the reserve -> floor, WITH a warning. --------------
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(312, 4, kAGT, 30));
        CHECK(p.affordable == 0, "affordable %lld", p.affordable);
        CHECK(p.batch == 16, "batch %d", p.batch);
        CHECK(hasWarning(p, "below the floor"), "floor use must warn");
    }
    // Small host (MacBook class, 64 GiB total, 40 GiB free): reserve is the
    // 6.4 GiB fraction; budget 8.4 GiB affords 69 -> 5 batches of 63 -> 64.
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(312, 4, kAGT, 40, 64));
        CHECK(p.affordable == 69, "affordable %lld", p.affordable);
        CHECK(p.batch == 64 && p.n_dispatch == 5, "batch %d dispatch %d", p.batch, p.n_dispatch);
    }
    // Reserve minimum (4 GiB) binds on a 16 GiB host.
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(312, 2, kAGT, 10, 16));
        CHECK(p.reserve_bytes == 4 * GiB, "reserve %zu", p.reserve_bytes);
    }

    // --- Sweep cap / split rule on a small mesh (W6h: 1080 steps). ----------
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(1080, 4, kW6h, 300));
        CHECK(p.batch == 540 && p.n_dispatch == 2, "batch %d dispatch %d", p.batch, p.n_dispatch);
    }
    {   // under 2 x 150 steps: whole sweep, one dispatch
        PoolBatchPlan p = planPoolBatchSize(autoIn(200, 4, kW6h, 300));
        CHECK(p.batch == 200 && p.n_dispatch == 1, "batch %d dispatch %d", p.batch, p.n_dispatch);
    }
    {   // 299 steps: halves of 150 -> split; 150 rounds up to 152 (device multiple)
        PoolBatchPlan p = planPoolBatchSize(autoIn(299, 4, kW6h, 300));
        CHECK(p.batch == 152 && p.n_dispatch == 2, "batch %d dispatch %d", p.batch, p.n_dispatch);
    }
    {   // 298 steps: halves of 149 -> whole sweep
        PoolBatchPlan p = planPoolBatchSize(autoIn(298, 4, kW6h, 300));
        CHECK(p.batch == 298 && p.n_dispatch == 1, "batch %d dispatch %d", p.batch, p.n_dispatch);
    }

    // --- Device counts 1..4 on the 30-day window --------------------------
    {
        const int expect[5] = {0, 488, 488, 489, 488};
        for (int d = 1; d <= 4; ++d) {
            PoolBatchPlan p = planPoolBatchSize(autoIn(975, d, kAGT, 300));
            CHECK(p.batch == expect[d] && p.batch % d == 0 && p.n_dispatch == 2,
                  "d=%d batch %d dispatch %d", d, p.batch, p.n_dispatch);
        }
        // tight, 3 devices: affords 44 -> 23 batches, even 43; 45 would break the
        // budget, so 43 stays (uneven, rotated), NOT 42 (= a 24th batch and pause)
        PoolBatchPlan p = planPoolBatchSize(autoIn(975, 3, kAGT, 60));
        CHECK(p.affordable == 44 && p.batch == 43 && p.n_dispatch == 23,
              "batch %d aff %lld dispatch %d", p.batch, p.affordable, p.n_dispatch);
    }

    // --- Explicit overrides (PTM-829) ---------------------------------------
    {
        PoolBatchInputs in = autoIn(312, 4, kAGT, 300); in.configured = 96;
        PoolBatchPlan p = planPoolBatchSize(in);
        CHECK(p.batch == 96 && p.warnings.empty(), "batch %d", p.batch);
        CHECK(p.why == "radiation.gpu_batch_size", "%s", p.why.c_str());
        in.configured = 5000; p = planPoolBatchSize(in);
        CHECK(p.batch == 312 && p.n_dispatch == 1, "clamp to sweep: %d", p.batch);
        in.configured = 49; in.host_available_bytes = 1;   // explicit ignores RAM
        p = planPoolBatchSize(in);
        CHECK(p.batch == 49 && p.warnings.empty(), "explicit must not consult RAM: %d", p.batch);
        in.configured = 7; p = planPoolBatchSize(in);
        CHECK(p.batch == 16, "floor raise: %d", p.batch);
        CHECK(hasWarning(p, "below the pool floor"), "PTM-829: the raise must warn");
        CHECK(p.why.find("RAISED") != std::string::npos, "announce must not credit the config: %s", p.why.c_str());
        in.n_steps = 10; p = planPoolBatchSize(in);
        CHECK(p.batch == 10 && hasWarning(p, "below the pool floor"), "floor = sweep: %d", p.batch);
    }

    // --- Shape: whole steps per device only when it keeps count + budget -----
    {   // main's own case: the budget affords 49 of 312 on 4 devices -> 7 batches,
        // evened 45, rounded up to 48 (whole steps per device, same count, <= 49)
        PoolBatchInputs in = autoIn(312, 4, kAGT, 300); in.env_budget_mb = 6048;
        PoolBatchPlan p = planPoolBatchSize(in);
        CHECK(p.affordable == 49 && p.batch == 48 && p.n_dispatch == 7,
              "aff %lld batch %d dispatch %d", p.affordable, p.batch, p.n_dispatch);
        CHECK(p.why.find("whole steps per device") != std::string::npos, "%s", p.why.c_str());
        // affords 45: rounding 45 up to 48 would break the budget -> stay at 45,
        // uneven, 7 batches (NOT 44 = 8 batches); solveBatch rotates the remainder.
        in.env_budget_mb = 5555; p = planPoolBatchSize(in);
        CHECK(p.affordable == 45 && p.batch == 45 && p.n_dispatch == 7,
              "aff %lld batch %d dispatch %d", p.affordable, p.batch, p.n_dispatch);
        CHECK(p.why.find("uneven") != std::string::npos, "%s", p.why.c_str());
    }

    // --- Budget sources ------------------------------------------------------
    {
        PoolBatchInputs in = autoIn(312, 4, kAGT, 300); in.env_budget_mb = 4096;
        PoolBatchPlan p = planPoolBatchSize(in);
        CHECK(p.budget_source == "PTM_RAD_POOL_BATCH_BUDGET_MB", "%s", p.budget_source.c_str());
        // 4096 MiB / (72 B x 1.8 M) = 33 steps (it was 49 under the 48 B count)
        CHECK(p.affordable == 33 && p.batch == 32, "aff %lld batch %d", p.affordable, p.batch);
        PoolBatchInputs u = autoIn(312, 4, kAGT, 0, 0);   // query failed
        p = planPoolBatchSize(u);
        CHECK(p.budget_source == "fixed fallback" && hasWarning(p, "fixed 4096 MiB"),
              "unknown memory must warn: %s", p.budget_source.c_str());
    }

    // --- Edge cases -----------------------------------------------------------
    {
        PoolBatchPlan p = planPoolBatchSize(autoIn(0, 4, kAGT, 300));
        CHECK(p.batch == 1, "empty sweep batch %d", p.batch);
        p = planPoolBatchSize(autoIn(3, 4, 10, 300));
        CHECK(p.batch == 3 && p.n_dispatch == 1, "tiny sweep %d", p.batch);
        p = planPoolBatchSize(autoIn(1000, 4, 0, 300));   // zero faces: no division by 0
        CHECK(p.batch >= 1 && p.batch <= 1000, "zero faces %d", p.batch);
    }

    // --- Invariants over a grid ---------------------------------------------
    {
        const int steps[] = {1, 5, 16, 17, 100, 149, 300, 301, 312, 720, 975, 1080, 8760};
        const size_t faces[] = {100, 4990, 509337, kAGT, 2400000, 12000000};
        const size_t avail[] = {0, 2, 8, 40, 60, 120, 300};
        int n = 0;
        for (int N : steps) for (size_t F : faces) for (size_t A : avail) for (int D = 1; D <= 4; ++D) {
            PoolBatchPlan p = planPoolBatchSize(autoIn(N, D, F, A));
            const int fl = std::min(N, 4 * D);
            bool ok = p.batch >= fl && p.batch <= N && p.n_dispatch == (N + p.batch - 1) / p.batch;
            const bool floored = hasWarning(p, "below the floor");
            if (!floored) ok = ok && (static_cast<size_t>(p.batch) * p.bytes_per_step <= p.budget_bytes);
            if (!floored) {   // fewest batches the budget allows (+1 only for the split)
                const long long mb = std::min<long long>(N, p.affordable);
                long long nb = (N + mb - 1) / mb;
                if (nb == 1 && (N + 1) / 2 >= radiation::kPoolOverlapMinSteps) nb = 2;
                ok = ok && (p.n_dispatch == nb || p.batch == fl);
                // a device multiple unless rounding up would break the budget
                const int up = ((p.batch + D - 1) / D) * D;
                if (p.batch != N && p.batch % D != 0)
                    ok = ok && (static_cast<long long>(up) > p.affordable
                                || (N + up - 1) / up != p.n_dispatch);
            }
            if (!ok) { CHECK(ok, "N=%d F=%zu A=%zu D=%d -> batch %d", N, F, A, D, p.batch); }
            ++n;
        }
        CHECK(n > 0, "grid ran");
        std::printf("invariants checked over %d input combinations\n", n);
    }

    std::printf("test_pool_batch_size: %d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
