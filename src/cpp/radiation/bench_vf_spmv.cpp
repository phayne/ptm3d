// bench_vf_spmv -- time the engine's own CSR view-factor products on a referee CSR.
//
// Purpose (history marcher M0, 2026-11-09): the streamed history marcher does one
// width-2 product (IR + visible right-hand sides) per history step, sequentially
// (an IR lag of 2 steps already costs 0.24 K), so its cost is the per-product wall
// of SparseVFOperator::matvec2 on the production operator. This tool measures that
// wall with the production kernel -- not a re-implementation -- and nothing else.
//
// Read-only on the CSR. No mesh, no geometry-hash check (there is no mesh here):
// this tool must never be used to APPLY a VF, only to time one.
//
// usage: bench_vf_spmv <file.csr> [reps=10] [warmup=2] [threads=t1,t2,...]
//   Threads: the list (one load, swept), else OMP_NUM_THREADS. Kernel: PTM_SPMV_CANONICAL (0 legacy, 1 canonical;
//   unset = the ISA default, which is what a production run uses).
// Output: one line per product kind with min/median/max seconds per product and the
// effective CSR stream rate (cols + vals + row pointers, the bytes that cannot be cached),
// plus an order-independent checksum of y so kernel A/Bs are visibly the same products.
// interl2 / interl3: SparseVFOperator::matvecInterleavedCanonical (the history marcher's
// entry point, canonical on every ISA) at widths 2 and 3; their checksum sums all W outputs.
#include "radiation_exchange.h"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using radiation::SparseVFOperator;

namespace {
double nowSec() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
void stats(std::vector<double> t, double& mn, double& md, double& mx) {
    std::sort(t.begin(), t.end());
    mn = t.front(); mx = t.back(); md = t[t.size() / 2];
}
double checksum(const std::vector<double>& y) {
    long double s = 0;   // order-fixed serial sum: identical products give identical checksums
    for (double v : y) s += v;
    return static_cast<double>(s);
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: bench_vf_spmv <file.csr> [reps=10] [warmup=2] [threads=t1,t2,...]\n");
        return 1;
    }
    const std::string path = argv[1];
    const int reps = argc > 2 ? std::atoi(argv[2]) : 10;
    const int warm = argc > 3 ? std::atoi(argv[3]) : 2;
    std::vector<int> tlist;
    if (argc > 4) {
        std::string a = argv[4];
        for (size_t p = 0; p <= a.size();) {
            const size_t q = std::min(a.find(',', p), a.size());
            const int t = std::atoi(a.substr(p, q - p).c_str());
            if (t < 1) { std::fprintf(stderr, "FATAL: bad thread list '%s'\n", argv[4]); return 1; }
            tlist.push_back(t); p = q + 1;
        }
    } else {
        tlist.push_back(omp_get_max_threads());
    }
    if (reps < 1 || warm < 0) { std::fprintf(stderr, "FATAL: reps >= 1, warmup >= 0\n"); return 1; }

    const double t0 = nowSec();
    std::ifstream f(path, std::ios::binary);
    int32_t n = 0; double gh = 0.0;
    f.read(reinterpret_cast<char*>(&n), 4); f.read(reinterpret_cast<char*>(&gh), 8);
    if (!f || n <= 0) { std::fprintf(stderr, "FATAL: cannot read CSR header from %s\n", path.c_str()); return 2; }
    std::vector<int64_t> ip(static_cast<size_t>(n) + 1);
    f.read(reinterpret_cast<char*>(ip.data()), static_cast<std::streamsize>(ip.size() * 8));
    int64_t nnz = 0; f.read(reinterpret_cast<char*>(&nnz), 8);
    if (!f || nnz != ip.back() || nnz <= 0) {
        std::fprintf(stderr, "FATAL: CSR nnz %lld inconsistent with row pointers (%lld)\n",
                     static_cast<long long>(nnz), static_cast<long long>(ip.back()));
        return 2;
    }
    std::vector<int32_t> cols(static_cast<size_t>(nnz));
    std::vector<float> vals(static_cast<size_t>(nnz));
    f.read(reinterpret_cast<char*>(cols.data()), static_cast<std::streamsize>(nnz) * 4);
    f.read(reinterpret_cast<char*>(vals.data()), static_cast<std::streamsize>(nnz) * 4);
    if (!f) { std::fprintf(stderr, "FATAL: CSR read truncated\n"); return 2; }
    const double tload = nowSec() - t0;
    const double gb = (static_cast<double>(nnz) * 8.0 + (static_cast<double>(n) + 1.0) * 8.0) / 1e9;
    std::printf("csr n=%d nnz=%lld (%.1f/row) stream %.3f GB, load %.1f s\n", n,
                static_cast<long long>(nnz), static_cast<double>(nnz) / n, gb, tload);

    SparseVFOperator op(static_cast<size_t>(n), std::move(ip), std::move(cols), std::move(vals));
    // Deterministic, physically scaled inputs: IR source ~ eps sigma T^4 over 40-390 K, visible 0-1300 W/m2.
    std::vector<double> x1(static_cast<size_t>(n)), x2(static_cast<size_t>(n)), y1, y2;
    for (int32_t i = 0; i < n; ++i) {
        const double u = std::fmod(0.6180339887498949 * (i + 1), 1.0);
        const double T = 40.0 + 350.0 * u;
        x1[static_cast<size_t>(i)] = 0.95 * 5.670374419e-8 * T * T * T * T;
        x2[static_cast<size_t>(i)] = 1300.0 * std::fmod(0.7548776662466927 * (i + 1), 1.0);
    }
    auto run = [&](const char* kind, auto&& fn, const std::vector<double>& y) {
        for (int k = 0; k < warm; ++k) fn();
        std::vector<double> t(static_cast<size_t>(reps));
        for (int k = 0; k < reps; ++k) { const double a = nowSec(); fn(); t[static_cast<size_t>(k)] = nowSec() - a; }
        double mn, md, mx; stats(t, mn, md, mx);
        std::printf("threads %2d %-8s s/product min %.4f med %.4f max %.4f | stream %.1f GB/s (med) | checksum %.17g\n",
                    omp_get_max_threads(), kind, mn, md, mx, gb / md, checksum(y));
        std::fflush(stdout);
    };
    for (int t : tlist) {
        omp_set_num_threads(t);
        run("matvec2", [&] { op.matvec2(x1, x2, y1, y2); }, y1);
        std::printf("threads %2d matvec2 y2 checksum %.17g\n", t, checksum(y2));
        run("matvec", [&] { op.matvec(x1, y1); }, y1);
        // The marcher's entry point, canonical on every ISA: width 2 (IR + visible) and
        // width 3 (IR + visible + second visible bounce).
        for (int W : {2, 3}) {
            std::vector<double> X(static_cast<size_t>(n) * W), Y(static_cast<size_t>(n) * W);
            for (int32_t i = 0; i < n; ++i)
                for (int j = 0; j < W; ++j) X[static_cast<size_t>(i) * W + j] = (j == 1 ? x2 : x1)[static_cast<size_t>(i)] * (j == 2 ? 0.12 : 1.0);
            run(W == 2 ? "interl2" : "interl3", [&] { op.matvecInterleavedCanonical(W, X.data(), Y.data()); }, Y);
        }
    }
    return 0;
}
