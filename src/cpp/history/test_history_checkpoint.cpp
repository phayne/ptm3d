// test_history_checkpoint -- the byte-identical pause/resume gate for the history marcher (M1).
//
// A toy marcher made of the engine's own pieces -- SparseVFOperator::matvecInterleavedCanonical
// (width 3: IR source, next first visible bounce, this step's second bounce) and the production
// CPU column (ThermalBatchSoA::advanceTimestep, per-step load / n_sub substeps / store, REC cadence:
// 1 substep on old history, 7 on the last third) -- on a synthetic 4096-face scene.
//   G1  uninterrupted run, 8 threads                                        -> reference bytes
//   G2  pause by SIGNAL flag at epoch 37 (3 threads); resume with 11 threads; pause by PAUSE FILE at
//       epoch 251; resume REFUSED while PAUSE exists; remove it, resume with 5 threads to the end
//       -> T and the carried bounce vector must be byte-identical to G1
//   G3  a checkpoint from another job (one identity field changed) is refused, never used
//   G4  a corrupt latest checkpoint falls back to .prev LOUDLY, and the result is still identical
// Restores T alone (k, cp, kappa are rebuilt from T inside advanceTimestep; kc, rho are static),
// which is the claim that keeps AGT checkpoints at ~0.33 GB instead of ~2 GB.
#include "history/history_checkpoint.h"
#include "radiation/radiation_exchange.h"
#include "thermal/thermal_model.h"

#include <omp.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

using namespace ptm::history;
static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("  FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

namespace {
constexpr int BS = thermal::THERMAL_BATCH_SIZE;
constexpr size_t NF = 4096;
constexpr uint64_t NE = 360, FINE_FROM = 240;
constexpr int NSUB_FINE = 7;
constexpr double EPS = 0.95, SIG = 5.670374419e-8, ABOND = 0.12;

std::unique_ptr<radiation::SparseVFOperator> makeVF() {
    std::mt19937_64 rng(11); std::uniform_real_distribution<double> U(0, 1);
    std::vector<int64_t> ip(NF + 1, 0); std::vector<int32_t> cols; std::vector<float> vals;
    for (size_t i = 0; i < NF; ++i) {
        std::vector<int32_t> c; while (c.size() < 64) { int32_t j = (int32_t)(rng() % NF); if ((size_t)j != i && std::find(c.begin(), c.end(), j) == c.end()) c.push_back(j); }
        std::sort(c.begin(), c.end());
        std::vector<double> v(64); double s = 0; for (auto& x : v) { x = U(rng); s += x; }
        const double rs = 0.05 + 0.45 * U(rng);
        for (int k = 0; k < 64; ++k) { cols.push_back(c[k]); vals.push_back((float)(v[k] / s * rs)); }
        ip[i + 1] = (int64_t)cols.size();
    }
    return std::make_unique<radiation::SparseVFOperator>(NF, std::move(ip), std::move(cols), std::move(vals));
}

double absorbed(uint64_t k, size_t i) {   // direct absorbed (W/m2), synthetic diurnal, 30% never lit
    const double u = std::fmod(0.6180339887498949 * (double)(i + 1), 1.0);
    if (u < 0.3) return 0.0;
    return std::max(0.0, 1100.0 * u * std::sin(2 * M_PI * ((double)k / 120.0 + u)));
}
double reflected(uint64_t k, size_t i) { return absorbed(k, i) * 0.15 / 0.85; }

struct Toy {
    const radiation::SparseVFOperator& M;
    thermal::Planet planet{"Moon"}; thermal::Config config; std::unique_ptr<thermal::Profile> ref;
    int nl = 0; size_t nb = NF / BS;
    std::vector<double> T;    // per batch: [layer * BS + f], nb x nl x BS
    std::vector<double> b1;   // first visible bounce incident for the NEXT epoch (carried)
    uint64_t k = 0;
    explicit Toy(const radiation::SparseVFOperator& m) : M(m) {
        ref = std::make_unique<thermal::Profile>(planet, 0.0, config, false); nl = ref->nlayers;
    }
    Identity identity() const {
        Identity id; id.set("scene", "toy4096-seed11"); id.set("forcing", "synthetic-v1");
        id.set("epochs", std::to_string(NE) + "/" + std::to_string(FINE_FROM)); id.set("grid_layers", std::to_string(nl));
        return id;
    }
    void fresh() {
        T.assign(nb * nl * BS, 0.0);
        for (size_t bi = 0; bi < nb; ++bi) for (int f = 0; f < BS; ++f) {
            const double u = std::fmod(0.6180339887498949 * (double)(bi * BS + f + 1), 1.0);
            for (int l = 0; l < nl; ++l) T[(bi * nl + l) * BS + f] = 60.0 + 200.0 * u * std::exp(-l / 5.0);
        }
        std::vector<double> X(NF), Y(NF);
        for (size_t i = 0; i < NF; ++i) X[i] = reflected(0, i);
        M.matvecInterleavedCanonical(1, X.data(), Y.data()); b1 = Y; k = 0;
    }
    Checkpoint save(uint64_t next) const { Checkpoint c; c.id = identity(); c.next_epoch = next; c.blobs = {Blob::of("T", T), Blob::of("b1", b1)}; return c; }
    bool restore(const Checkpoint& c) {
        const Blob* t = c.find("T"); const Blob* b = c.find("b1");
        return t && b && t->to(T) && b->to(b1) && T.size() == nb * nl * BS && b1.size() == NF && (k = c.next_epoch, true);
    }
    void step() {
        std::vector<double> X(3 * NF), Y(3 * NF);
        for (size_t i = 0; i < NF; ++i) {
            const size_t bi = i / BS, f = i % BS; const double T0 = T[(bi * nl + 0) * BS + f];
            X[3 * i] = EPS * SIG * T0 * T0 * T0 * T0; X[3 * i + 1] = reflected(k + 1, i); X[3 * i + 2] = ABOND * b1[i];
        }
        M.matvecInterleavedCanonical(3, X.data(), Y.data());
        std::vector<double> flux(NF);
        for (size_t i = 0; i < NF; ++i) flux[i] = absorbed(k, i) + EPS * Y[3 * i] + (1 - ABOND) * (b1[i] + Y[3 * i + 2]);
        for (size_t i = 0; i < NF; ++i) b1[i] = Y[3 * i + 1];
        const int nsub = (k >= FINE_FROM) ? NSUB_FINE : 1; const double dt = planet.day / 120.0 / nsub;
        #pragma omp parallel
        {
            // Reused across batches ON PURPOSE: if any derived array were read before advanceTimestep
            // rebuilds it, the result would depend on which batches a thread saw -> the gate fails.
            thread_local std::unique_ptr<thermal::ThermalBatchSoA> sb;
            thread_local const void* owner = nullptr;
            if (!sb || owner != this) {
                sb = std::make_unique<thermal::ThermalBatchSoA>(); sb->initialize(nl, BS, planet, config);
                sb->copyGridFromProfile(*ref); owner = this;
                for (int f = 0; f < BS; ++f) { sb->emissivity[f] = planet.emissivity; sb->Qb[f] = planet.Qb; }
            }
            thermal::ThermalBatchSoA& b = *sb;
            #pragma omp for schedule(dynamic, 7)
            for (size_t bi = 0; bi < nb; ++bi) {
                b.n_faces = BS;
                std::memcpy(b.T, &T[bi * nl * BS], sizeof(double) * nl * BS);
                b.setFluxContiguous(&flux[bi * BS], BS);
                for (int s = 0; s < nsub; ++s) b.advanceTimestep(dt);
                std::memcpy(&T[bi * nl * BS], b.T, sizeof(double) * nl * BS);
            }
        }
        ++k;
    }
    // Runs to the end or to a pause. Returns 0 (done) or kExitPaused; -1 on error.
    int run(const std::string& dir, int threads, uint64_t signal_at, const char* pause_file_at_note, uint64_t file_at) {
        omp_set_num_threads(threads);
        PauseControl pc(dir); CheckpointPolicy pol; pol.every = 50;
        if (!dir.empty()) {
            if (pc.pauseFilePresent()) return kExitPaused + 100;   // refuse to resume while PAUSE exists
            Checkpoint c; bool found; std::string note;
            const std::string e = loadCheckpoint(dir, identity(), c, found, note);
            if (!e.empty()) { std::printf("    load refused: %s\n", e.c_str()); return -1; }
            if (!note.empty()) std::printf("    %s\n", note.c_str());
            if (found) { if (!restore(c)) return -1; std::printf("    resumed at epoch %llu with %d threads\n", (unsigned long long)k, threads); }
            else fresh();
        } else fresh();
        (void)pause_file_at_note;
        while (k < NE) {
            step();
            if (k == signal_at) PauseControl::requestForTest();
            if (k == file_at && !dir.empty()) { FILE* f = std::fopen((dir + "/PAUSE").c_str(), "w"); if (f) std::fclose(f); }
            if (!dir.empty() && (pc.requested() || pol.due(k))) {
                const std::string e = writeCheckpoint(dir, save(k));
                if (!e.empty()) { std::printf("    checkpoint failed: %s\n", e.c_str()); return -1; }
                if (pc.requested()) { PauseControl::clearForTest(); std::printf("    paused after epoch %llu\n", (unsigned long long)k); return kExitPaused; }
            }
        }
        if (!dir.empty()) { const std::string e = writeCheckpoint(dir, save(k)); if (!e.empty()) return -1; }
        return 0;
    }
};

bool same(const Toy& a, const Toy& b) {
    return a.T.size() == b.T.size() && a.b1.size() == b.b1.size() && a.k == b.k &&
           std::memcmp(a.T.data(), b.T.data(), 8 * a.T.size()) == 0 && std::memcmp(a.b1.data(), b.b1.data(), 8 * a.b1.size()) == 0;
}
std::string mkdir_tmp(const char* tag) {
    char t[256]; std::snprintf(t, sizeof t, "/tmp/test_history_checkpoint_%s_%d_XXXXXX", tag, (int)getpid());
    const char* d = mkdtemp(t); return d ? d : "";
}
void rmdir_all(const std::string& d) {
    for (const char* f : {"/history.ckpt", "/history.ckpt.prev", "/history.ckpt.tmp", "/PAUSE"}) std::remove((d + f).c_str());
    rmdir(d.c_str());
}
}  // namespace

int main() {
    auto M = makeVF();
    std::printf("G1 uninterrupted, 8 threads\n");
    Toy A(*M); CHECK(A.run("", 8, ~0ull, nullptr, ~0ull) == 0, "reference run failed");
    double tmax = 0, tmin = 1e9; for (size_t bi = 0; bi < A.nb; ++bi) for (int f = 0; f < BS; ++f) { tmax = std::max(tmax, A.T[bi * A.nl * BS + f]); tmin = std::min(tmin, A.T[bi * A.nl * BS + f]); }
    std::printf("  %llu epochs, %d layers, surface T %.1f..%.1f K, digest T %s b1 %s\n", (unsigned long long)A.k, A.nl, tmin, tmax,
                hex(digest(A.T.data(), 8 * A.T.size())).c_str(), hex(digest(A.b1.data(), 8 * A.b1.size())).c_str());

    std::printf("G2 pause (signal @37, 3 thr) -> resume (11 thr) -> pause (PAUSE file @251) -> refused while PAUSE -> resume (5 thr)\n");
    const std::string d = mkdir_tmp("g2");
    { Toy B(*M); CHECK(B.run(d, 3, 37, nullptr, ~0ull) == kExitPaused, "signal pause did not pause"); }
    { Toy B(*M); CHECK(B.run(d, 11, ~0ull, nullptr, 251) == kExitPaused, "PAUSE file did not pause"); }
    { Toy B(*M); const int r = B.run(d, 5, ~0ull, nullptr, ~0ull); CHECK(r == kExitPaused + 100, "resume ran with PAUSE present (rc %d)", r);
      std::printf("    resume with PAUSE present -> %s\n", r == kExitPaused + 100 ? "refused" : "RAN"); }
    std::remove((d + "/PAUSE").c_str());
    Toy B(*M); CHECK(B.run(d, 5, ~0ull, nullptr, ~0ull) == 0, "final resume failed");
    std::printf("  final digest T %s b1 %s -> %s\n", hex(digest(B.T.data(), 8 * B.T.size())).c_str(),
                hex(digest(B.b1.data(), 8 * B.b1.size())).c_str(), same(A, B) ? "BYTE-IDENTICAL to G1" : "DIFFERENT");
    CHECK(same(A, B), "resumed history differs from the uninterrupted one");

    std::printf("G3 checkpoint from another job\n");
    {
        Toy C(*M); Identity other = C.identity(); other.set("forcing", "synthetic-v2");
        Checkpoint c; bool found; std::string note;
        const std::string e = loadCheckpoint(d, other, c, found, note);
        std::printf("    %s\n", e.empty() ? "ACCEPTED" : ("refused: " + e.substr(0, e.find('\n')) + " / " + e.substr(e.find('\n') + 1, 60)).c_str());
        CHECK(!e.empty() && !found, "a checkpoint with a different identity was accepted");
    }
    rmdir_all(d);

    std::printf("G4 corrupt latest checkpoint -> previous one, loudly\n");
    const std::string d4 = mkdir_tmp("g4");
    { Toy D(*M); CHECK(D.run(d4, 4, 180, nullptr, ~0ull) == kExitPaused, "G4 pause failed"); }   // ckpts at 150 (policy) and 180 (pause)
    {   // flip a byte in the middle of history.ckpt
        FILE* f = std::fopen((d4 + "/history.ckpt").c_str(), "r+b"); CHECK(f != nullptr, "no checkpoint");
        if (f) { std::fseek(f, 0, SEEK_END); const long n = std::ftell(f); std::fseek(f, n / 2, SEEK_SET); int c = std::fgetc(f);
                 std::fseek(f, n / 2, SEEK_SET); std::fputc(c ^ 0x10, f); std::fclose(f); }
    }
    Toy D(*M); CHECK(D.run(d4, 6, ~0ull, nullptr, ~0ull) == 0, "G4 resume failed");
    std::printf("  -> %s\n", same(A, D) ? "BYTE-IDENTICAL to G1" : "DIFFERENT");
    CHECK(same(A, D), "fallback resume differs");
    rmdir_all(d4);

    std::printf("%s (%d failure%s)\n", g_fail ? "FAIL" : "PASS", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
