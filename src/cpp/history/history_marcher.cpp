// history_marcher.cpp -- see history_marcher.h.
#include "history_marcher.h"

#include <omp.h>

#include <algorithm>
#include "thermal/eq_certificate.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace ptm {
namespace history {

namespace {
constexpr int BS = thermal::THERMAL_BATCH_SIZE;
constexpr double SIG = 5.670374419e-8;
double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

bool preadAll(int fd, void* buf, size_t n, uint64_t off) {
    char* p = static_cast<char*>(buf);
    while (n) { const ssize_t r = ::pread(fd, p, n, (off_t)off); if (r <= 0) return false; p += r; n -= (size_t)r; off += (uint64_t)r; }
    return true;
}
uint64_t fileDigest(int fd, uint64_t size) {   // word hash of the whole file, streamed
    std::vector<unsigned char> buf(1 << 24); uint64_t h = 1469598103934665603ULL, off = 0;
    std::vector<unsigned char> tail;
    while (off < size) {
        const size_t n = (size_t)std::min<uint64_t>(buf.size(), size - off);
        if (!preadAll(fd, buf.data(), n, off)) return 0;
        size_t i = 0;
        for (; i + 8 <= n; i += 8) { uint64_t w; std::memcpy(&w, buf.data() + i, 8); h ^= w; h *= 1099511628211ULL; }
        if (i < n) { uint64_t w = 0; std::memcpy(&w, buf.data() + i, n - i); h ^= w; h *= 1099511628211ULL; }   // only the last chunk can have a tail (chunks are 8-aligned)
        off += n;
    }
    return h;
}
}  // namespace

// ---- forcing -------------------------------------------------------------------------------
ForcingFile::~ForcingFile() { if (fd_ >= 0) ::close(fd_); }

std::string ForcingFile::open(const std::string& path) {
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) return "cannot open forcing file " + path;
    struct stat st; if (fstat(fd_, &st) != 0) return "cannot stat " + path;
    char mg[8];
    if (!preadAll(fd_, mg, 8, 0) || std::memcmp(mg, "PTMHF001", 8) != 0) return path + ": not a PTMHF001 forcing file";
    if (!preadAll(fd_, &ne_, 8, 8) || !preadAll(fd_, &nf_, 8, 16) || ne_ < 3 || nf_ == 0) return path + ": bad header";
    et_.resize(ne_);
    if (!preadAll(fd_, et_.data(), 8 * ne_, 24)) return path + ": truncated epoch list";
    data_off_ = 24 + 8 * ne_;
    const uint64_t want = data_off_ + ne_ * nf_ * 16;
    if ((uint64_t)st.st_size != want)
        return path + ": size " + std::to_string(st.st_size) + " != " + std::to_string(want) + " expected from its header";
    for (uint64_t k = 1; k < ne_; ++k) if (!(et_[k] > et_[k - 1])) return path + ": epochs not strictly increasing";
    digest_ = fileDigest(fd_, (uint64_t)st.st_size);
    return "";
}

bool ForcingFile::read(uint64_t k, double* a, double* r) const {
    if (k >= ne_) return false;
    const uint64_t off = data_off_ + k * nf_ * 16;
    return preadAll(fd_, a, 8 * nf_, off) && preadAll(fd_, r, 8 * nf_, off + 8 * nf_);
}

// ---- CSR -------------------------------------------------------------------------------------
std::unique_ptr<radiation::SparseVFOperator> loadCsr(const std::string& path, uint64_t n_expected, uint64_t& digest,
                                                     std::string& err) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { err = "cannot open CSR " + path; return nullptr; }
    struct stat st; fstat(fd, &st);
    int32_t n = 0; double gh = 0; int64_t nnz = 0;
    if (!preadAll(fd, &n, 4, 0) || !preadAll(fd, &gh, 8, 4) || n <= 0 || (uint64_t)n != n_expected) {
        ::close(fd); err = path + ": face count " + std::to_string(n) + " != forcing faces " + std::to_string(n_expected); return nullptr;
    }
    std::vector<int64_t> ip((size_t)n + 1);
    if (!preadAll(fd, ip.data(), 8 * ip.size(), 12) || !preadAll(fd, &nnz, 8, 12 + 8 * ip.size()) || nnz != ip.back() || nnz < 0) {
        ::close(fd); err = path + ": bad row pointers"; return nullptr;
    }
    const uint64_t co = 12 + 8 * ip.size() + 8, vo = co + 4 * (uint64_t)nnz;
    if ((uint64_t)st.st_size != vo + 4 * (uint64_t)nnz) { ::close(fd); err = path + ": size does not match nnz"; return nullptr; }
    std::vector<int32_t> cols((size_t)nnz); std::vector<float> vals((size_t)nnz);
    if (!preadAll(fd, cols.data(), 4 * (size_t)nnz, co) || !preadAll(fd, vals.data(), 4 * (size_t)nnz, vo)) { ::close(fd); err = path + ": truncated"; return nullptr; }
    digest = fileDigest(fd, (uint64_t)st.st_size);
    ::close(fd);
    size_t nover = 0; double neg = 0;
    for (int32_t i = 0; i < n; ++i) {
        double s = 0;
        for (int64_t k = ip[i]; k < ip[i + 1]; ++k) {
            s += vals[k]; neg = std::min(neg, (double)vals[k]);
            if (cols[k] < 0 || cols[k] >= n) { err = path + ": column index out of range"; return nullptr; }
        }
        if (s > 1.0 + 1e-6) ++nover;
    }
    if (nover || neg < -1e-9) { err = path + ": violates VF invariants (rows > 1: " + std::to_string(nover) + ", min entry " + std::to_string(neg) + ")"; return nullptr; }
    return std::make_unique<radiation::SparseVFOperator>((size_t)n, std::move(ip), std::move(cols), std::move(vals));
}

// ---- tau1 sizing ----------------------------------------------------------------------------
Tau1Report estimateTau1(const Column& col, const ForcingFile& F, uint64_t k0, uint64_t k1, double qb, double eps,
                        double period_s) {
    Tau1Report rep;
    const size_t nf = F.numFaces(); const int nl = col.nl;
    std::vector<double> mean(nf, 0.0), a(nf), r(nf);
    uint64_t ns = 0;
    for (uint64_t k = k0; k < k1; k += 12, ++ns) {
        if (!F.read(k, a.data(), r.data())) { rep.error = "forcing read failed"; return rep; }
        for (size_t i = 0; i < nf; ++i) mean[i] += a[i];
    }
    if (!ns) { rep.error = "empty span"; return rep; }
    std::vector<double> Rk(nl, 0.0);
    for (int l = 1; l < nl; ++l) Rk[l] = Rk[l - 1] + col.dz[l - 1] / col.kc[l - 1];
    std::vector<double> tau(nf, -1.0);
    const size_t nb = (nf + BS - 1) / BS;
    #pragma omp parallel
    {
        thermal::ThermalBatchSoA* b = new thermal::ThermalBatchSoA();
        b->initialize(nl, BS, col.planet, col.config); b->copyGridFromProfile(*col.ref);
        std::vector<double> T(nl), k(nl), rho(nl), cp(nl);
        #pragma omp for schedule(dynamic, 64)
        for (size_t bi = 0; bi < nb; ++bi) {
            for (int f = 0; f < BS; ++f) {
                const size_t i = std::min(bi * BS + f, nf - 1);
                const double T0 = std::pow((mean[i] / ns + qb) / (eps * SIG), 0.25);
                for (int l = 0; l < nl; ++l) b->T[l * BS + f] = T0 + qb * Rk[l];
            }
            thermal::updateMaterialPropertiesBatchSIMD(*b);   // the production k(T), cp(T), rho
            for (int f = 0; f < BS && bi * BS + f < nf; ++f) {
                for (int l = 0; l < nl; ++l) { T[l] = b->T[l * BS + f]; k[l] = b->k[l * BS + f]; rho[l] = b->rho[l * BS + f]; cp[l] = b->cp[l * BS + f]; }
                tau[bi * BS + f] = thermal::eqColumnTau1(col.z.data(), T.data(), k.data(), rho.data(), cp.data(), nl, eps, SIG);
            }
        }
        delete b;
    }
    std::vector<double> lit, dark;
    for (size_t i = 0; i < nf; ++i) { if (!(tau[i] > 0)) continue; (mean[i] > 0 ? lit : dark).push_back(tau[i] / period_s); }
    auto pct = [](std::vector<double> v, double q) { if (v.empty()) return 0.0; std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, (size_t)(q * (v.size() - 1) + 0.5))]; };
    rep.lit_p50 = pct(lit, 0.5); rep.lit_p99 = pct(lit, 0.99); rep.dark_p50 = pct(dark, 0.5); rep.dark_p99 = pct(dark, 0.99);
    rep.n_lit = lit.size(); rep.n_dark = dark.size();
    rep.tau1_lun = std::max(rep.lit_p99, rep.dark_p99);
    return rep;
}

// ---- marcher -------------------------------------------------------------------------------
Marcher::Marcher(const MarcherConfig& c, const Column& col, const radiation::SparseVFOperator& M, const ForcingFile& F,
                 const Identity& id)
    : c_(c), col_(col), M_(M), F_(F), id_(id) {
    nf_ = F_.numFaces(); nb_ = (nf_ + BS - 1) / BS;
}

double Marcher::tgeo() const { return std::pow(c_.qb / (c_.emissivity * SIG), 0.25); }

void Marcher::visBounces(const std::vector<double>& src, std::vector<double>& b1, std::vector<double>& vis) {
    std::vector<double> X(nf_), b2(nf_);
    b1.assign(nf_, 0.0); vis.assign(nf_, 0.0);
    M_.matvecInterleavedCanonical(1, src.data(), b1.data());
    for (size_t i = 0; i < nf_; ++i) X[i] = c_.a_bond * b1[i];
    M_.matvecInterleavedCanonical(1, X.data(), b2.data());
    for (size_t i = 0; i < nf_; ++i) vis[i] = (1.0 - c_.a_bond) * (b1[i] + b2[i]);
}

void Marcher::columnStep(const std::vector<double>& fa, const std::vector<double>& fb, double dt, double target) {
    int ns = (int)std::ceil(dt / target); if (ns < 1) ns = 1;
    const double h = dt / ns; const int nl = col_.nl;
    #pragma omp parallel
    {
        // One scratch batch per thread, reused across batches: advanceTimestep rebuilds k, cp, kappa
        // from T before use and kc, rho are static, so only T moves in and out (gated byte-identical).
        thread_local std::unique_ptr<thermal::ThermalBatchSoA> sb;
        thread_local const Column* owner = nullptr;
        thread_local double owner_qb = -1, owner_eps = -1;
        if (!sb || owner != &col_ || owner_qb != c_.qb || owner_eps != c_.emissivity) {
            sb = std::make_unique<thermal::ThermalBatchSoA>();
            sb->initialize(nl, BS, col_.planet, col_.config); sb->copyGridFromProfile(*col_.ref);
            for (int f = 0; f < BS; ++f) { sb->emissivity[f] = c_.emissivity; sb->Qb[f] = c_.qb; }
            owner = &col_; owner_qb = c_.qb; owner_eps = c_.emissivity;
        }
        thermal::ThermalBatchSoA& b = *sb;
        alignas(64) double fl[BS];
        #pragma omp for schedule(dynamic, 64)
        for (size_t bi = 0; bi < nb_; ++bi) {
            const int nfb = (int)std::min<size_t>(BS, nf_ - bi * BS);
            b.n_faces = nfb;
            std::memcpy(b.T, &T_[bi * nl * BS], sizeof(double) * nl * BS);
            for (int s = 0; s < ns; ++s) {
                const double al = (s + 0.5) / ns;
                for (int f = 0; f < BS; ++f) {
                    const size_t i = bi * BS + f;
                    fl[f] = (f < nfb) ? (1.0 - al) * fa[i] + al * fb[i] : 0.0;
                }
                b.setFluxContiguous(fl, BS);
                b.advanceTimestep(h);
            }
            std::memcpy(&T_[bi * nl * BS], b.T, sizeof(double) * nl * BS);
        }
    }
}

void Marcher::seed() {
    const int nl = col_.nl; const uint64_t P = (uint64_t)c_.samples_per_lunation;
    const uint64_t i1 = c_.i0 + (uint64_t)c_.seed_lunations * P;
    std::vector<uint64_t> idx;
    for (uint64_t k = c_.i0; k <= i1; k += (uint64_t)c_.seed_stride) idx.push_back(k);
    const size_t ns = idx.size();
    std::printf("[history] seed: %zu samples (every %d of epochs %llu..%llu), %d IR sweep(s) from radiative-equilibrium walls, "
                "replay x%d at P/%.0f%s\n", ns - 1, c_.seed_stride, (unsigned long long)c_.i0, (unsigned long long)i1,
                c_.seed_ir_sweeps, c_.seed_reps, c_.seed_target_div, c_.seed_aitken ? " + Aitken" : "");
    // The seed block's forcing is temperature-independent and replayed 3 times; it is held in FLOAT32
    // (481 x n_faces: 3.5 GB at 1.8M faces instead of 6.9 GB in float64). Relative rounding 6e-8 of the
    // absorbed flux -- far below the 0.3 K product bar (verified against the float64 numpy REC marcher).
    std::vector<std::vector<float>> Fb(ns, std::vector<float>(nf_));
    std::vector<double> a(nf_), r(nf_), b1, vis, IR(nf_), X(nf_);
    for (size_t j = 0; j < ns; ++j) {
        if (!F_.read(idx[j], a.data(), r.data())) { std::fprintf(stderr, "[history] FATAL: forcing read failed at epoch %llu\n", (unsigned long long)idx[j]); std::exit(3); }
        visBounces(r, b1, vis);
        std::vector<double> Fnv(nf_);
        for (size_t i = 0; i < nf_; ++i) Fnv[i] = a[i] + vis[i];
        std::fill(IR.begin(), IR.end(), 0.0);
        for (int it = 0; it < c_.seed_ir_sweeps; ++it) {
            for (size_t i = 0; i < nf_; ++i) X[i] = Fnv[i] + c_.qb + IR[i];
            M_.matvecInterleavedCanonical(1, X.data(), IR.data());
            for (size_t i = 0; i < nf_; ++i) IR[i] *= c_.emissivity;
        }
        for (size_t i = 0; i < nf_; ++i) Fb[j][i] = static_cast<float>(Fnv[i] + IR[i]);
    }
    // start state: mean forcing over the block (the closing sample excluded), conductive profile
    std::vector<double> Fm(nf_, 0.0);
    for (size_t j = 0; j + 1 < ns; ++j) for (size_t i = 0; i < nf_; ++i) Fm[i] += Fb[j][i];
    for (size_t i = 0; i < nf_; ++i) Fm[i] /= (double)(ns - 1);
    std::vector<double> Rk(nl, 0.0);
    for (int l = 1; l < nl; ++l) Rk[l] = Rk[l - 1] + col_.dz[l - 1] / col_.kc[l - 1];
    T_.assign(nb_ * nl * BS, 0.0);
    for (size_t bi = 0; bi < nb_; ++bi) for (int f = 0; f < BS; ++f) {
        const size_t i = std::min(bi * BS + f, nf_ - 1);   // padded lanes copy the last face (never output)
        const double T0 = std::pow((Fm[i] + c_.qb) / (c_.emissivity * SIG), 0.25);
        for (int l = 0; l < nl; ++l) T_[(bi * nl + l) * BS + f] = T0 + c_.qb * Rk[l];
    }
    // periodic replay of the block: row ns-1 is replaced by row 0 (the block closes on itself)
    const double target = c_.period_s / c_.seed_target_div;
    std::vector<double> h0 = T_, h1, h2, fa(nf_), fb(nf_);
    for (int rep = 0; rep < c_.seed_reps; ++rep) {
        for (size_t j = 0; j + 1 < ns; ++j) {
            const std::vector<float>& ra = Fb[j]; const std::vector<float>& rb = (j + 2 == ns) ? Fb[0] : Fb[j + 1];
            for (size_t i = 0; i < nf_; ++i) { fa[i] = ra[i]; fb[i] = rb[i]; }
            columnStep(fa, fb, F_.et()[idx[j + 1]] - F_.et()[idx[j]], target);
        }
        h2.swap(h1); h1.swap(h0); h0 = T_;   // h0 newest, h1 previous, h2 the one before
        if (c_.seed_aitken && rep == c_.seed_reps - 2 && c_.seed_reps >= 2) {
            const double tg = tgeo(); size_t used = 0;
            for (size_t bi = 0; bi < nb_; ++bi) for (int f = 0; f < BS; ++f) {
                const size_t deep = (bi * nl + (nl - 1)) * BS + f;
                const double d1 = h1[deep] - h2[deep], d2 = h0[deep] - h1[deep], rho = d2 / d1;
                const bool ok = std::isfinite(rho) && rho > 0 && rho < 0.9 && std::fabs(d1) > 1e-6;
                const double ext = ok ? rho / (1 - rho) : 0.0; used += ok;
                for (int l = 0; l < nl; ++l) {
                    const size_t q = (bi * nl + l) * BS + f;
                    T_[q] = std::max(h0[q] + ext * (h0[q] - h1[q]), tg);
                }
            }
            h0 = T_;
            std::printf("[history] seed: Aitken after pass %d on %zu of %zu columns (incl. padding lanes)\n", rep + 1, used, nb_ * BS);
        }
    }
    // visible pipeline state for the march: vis at i0, first bounce at i0 + 1
    if (!F_.read(c_.i0, a.data(), r.data())) std::exit(3);
    visBounces(r, b1, vis_k_);
    if (!F_.read(c_.i0 + 1, a.data(), r.data())) std::exit(3);
    b1_next_.assign(nf_, 0.0); M_.matvecInterleavedCanonical(1, r.data(), b1_next_.data());
    k_ = c_.i0;
}

void Marcher::step(uint64_t k) {
    const int nl = col_.nl; const double A = c_.a_bond, eps = c_.emissivity;
    std::vector<double> ak(nf_), ak1(nf_), rk(nf_), r2(nf_, 0.0), tmp(nf_);
    if (!F_.read(k, ak.data(), rk.data()) || !F_.read(k + 1, ak1.data(), tmp.data())) {
        std::fprintf(stderr, "[history] FATAL: forcing read failed at epoch %llu\n", (unsigned long long)k); std::exit(3);
    }
    if (k + 2 < F_.numEpochs()) { if (!F_.read(k + 2, tmp.data(), r2.data())) std::exit(3); }
    else {   // the history ends at least 2 epochs before the forcing does (checked by the caller)
        std::fprintf(stderr, "[history] FATAL: epoch %llu needs forcing at k+2 beyond the file\n", (unsigned long long)k); std::exit(3);
    }
    std::vector<double> X(3 * nf_), Y(3 * nf_);
    for (size_t i = 0; i < nf_; ++i) {
        const size_t bi = i / BS, f = i % BS; const double T0 = T_[(bi * nl + 0) * BS + f];
        X[3 * i] = eps * SIG * T0 * T0 * T0 * T0; X[3 * i + 1] = r2[i]; X[3 * i + 2] = A * b1_next_[i];
    }
    M_.matvecInterleavedCanonical(3, X.data(), Y.data());
    std::vector<double> fa(nf_), fb(nf_), vis1(nf_);
    for (size_t i = 0; i < nf_; ++i) {
        const double IR = eps * Y[3 * i];
        vis1[i] = (1.0 - A) * (b1_next_[i] + Y[3 * i + 2]);
        fa[i] = ak[i] + vis_k_[i] + IR;
        fb[i] = ak1[i] + vis1[i] + IR;
    }
    const double target = c_.period_s / ((k < c_.fine_from) ? c_.old_target_div : c_.fine_target_div);
    columnStep(fa, fb, F_.et()[k + 1] - F_.et()[k], target);
    vis_k_.swap(vis1);
    for (size_t i = 0; i < nf_; ++i) b1_next_[i] = Y[3 * i + 1];
}

Checkpoint Marcher::save(uint64_t next, const char* phase) const {
    Checkpoint c; c.id = id_; c.next_epoch = next;
    std::vector<char> ph(phase, phase + std::strlen(phase));
    c.blobs = {Blob::of("T", T_), Blob::of("b1_next", b1_next_), Blob::of("vis_k", vis_k_), Blob::of("phase", ph)};
    return c;
}

bool Marcher::restore(const Checkpoint& c) {
    const Blob* t = c.find("T"); const Blob* b = c.find("b1_next"); const Blob* v = c.find("vis_k");
    if (!t || !b || !v || !t->to(T_) || !b->to(b1_next_) || !v->to(vis_k_)) return false;
    if (T_.size() != nb_ * col_.nl * BS || b1_next_.size() != nf_ || vis_k_.size() != nf_) return false;
    k_ = c.next_epoch;
    return k_ >= c_.i0 && k_ <= c_.ti;
}

int Marcher::run() {
    PauseControl pc(c_.job_dir);
    if (pc.pauseFilePresent()) {
        std::fprintf(stderr, "[history] %s/PAUSE exists: this job is paused. Remove the file to resume.\n", c_.job_dir.c_str());
        return kExitPaused;
    }
    Checkpoint ck; bool found = false; std::string note;
    const std::string e = loadCheckpoint(c_.job_dir, id_, ck, found, note);
    if (!e.empty()) { std::fprintf(stderr, "[history] FATAL: %s\n", e.c_str()); return 4; }
    if (!note.empty()) std::fprintf(stderr, "[history] %s\n", note.c_str());
    const double t0 = now();
    if (found) {
        if (!restore(ck)) { std::fprintf(stderr, "[history] FATAL: checkpoint content does not fit this job\n"); return 4; }
        std::printf("[history] resumed from checkpoint at epoch %llu (handoff at %llu), %d threads\n",
                    (unsigned long long)k_, (unsigned long long)c_.ti, omp_get_max_threads());
    } else {
        seed();
        const std::string w = writeCheckpoint(c_.job_dir, save(k_, "seeded"));
        if (!w.empty()) { std::fprintf(stderr, "[history] FATAL: checkpoint: %s\n", w.c_str()); return 3; }
        std::printf("[history] seed done in %.1f s; marching epochs %llu..%llu (fine target from %llu), %d threads\n", now() - t0,
                    (unsigned long long)c_.i0, (unsigned long long)c_.ti, (unsigned long long)c_.fine_from, omp_get_max_threads());
    }
    CheckpointPolicy pol; pol.every = c_.checkpoint_every;
    while (k_ < c_.ti) {
        step(k_); ++k_;
        const bool pause = pc.requested();
        if (pause || pol.due(k_ - c_.i0) || k_ == c_.ti) {
            const std::string w = writeCheckpoint(c_.job_dir, save(k_, k_ == c_.ti ? "handoff" : "march"));
            if (!w.empty()) { std::fprintf(stderr, "[history] FATAL: checkpoint: %s\n", w.c_str()); return 3; }
        }
        if ((k_ - c_.i0) % (uint64_t)c_.samples_per_lunation == 0 || k_ == c_.ti) {
            double s = 0; for (size_t i = 0; i < nf_; ++i) s += T_[((i / BS) * col_.nl) * BS + i % BS];
            std::printf("[history] epoch %llu/%llu  ET %.3f  mean surface %.3f K  %.1f s\n", (unsigned long long)k_,
                        (unsigned long long)c_.ti, F_.et()[k_], s / nf_, now() - t0);
            std::fflush(stdout);
        }
        if (pause && k_ < c_.ti) {
            std::printf("[history] PAUSED after epoch %llu (checkpoint written). Resume: remove %s/PAUSE if present and rerun the "
                        "same command.\n", (unsigned long long)k_, c_.job_dir.c_str());
            return kExitPaused;
        }
    }
    return 0;
}

std::vector<double> Marcher::stateFaceMajor() const {
    const int nl = col_.nl; std::vector<double> out(nf_ * nl);
    for (size_t i = 0; i < nf_; ++i) for (int l = 0; l < nl; ++l) out[i * nl + l] = T_[((i / BS) * nl + l) * BS + i % BS];
    return out;
}

}  // namespace history
}  // namespace ptm
