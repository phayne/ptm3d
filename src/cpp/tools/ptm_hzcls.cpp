// ptm_hzcls -- build / apply a per-site solar horizon-classifier envelope (.hzc).
// History marcher M0 (2026-11-09). Engine-tree port of the 2026-11-08 measurement tool
// hzcls.cpp (rule v2; not in the repo). Same arguments, same HZ_*
// environment, same .cls.u8 / .faces.u32 / .envrays.u64 / .nmulti.u32 outputs, so its output
// can be compared byte for byte with the tool's (test_horizon_classifier --parity).
// Not used by ptm_thermal_sim.
//
// usage: ptm_hzcls core.ply ring.ply|- track.bin out_prefix [threads]
//   track.bin: int64 N, then N x 7 float64 (et, sun xyz km, earth xyz km), IAU_MOON, LT+S, observer MOON.
// env (as the tool): HZ_DC (0.125) HZ_DEL (0.05) HZ_NBIS (8) HZ_MARGIN (0.02 deg) HZ_PLANE_GUARD (1)
//      HZ_NFACES (0 = all) HZ_SEED (1) HZ_EPOCH_LO/HI (all) HZ_EARTH_R (6378.137); PTM_ILLUM_TRUE_NORMALS as the engine.
// env (new):
//      HZ_WRITE_HZC=1      also write out_prefix.hzc (the per-site envelope, keyed)
//      HZ_ROUNDTRIP=1      write .hzc, read it back (key verified), re-classify from the file; FATAL (rc 5) if any
//                          face-epoch gained a certification (16-bit codes round toward STRADDLE)
//      HZ_BUILD_ONLY=1     with HZ_WRITE_HZC: stream the envelope tile by tile (bounded RAM), print the census,
//                          write no per-face-epoch classes (AGT scale)
//      HZ_CENSUS_ONLY=1    any mode: print the census, keep/write no per-face-epoch classes
//      HZ_TILE=n           faces per tile (32768)
//      HZ_CLASSIFY_FROM=f  do not build: load envelope f (key must match this scene + predicate) and classify
//                          the track with it -- the cross-window reuse path; fall-throughs are counted
//      HZ_ALLOW_NONTRUE_NORMALS=1  permit a build without PTM_ILLUM_TRUE_NORMALS (else FATAL: the production
//                          chord runs with true normals; a silent winding_majority mismatch is what this guards)
#include "illumination/horizon_classifier.h"
#include "illumination/horizon_classifier_scene.h"
#include "mesh/face_orientation.h"

#include <omp.h>
#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace ptm::hzc;

static double envd(const char* k, double d) { const char* v = std::getenv(k); return v && *v ? std::atof(v) : d; }
static long envl(const char* k, long d) { const char* v = std::getenv(k); return v && *v ? std::atol(v) : d; }

// Binary little-endian PLY, double x/y/z first, uchar/int32 triangle list (the tool's reader).
static void readPly(const std::string& p, std::vector<V3>& V, std::vector<F3>& F) {
    std::ifstream f(p, std::ios::binary);
    if (!f) { std::cerr << "FATAL: cannot open " << p << "\n"; std::exit(2); }
    std::string line; size_t nv = 0, nf = 0; int cur = 0;
    std::vector<std::pair<std::string, int>> vprops;
    std::vector<int> fextra;
    bool list_first = true, seen_list = false;
    auto sz = [](const std::string& t) {
        if (t == "double" || t == "float64") return 8;
        if (t == "float" || t == "float32" || t == "int" || t == "int32" || t == "uint" || t == "uint32") return 4;
        if (t == "uchar" || t == "uint8" || t == "char" || t == "int8") return 1;
        if (t == "short" || t == "ushort" || t == "int16" || t == "uint16") return 2;
        return -1; };
    bool binle = false;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ss(line); std::string w; ss >> w;
        if (w == "format") { std::string fm; ss >> fm; binle = (fm == "binary_little_endian"); }
        else if (w == "element") { std::string e; size_t n; ss >> e >> n;
            if (e == "vertex") { nv = n; cur = 1; } else if (e == "face") { nf = n; cur = 2; } else cur = 3; }
        else if (w == "property") { std::string t; ss >> t;
            if (cur == 1) { std::string nm; ss >> nm; vprops.push_back({nm, sz(t)}); }
            else if (cur == 2) {
                if (t == "list") { std::string a, b; ss >> a >> b;
                    if (sz(a) != 1 || sz(b) != 4) { std::cerr << "FATAL: unsupported face list in " << p << "\n"; std::exit(2); }
                    seen_list = true; }
                else { std::string nm; ss >> nm; if (!seen_list) list_first = false; else fextra.push_back(sz(t)); }
            } }
        else if (w == "end_header") break;
    }
    if (!binle) { std::cerr << "FATAL: " << p << " is not binary_little_endian PLY\n"; std::exit(2); }
    if (!list_first) { std::cerr << "FATAL: face props before list unsupported\n"; std::exit(2); }
    int vstride = 0, ox = -1, oy = -1, oz = -1;
    for (auto& pr : vprops) { if (pr.first == "x") ox = vstride; if (pr.first == "y") oy = vstride;
        if (pr.first == "z") oz = vstride; vstride += pr.second; }
    if (vprops.empty() || vprops[0].second != 8 || ox < 0 || oy < 0 || oz < 0) { std::cerr << "FATAL: expect double x y z in " << p << "\n"; std::exit(2); }
    std::vector<char> buf(vstride);
    V.resize(nv);
    for (size_t i = 0; i < nv; ++i) { f.read(buf.data(), vstride);
        std::memcpy(&V[i][0], buf.data() + ox, 8); std::memcpy(&V[i][1], buf.data() + oy, 8); std::memcpy(&V[i][2], buf.data() + oz, 8); }
    const int extra = std::accumulate(fextra.begin(), fextra.end(), 0);
    F.resize(nf);
    for (size_t i = 0; i < nf; ++i) { unsigned char n; f.read((char*)&n, 1);
        if (n != 3) { std::cerr << "FATAL: non-triangle face in " << p << "\n"; std::exit(2); }
        f.read((char*)F[i].data(), 12); if (extra) f.ignore(extra); }
    if (!f) { std::cerr << "FATAL: short read " << p << "\n"; std::exit(2); }
}

static void wr(const std::string& p, const void* d, size_t b) {
    std::ofstream o(p, std::ios::binary); o.write((const char*)d, (std::streamsize)b); o.close();
    if (!o) { std::cerr << "FATAL: write failed: " << p << "\n"; std::exit(3); }
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: ptm_hzcls core.ply ring.ply|- track.bin out_prefix [threads]\n"; return 1; }
    const std::string corep = argv[1], ringp = argv[2], trackp = argv[3], outp = argv[4];
    if (argc > 5) omp_set_num_threads(std::atoi(argv[5]));
    Params P;
    P.dc_deg = envd("HZ_DC", 0.125); P.del_deg = envd("HZ_DEL", 0.05); P.nbis = (int)envl("HZ_NBIS", 8);
    P.margin_deg = envd("HZ_MARGIN", 0.02); P.plane_guard = envl("HZ_PLANE_GUARD", 1) != 0;
    P.earth_radius_km = envd("HZ_EARTH_R", 6378.137);
    const size_t NSAMP = envl("HZ_NFACES", 0); const unsigned SEED = (unsigned)envl("HZ_SEED", 1);
    const char* from = std::getenv("HZ_CLASSIFY_FROM");
    const bool roundtrip = envl("HZ_ROUNDTRIP", 0) != 0;
    const bool write_hzc = roundtrip || envl("HZ_WRITE_HZC", 0) != 0;
    const bool tn = ptm::orient::trueNormalsRequestedFromEnv();
    if (!tn && envl("HZ_ALLOW_NONTRUE_NORMALS", 0) == 0) {
        std::cerr << "FATAL: PTM_ILLUM_TRUE_NORMALS is not set. The production chord runs with true normals; an envelope "
                     "built with the winding_majority rule certifies a different predicate. Set it, or set "
                     "HZ_ALLOW_NONTRUE_NORMALS=1 to build one deliberately (the key records the rule either way).\n";
        return 1;
    }

    std::vector<V3> Vc, Vr; std::vector<F3> Fc, Fr;
    readPly(corep, Vc, Fc); if (ringp != "-") readPly(ringp, Vr, Fr);
    ptm::orient::FaceOrientation ori = ptm::orient::computeFaceOrientation(Vc, Fc, tn);
    std::printf("orientation: requested %s applied %s reversed %zu/%zu\n", ori.rule_requested.c_str(), ori.rule.c_str(),
                ori.n_reversed, ori.n_faces);
    ClassifierScene S(Vc, Fc, Vr, Fr, ori.sign, ori.rule);
    const SceneFacts& sf = S.facts();
    std::printf("core %zu faces (%zu boundary edges, sphere %s), ring %zu, curtain %zu, R_occ %.2f\n", S.numCoreFaces(),
                S.numBoundaryEdges(), sf.sphere_on ? "on" : "off", S.numRingFaces(), S.numCurtainFaces(), sf.R_occ);

    std::ifstream tf(trackp, std::ios::binary); int64_t NE = 0; tf.read((char*)&NE, 8);
    std::vector<double> T(7 * NE); tf.read((char*)T.data(), 8 * 7 * NE);
    if (!tf || NE <= 0) { std::cerr << "FATAL: short track " << trackp << "\n"; return 2; }
    const int64_t E0 = envl("HZ_EPOCH_LO", 0), E1 = std::min<int64_t>(envl("HZ_EPOCH_HI", NE), NE); const int64_t NEu = E1 - E0;
    if (E0 < 0 || NEu <= 0) { std::cerr << "FATAL: empty epoch range\n"; return 2; }
    std::vector<Epoch> ep(NEu);
    for (int64_t e = 0; e < NEu; ++e) { const double* tr = &T[7 * (E0 + e)];
        ep[e].et = tr[0]; ep[e].sun_m = {tr[1] * 1000, tr[2] * 1000, tr[3] * 1000}; ep[e].earth_m = {tr[4] * 1000, tr[5] * 1000, tr[6] * 1000}; }

    HzcReader RD;
    size_t nf = Fc.size(); std::vector<size_t> idx(nf); std::iota(idx.begin(), idx.end(), 0);
    if (from) {
        // The file's params win (they are part of what was certified); the run's key must match.
        HzcReader probe; std::string err = probe.open(from, nullptr, false);
        if (!err.empty()) { std::cerr << "FATAL: " << err << "\n"; return 4; }
        const Key want = S.key(probe.params());
        err = RD.open(from, &want, true);
        if (!err.empty()) { std::cerr << "FATAL: " << err << "\n"; return 4; }
        P = RD.params();
        idx.resize(RD.numFaces()); for (uint64_t i = 0; i < RD.numFaces(); ++i) idx[i] = RD.faceAt(i);
        std::printf("envelope %s: key scene %016" PRIx64 " predicate %016" PRIx64 " OK; %zu faces; note: %s\n", from,
                    RD.key().scene, RD.key().predicate, idx.size(), RD.note().c_str());
    } else if (NSAMP && NSAMP < nf) {
        std::mt19937_64 rng(SEED); std::shuffle(idx.begin(), idx.end(), rng); idx.resize(NSAMP); std::sort(idx.begin(), idx.end());
    }
    const size_t NS = idx.size();
    const Key key_here = S.key(P);   // after P is final (classify mode takes the file's params)
    const bool build_only = envl("HZ_BUILD_ONLY", 0) != 0;
    const bool census_only = build_only || envl("HZ_CENSUS_ONLY", 0) != 0;   // no per-face-epoch classes kept or written
    const size_t TILE = (size_t)std::max<long>(1, envl("HZ_TILE", 32768));
    if (build_only && (from || !write_hzc)) { std::cerr << "FATAL: HZ_BUILD_ONLY needs HZ_WRITE_HZC=1 and no HZ_CLASSIFY_FROM\n"; return 1; }
    if (census_only && roundtrip) { std::cerr << "FATAL: HZ_ROUNDTRIP needs the per-face-epoch classes (not HZ_BUILD_ONLY)\n"; return 1; }
    std::printf("faces %zu, epochs %" PRId64 " [%" PRId64 ",%" PRId64 "), DC %.4f DEL %.4f NBIS %d margin %.4f deg, plane guard %d, %s, tile %zu\n",
                NS, NEu, E0, E1, P.dc_deg, P.del_deg, P.nbis, P.margin_deg, (int)P.plane_guard,
                from ? "classify from envelope" : (build_only ? "build only (census, no per-face-epoch output)" : "build"), TILE);

    // Per-face-epoch classes are kept only when they are written (small runs); a build-only
    // run holds one tile of envelopes at a time and streams it to the .hzc.
    std::vector<uint8_t> cls(census_only ? 0 : NS * NEu);
    std::vector<uint64_t> env_rays(NS, 0); std::vector<uint32_t> nmulti(NS, 0);
    ClassCensus census; uint64_t nc = 0;
    HzcWriter W;
    if (!from && write_hzc) {
        char note[256]; std::snprintf(note, sizeof note, "ptm_hzcls build: %" PRId64 " epochs ET %.3f..%.3f", NEu, ep.front().et, ep.back().et);
        const std::string err = W.open(outp + ".hzc", key_here, P, nf, note);
        if (!err.empty()) { std::cerr << "FATAL: " << err << "\n"; return 3; }
    }
    auto t0 = std::chrono::steady_clock::now();
    std::vector<FaceEnvelope> tile_env;
    for (size_t s0 = 0; s0 < NS; s0 += TILE) {
        const size_t s1 = std::min(NS, s0 + TILE);
        tile_env.assign(s1 - s0, FaceEnvelope());
        #pragma omp parallel
        {
            ClassCensus local; std::vector<EpochGeom> g(NEu);
            #pragma omp for schedule(dynamic, 16)
            for (size_t s = s0; s < s1; ++s) {
                const size_t fi = idx[s];
                const FaceFrame fr = S.frame(fi);
                for (int64_t e = 0; e < NEu; ++e) g[e] = epochGeom(fr, ep[e], sf, P);
                FaceEnvelope& E = tile_env[s - s0];
                if (from) {
                    if (!RD.get((uint32_t)fi, E)) E.cols.clear();   // unreadable -> every epoch falls through
                } else {
                    E = buildEnvelope(fr, g, sf, P, [&](const V3& d) { return S.occluded(fi, fr, d); }, &env_rays[s]);
                    for (const Column& c : E.cols) nmulti[s] += c.multi;
                }
                for (int64_t e = 0; e < NEu; ++e) {
                    uint8_t r; const uint8_t k = classify(E, g[e], P, &r); local.add(k, r);
                    if (!census_only) cls[s * NEu + e] = k;
                }
            }
            #pragma omp critical
            census.merge(local);
        }
        for (size_t s = s0; s < s1; ++s) {
            nc += tile_env[s - s0].cols.size();
            if (!from && write_hzc) {
                const std::string err = W.append((uint32_t)idx[s], tile_env[s - s0]);
                if (!err.empty()) { std::cerr << "FATAL: " << err << "\n"; return 3; }
            }
        }
        if (census_only && (s0 / TILE) % 8 == 0)
            std::printf("  tile %zu/%zu: %zu faces done, %.0f s, %.2f GB written\n", s0 / TILE + 1, (NS + TILE - 1) / TILE, s1,
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), W.bytesWritten() / 1e9);
    }
    std::vector<FaceEnvelope>().swap(tile_env);
    if (!from && write_hzc) {
        const std::string err = W.close();
        if (!err.empty()) { std::cerr << "FATAL: " << err << "\n"; return 3; }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    uint64_t er = 0, nm = 0;
    for (size_t s = 0; s < NS; ++s) { er += env_rays[s]; nm += nmulti[s]; }
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    const double rss_gb = ru.ru_maxrss / 1e9;      // bytes on macOS
#else
    const double rss_gb = ru.ru_maxrss * 1024.0 / 1e9;   // KiB on Linux
#endif
    std::printf("wall %.1f s; envelope Embree rays %" PRIu64 " (%.0f /face; %.1f cols/face, multi cols %.4f%%); peak RSS %.2f GB\n",
                secs, er, (double)er / NS, (double)nc / NS, 100.0 * nm / std::max<uint64_t>(1, nc), rss_gb);
    std::printf("%s\n", census.summary().c_str());
    if (census.fallthrough())
        std::printf("NOTE: %" PRIu64 " face-epochs fell through to the exact chord (outside the built ladder or an unvisited "
                    "column); they cost chord rays, not accuracy\n", census.fallthrough());
    std::printf("key scene %016" PRIx64 " predicate %016" PRIx64 "\n", key_here.scene, key_here.predicate);
    if (!from && write_hzc)
        std::printf("wrote %s.hzc: %.3f GB, %" PRIu64 " columns (%.2f B/column incl. runs and index)\n", outp.c_str(),
                    W.bytesWritten() / 1e9, W.columns(), (double)W.bytesWritten() / std::max<uint64_t>(1, W.columns()));

    { std::ofstream o(outp + ".meta"); if (!o) { std::cerr << "FATAL: cannot open " << outp << ".meta\n"; return 3; }
      o << "nfaces " << NS << "\nnepochs " << NEu << "\nepoch_lo " << E0 << "\nDC " << P.dc_deg << "\nDEL " << P.del_deg
        << "\nNBIS " << P.nbis << "\nMARGIN_deg " << P.margin_deg << "\nPLANE_GUARD " << P.plane_guard << "\nwall_s " << secs
        << "\nenv_rays " << er << "\npeak_rss_gb " << rss_gb << "\nkey_scene " << ptm::orient::hex16(key_here.scene)
        << "\nkey_predicate " << ptm::orient::hex16(key_here.predicate) << "\nfallthrough " << census.fallthrough() << "\n";
      o.close(); if (!o) { std::cerr << "FATAL: write failed: " << outp << ".meta\n"; return 3; } }
    std::vector<uint32_t> id32(idx.begin(), idx.end());
    wr(outp + ".faces.u32", id32.data(), 4 * NS);
    if (!census_only) wr(outp + ".cls.u8", cls.data(), cls.size());
    if (!from) { wr(outp + ".envrays.u64", env_rays.data(), 8 * NS); wr(outp + ".nmulti.u32", nmulti.data(), 4 * NS); }

    if (roundtrip) {
        // Conservativeness gate: classes from the quantized file must be a SUBSET of the double
        // envelope's certifications -- a face-epoch may lose CLEAR/BLOCKED (-> STRADDLE), never gain one.
        HzcReader R; std::string err = R.open(outp + ".hzc", &key_here, true);
        if (!err.empty()) { std::cerr << "FATAL: round trip: " << err << "\n"; return 4; }
        uint64_t changed = 0, unsafe = 0, missing = 0; ClassCensus rc;
        #pragma omp parallel
        {
            ClassCensus local; uint64_t ch = 0, un = 0, mi = 0; std::vector<EpochGeom> g(NEu); FaceEnvelope E;
            #pragma omp for schedule(dynamic, 16)
            for (size_t s = 0; s < NS; ++s) {
                const FaceFrame fr = S.frame(idx[s]);
                if (!R.get((uint32_t)idx[s], E)) { ++mi; E.cols.clear(); }
                for (int64_t e = 0; e < NEu; ++e) {
                    g[e] = epochGeom(fr, ep[e], sf, R.params()); uint8_t r;
                    const uint8_t k = classify(E, g[e], R.params(), &r); local.add(k, r);
                    const uint8_t k0 = cls[s * NEu + e];
                    if (k != k0) { ++ch; if (k != STRADDLE) ++un; }
                }
            }
            #pragma omp critical
            { rc.merge(local); changed += ch; unsafe += un; missing += mi; }
        }
        std::printf("round trip (16-bit codes): %" PRIu64 " face-epochs changed class, %" PRIu64 " NOT to STRADDLE, %" PRIu64
                    " faces missing\n%s\n", changed, unsafe, missing, rc.summary().c_str());
        if (unsafe || missing) { std::cerr << "FATAL: the quantized envelope certified a face-epoch the double envelope did not\n"; return 5; }
    }
    return 0;
}
