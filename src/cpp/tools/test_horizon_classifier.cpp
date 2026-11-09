// test_horizon_classifier -- gate for the horizon-classifier library (history marcher M0).
//
// Part A (always; synthetic, no data): a flat plate on a sphere of radius 1 km with a 20 m
// radial wall 30 m north of the test face.
//   A1  every CLEAR / BLOCKED certification on the build track is confirmed by brute force
//       (dense disk sampling with the same predicate), and the build track never falls through;
//   A2  a track outside the build (higher/lower elevations, unvisited azimuths) is never certified:
//       every such face-epoch falls through, and the census says why;
//   A3  a disk crossing the face plane always straddles (rule v2 guard);
//   A4  .hzc (v2, 16-bit codes): round trip (only -> STRADDLE), key mismatch on changed params, on a 1 mm vertex
//       move at the SAME triangle count (PTM-800), and a corrupt byte, all refused;
//   A5  audit selection: deterministic, rate 0.1% within 5 sigma; AuditCensus flags wrong classes;
//   A6  16-bit elevation codes: always rounded toward STRADDLE (up / down as asked), minimally, sentinels kept.
// Part B (--parity REF NEW): byte comparison of the measurement tool's outputs (REF prefix) with
//   ptm_hzcls's (NEW prefix): .faces.u32 .cls.u8 .envrays.u64 .nmulti.u32. Exit 1 on any difference,
//   with the class confusion matrix.
#include "illumination/horizon_classifier.h"
#include "illumination/horizon_classifier_scene.h"
#include "mesh/face_orientation.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

using namespace ptm::hzc;

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { ++g_fail; std::printf("  FAIL %s:%d: ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static std::vector<char> slurp(const std::string& p, bool& ok) {
    std::ifstream f(p, std::ios::binary); ok = (bool)f;
    return std::vector<char>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

static int parity(const std::string& ref, const std::string& neu) {
    int bad = 0;
    std::vector<char> cr, cn;
    for (const char* ext : {".faces.u32", ".cls.u8", ".envrays.u64", ".nmulti.u32"}) {
        bool a, b; auto R = slurp(ref + ext, a); auto N = slurp(neu + ext, b);
        if (!a || !b) { std::printf("  %-12s MISSING (%s %s)\n", ext, a ? "" : (ref + ext).c_str(), b ? "" : (neu + ext).c_str()); ++bad; continue; }
        size_t nd = 0; const size_t n = std::min(R.size(), N.size());
        for (size_t i = 0; i < n; ++i) nd += (R[i] != N[i]);
        const bool same = R.size() == N.size() && nd == 0;
        std::printf("  %-12s %s (%zu vs %zu bytes, %zu differ)\n", ext, same ? "IDENTICAL" : "DIFFERENT", R.size(), N.size(), nd);
        if (!same) ++bad;
        if (std::strcmp(ext, ".cls.u8") == 0) { cr = R; cn = N; }
    }
    if (!cr.empty() && cr.size() == cn.size()) {
        unsigned long long M[6][6] = {};
        for (size_t i = 0; i < cr.size(); ++i) { const int a = (uint8_t)cr[i], b = (uint8_t)cn[i]; if (a < 6 && b < 6) ++M[a][b]; }
        const char* nm[6] = {"free-plane", "free-cone", "CLEAR", "BLOCKED", "STRADDLE", "ECLIPSE"};
        std::printf("  class counts (ref -> new, off-diagonal = disagreement):\n");
        for (int a = 0; a < 6; ++a) { std::printf("    %-10s", nm[a]); for (int b = 0; b < 6; ++b) std::printf(" %11llu", M[a][b]); std::printf("\n"); }
    }
    std::printf("PARITY %s\n", bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}

// ---- synthetic scene ----------------------------------------------------------------------
struct Synth {
    std::vector<V3> V; std::vector<F3> F; size_t test_face = 0;
    explicit Synth(double wall_shift_m = 0.0) {
        const double R = 1000.0; const int N = 10; const double L = 50.0;
        // plate: x = R, y,z in [-L, L], normal +x (winding chosen so cross(B-A, C-A) = +x)
        for (int i = 0; i <= N; ++i) for (int j = 0; j <= N; ++j) V.push_back({R, -L + 2 * L * i / N, -L + 2 * L * j / N});
        auto id = [&](int i, int j) { return i * (N + 1) + j; };
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
            F.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            F.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
        // the face nearest the plate center
        test_face = 2 * ((N / 2) * N + (N / 2));
        // wall: z = 30 (+ shift), y in [-60, 60], x in [R, R + 20]
        const double zw = 30.0 + wall_shift_m; const int b = (int)V.size();
        V.push_back({R, -60, zw}); V.push_back({R, 60, zw}); V.push_back({R + 20, 60, zw}); V.push_back({R + 20, -60, zw});
        F.push_back({b, b + 1, b + 2}); F.push_back({b, b + 2, b + 3});
    }
};

static Epoch sunAt(const FaceFrame& f, double az_deg, double el_deg) {
    const double az = az_deg * M_PI / 180, el = el_deg * M_PI / 180, D = 1.495978707e11;
    const V3 d = {std::cos(el) * (std::cos(az) * f.north[0] + std::sin(az) * f.east[0]) + std::sin(el) * f.up[0],
                  std::cos(el) * (std::cos(az) * f.north[1] + std::sin(az) * f.east[1]) + std::sin(el) * f.up[1],
                  std::cos(el) * (std::cos(az) * f.north[2] + std::sin(az) * f.east[2]) + std::sin(el) * f.up[2]};
    Epoch e; e.et = az_deg * 1000 + el_deg;
    e.sun_m = {f.fc[0] + D * d[0], f.fc[1] + D * d[1], f.fc[2] + D * d[2]};
    e.earth_m = {f.fc[0] - 3.84e8 * d[0], f.fc[1] - 3.84e8 * d[1], f.fc[2] - 3.84e8 * d[2]};   // behind: never eclipses
    return e;
}

// Brute-force visible fraction of the disk (uniform grid in the disk, same predicate).
static double bruteFrac(const ClassifierScene& S, size_t fi, const FaceFrame& f, const EpochGeom& q) {
    const V3 sv = q.s; const double us = sv[0] * f.up[0] + sv[1] * f.up[1] + sv[2] * f.up[2];
    V3 v = {f.up[0] - us * sv[0], f.up[1] - us * sv[1], f.up[2] - us * sv[2]};
    const double lv = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); v = {v[0] / lv, v[1] / lv, v[2] / lv};
    const V3 u = {v[1] * sv[2] - v[2] * sv[1], v[2] * sv[0] - v[0] * sv[2], v[0] * sv[1] - v[1] * sv[0]};
    const double ta = std::tan(q.alpha); int n = 0, clear = 0; const int G = 24;
    for (int a = 0; a < G; ++a) for (int b = 0; b < G; ++b) {
        const double x = -1 + 2 * (a + 0.5) / G, y = -1 + 2 * (b + 0.5) / G; if (x * x + y * y > 1) continue;
        V3 d = {sv[0] + ta * (x * u[0] + y * v[0]), sv[1] + ta * (x * u[1] + y * v[1]), sv[2] + ta * (x * u[2] + y * v[2])};
        const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]); d = {d[0] / l, d[1] / l, d[2] / l};
        ++n;
        if (d[0] * f.fn[0] + d[1] * f.fn[1] + d[2] * f.fn[2] > 0 && !S.occluded(fi, f, d)) ++clear;
    }
    return n ? (double)clear / n : 0.0;
}

static void partA() {
    std::printf("Part A (synthetic)\n");
    Synth G;
    auto ori = ptm::orient::computeFaceOrientation(G.V, G.F, false);
    ClassifierScene S(G.V, G.F, {}, {}, ori.sign, ori.rule);
    const size_t fi = G.test_face; const FaceFrame f = S.frame(fi);
    CHECK(f.fn[0] > 0.999, "test face normal not +x (%g %g %g)", f.fn[0], f.fn[1], f.fn[2]);
    Params P; P.dc_deg = 0.5; P.del_deg = 0.25; P.nbis = 8; P.margin_deg = 0.1;
    const SceneFacts& sf = S.facts();

    // A1: build track = az {0 (north, wall), 180 (south, open)} x el 3..60 deg step 0.5
    std::vector<Epoch> A;
    for (double az : {0.0, 180.0}) for (double el = 3.0; el <= 60.0; el += 0.5) A.push_back(sunAt(f, az, el));
    std::vector<EpochGeom> g; for (auto& e : A) g.push_back(epochGeom(f, e, sf, P));
    uint64_t rays = 0;
    const FaceEnvelope env = buildEnvelope(f, g, sf, P, [&](const V3& d) { return S.occluded(fi, f, d); }, &rays);
    ClassCensus cA; int ncert = 0, nbad = 0;
    for (auto& q : g) {
        uint8_t r; const uint8_t k = classify(env, q, P, &r); cA.add(k, r);
        if (k == CLEAR || k == BLOCKED) {
            ++ncert; const double fr = bruteFrac(S, fi, f, q);
            if ((k == CLEAR && fr < 1.0) || (k == BLOCKED && fr > 0.0)) ++nbad;
        }
    }
    std::printf("  A1 build track: %s\n", cA.summary().c_str());
    CHECK(cA.fallthrough() == 0, "build track fell through %llu times", (unsigned long long)cA.fallthrough());
    CHECK(cA.by_class[CLEAR] > 0 && cA.by_class[BLOCKED] > 0 && cA.by_class[STRADDLE] > 0,
          "expected all of CLEAR/BLOCKED/STRADDLE on the build track");
    CHECK(nbad == 0, "%d of %d certifications contradicted by brute force", nbad, ncert);
    std::printf("  A1 %d certifications, %d contradicted by brute force; %llu envelope rays\n", ncert, nbad, (unsigned long long)rays);

    // A2: outside the build: el 70..85 (above every band), az 90 / 270 (never visited at low el)
    ClassCensus cB; int cert_out = 0;
    for (double az : {0.0, 90.0, 180.0, 270.0}) for (double el = (az == 0.0 || az == 180.0) ? 70.0 : 5.0; el <= 85.0; el += 1.0) {
        const EpochGeom q = epochGeom(f, sunAt(f, az, el), sf, P); uint8_t r; const uint8_t k = classify(env, q, P, &r);
        cB.add(k, r); if (k == CLEAR || k == BLOCKED) ++cert_out;
    }
    // below the built band in visited columns (the build started at el 3)
    for (double az : {0.0, 180.0}) for (double el = 1.0; el <= 2.5; el += 0.5) {
        const EpochGeom q = epochGeom(f, sunAt(f, az, el), sf, P); uint8_t r; const uint8_t k = classify(env, q, P, &r);
        cB.add(k, r); if (k == CLEAR || k == BLOCKED) ++cert_out;
    }
    std::printf("  A2 foreign track: %s\n", cB.summary().c_str());
    CHECK(cert_out == 0, "%d foreign face-epochs certified", cert_out);
    CHECK(cB.by_reason[R_UNVISITED] > 0 && cB.by_reason[R_OUT_OF_LADDER] > 0, "both fall-through reasons expected");

    // A3: plane guard -- sun at el 0.1 deg south (disk crosses the plate plane)
    {
        const EpochGeom q = epochGeom(f, sunAt(f, 180.0, 0.1), sf, P);
        std::vector<EpochGeom> g3{q}; uint64_t r3 = 0;
        const FaceEnvelope e3 = buildEnvelope(f, g3, sf, P, [&](const V3& d) { return S.occluded(fi, f, d); }, &r3);
        uint8_t r; const uint8_t k = classify(e3, q, P, &r);
        Params P1 = P; P1.plane_guard = false; uint8_t r1; const uint8_t k1 = classify(e3, q, P1, &r1);
        std::printf("  A3 disk across the face plane: v2 class %d reason %d; v1 class %d\n", k, r, k1);
        CHECK(k == STRADDLE && r == R_PLANE_GUARD, "plane guard did not force STRADDLE (class %d reason %d)", k, r);
    }

    // A4: .hzc
    {
        const std::string p = "test_horizon_classifier.tmp.hzc";
        HzcFile H; H.key = S.key(P); H.params = P; H.nfaces_scene = G.F.size(); H.faces = {(uint32_t)fi}; H.env = {env};
        H.build_note = "synthetic";
        CHECK(writeHzc(p, H).empty(), "write failed");
        HzcFile R; const Key k = S.key(P);
        std::string err = readHzc(p, R, &k);
        CHECK(err.empty(), "read with the right key failed: %s", err.c_str());
        int changed = 0, unsafe = 0;
        if (err.empty()) for (auto& q : g) {
            const uint8_t a = classify(env, q, P, nullptr), b = classify(R.env[0], q, R.params, nullptr);
            if (a != b) { ++changed; if (b != STRADDLE) ++unsafe; }
        }
        std::printf("  A4 round trip: %d changed, %d not to STRADDLE\n", changed, unsafe);
        CHECK(unsafe == 0, "round trip certified something the double envelope did not");
        Params P2 = P; P2.margin_deg = 0.05; const Key k2 = S.key(P2);
        err = readHzc(p, R, &k2); CHECK(!err.empty(), "changed params were accepted");
        std::printf("  A4 changed margin -> %s\n", err.empty() ? "ACCEPTED" : "refused");
        Synth G2(0.001);   // move the wall by 1 mm: same vertex and triangle counts
        auto o2 = ptm::orient::computeFaceOrientation(G2.V, G2.F, false);
        ClassifierScene S2(G2.V, G2.F, {}, {}, o2.sign, o2.rule);
        const Key k3 = S2.key(P);
        CHECK(G2.F.size() == G.F.size() && G2.V.size() == G.V.size(), "synthetic counts differ");
        err = readHzc(p, R, &k3); CHECK(!err.empty(), "1 mm-moved scene with the same triangle count was accepted");
        std::printf("  A4 1 mm wall move, same triangle count -> %s\n", err.empty() ? "ACCEPTED" : "refused");
        auto flip = [&](bool middle) {   // flip one bit: a payload byte (middle) or a footer offset byte
            std::fstream io(p, std::ios::binary | std::ios::in | std::ios::out); io.seekg(0, std::ios::end); const auto n = io.tellg();
            const std::streamoff at = middle ? (std::streamoff)n / 2 : (std::streamoff)n - 20;
            char c; io.seekg(at); io.read(&c, 1); c ^= 0x40; io.seekp(at); io.write(&c, 1);
        };
        flip(true);
        err = readHzc(p, R, &k); CHECK(!err.empty(), "corrupt payload byte was accepted");
        std::printf("  A4 corrupt payload byte -> %s\n", err.empty() ? "ACCEPTED" : ("refused: " + err.substr(err.find(':') + 2)).c_str());
        flip(true); flip(false);
        err = readHzc(p, R, &k); CHECK(!err.empty(), "corrupt footer byte was accepted");
        std::printf("  A4 corrupt footer byte -> %s\n", err.empty() ? "ACCEPTED" : ("refused: " + err.substr(err.find(':') + 2)).c_str());
        std::remove(p.c_str());
    }

    // A6: 16-bit elevation codes round toward STRADDLE, minimally, and keep the sentinels
    {
        uint64_t n = 0, bad_dir = 0, bad_min = 0; std::mt19937_64 rng(7);
        std::uniform_real_distribution<double> U(-M_PI / 2 - 0.05, M_PI / 2 + 0.05);
        std::vector<double> xs = {-M_PI / 2, M_PI / 2, 0.0, -0.0, 1e-300, -1e-300, qDecode(1), qDecode(2), qDecode(kQMaxFinite),
                                  std::nextafter(qDecode(100), 0.0), std::nextafter(qDecode(100), 1.0), 89.9 * M_PI / 180};
        for (int i = 0; i < 1000000; ++i) xs.push_back(U(rng));
        for (double x : xs) {
            ++n;
            const uint16_t u = qUp(x), d = qDown(x);
            if (!(qDecode(u) >= x) || !(qDecode(d) <= x)) ++bad_dir;
            if (u > 1 && u <= kQMaxFinite && qDecode((uint16_t)(u - 1)) >= x) ++bad_min;
            if (d >= 1 && d < kQMaxFinite && qDecode((uint16_t)(d + 1)) <= x) ++bad_min;
        }
        const bool sent = qUp(-1e9) == kQNegInf && qDown(1e9) == kQPosInf && qUp(1e9) == kQPosInf && qDown(-1e9) == kQNegInf;
        std::printf("  A6 codes: %llu values, %llu wrong-direction, %llu non-minimal; sentinels %s; step %.3g deg\n",
                    (unsigned long long)n, (unsigned long long)bad_dir, (unsigned long long)bad_min, sent ? "OK" : "WRONG",
                    (qDecode(2) - qDecode(1)) * 180 / M_PI);
        CHECK(bad_dir == 0 && bad_min == 0 && sent, "16-bit elevation codes are not conservative/minimal");
    }

    // A5: audit
    {
        const long N = 2000000; long n = 0, rep = 0;
        for (long i = 0; i < N; ++i) {
            const uint32_t face = (uint32_t)(i % 100003); const double et = 1e9 + 5.906 * 3600 * (double)(i / 100003);
            const bool a = auditPick(face, et); n += a; rep += (a == auditPick(face, et));
        }
        const double mean = N * kAuditRate, sd = std::sqrt(mean * (1 - kAuditRate));
        std::printf("  A5 audit picks %ld of %ld (expect %.0f +- %.0f), deterministic %ld/%ld\n", n, N, mean, sd, rep, N);
        CHECK(std::fabs(n - mean) < 5 * sd, "audit rate off");
        CHECK(rep == N, "audit selection not deterministic");
        AuditCensus ac; ac.add(CLEAR, 1.0); ac.add(BLOCKED, 0.0); ac.add(CLEAR, 0.9); ac.add(BLOCKED, 0.01);
        CHECK(ac.wrong_clear == 1 && ac.wrong_blocked == 1 && !ac.clean(), "AuditCensus did not flag");
        std::printf("  A5 %s\n", ac.summary().c_str());
    }
}

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--parity") {
        std::printf("Part B: parity %s (measurement tool) vs %s (ptm_hzcls)\n", argv[2], argv[3]);
        return parity(argv[2], argv[3]);
    }
    partA();
    std::printf("%s (%d failure%s)\n", g_fail ? "FAIL" : "PASS", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
