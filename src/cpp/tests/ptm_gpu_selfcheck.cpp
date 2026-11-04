// ptm_gpu_selfcheck -- is the GPU shadow-ray tracer right on THIS hardware?
//
// GPU illumination (illumination.use_gpu) was validated against CPU Embree on
// one GPU model (AMD Radeon Pro W6800X, macOS, 2026-11-03; PTM-794). A
// different GPU, driver or OpenCL compiler can still mis-execute the kernel:
// the macOS OpenCL CPU device does (PTM-797). This program answers the
// question for the devices it finds, in a few seconds, with no data files.
//
// What it does, per enumerated OpenCL GPU:
//   1. builds two synthetic crater meshes in memory (one with ~20 m facets,
//      one with ~2 mm facets, the two scales the lunar work spans) and a
//      259,200-triangle crater whose BVH is deep (>= 17 levels);
//   2. traces shadow rays toward 8 Sun directions from 0.5 to 60 degrees
//      elevation through the PRODUCTION entry point
//      (BVHRayTracer::traceShadowRaysEx, production BVH settings), in the ray
//      classes production uses: rays lifted off the face, and at-face rays
//      that skip their own facet and apply the rise-epsilon filter (with
//      h_eps = 0.2 facet, and with h_eps = 0);
//   3. compares every answer with a float64 brute-force ground truth computed
//      on exactly the float32 geometry and rays the device was given.
// Rays whose float64 answer is within a small tolerance of a tie (an edge
// grazed, a hit at the near limit) are excluded and counted: either answer is
// correct for them. Every other ray must agree. The kernel compiled for the
// host (the engine's PTM_BVH_HOST_KERNEL reference) is checked the same way;
// if THAT disagrees, the self-check itself is broken and says so.
//
// Exit status: 0 PASS (every device), 1 FAIL (a device gave a wrong answer),
// 2 SKIP (no OpenCL GPU support in this build, or no device), 3 INVALID (the
// check could not be trusted: reference disagreed, too many ties, or a
// vacuous scene). Usage: ptm_gpu_selfcheck [--device all|N] [--pause-ms N]
// [--rays-per-direction N]. Also reachable as `ptm gpu-selfcheck`.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <thread>
#include <vector>

#ifdef PTM_GPU_OPENCL
#include "bvh_opencl.h"
#include "opencl_config.h"
#endif

namespace {

using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
double f32(double x) { return static_cast<double>(static_cast<float>(x)); }

struct Scene {
    std::string name;
    double edge = 0.0;                          // nominal facet edge, m
    std::vector<Vec3> V;                        // float32-representable
    std::vector<std::array<uint32_t, 3>> F;
    std::vector<Vec3> N, C;                     // unit normal (up), centroid
};

// A bowl crater (depth/diameter 0.2) with a raised rim and deterministic
// roughness, on an n x n vertex grid. Coordinates are rounded to float32 at
// construction, so the device and the ground truth see the same geometry.
Scene makeCrater(const std::string& name, double R, int n) {
    Scene s;
    s.name = name;
    const double L = 1.6 * R;
    s.edge = 2.0 * L / (n - 1);
    auto height = [R](double x, double y) {
        const double r = std::sqrt(x * x + y * y) / R;
        double z = (r < 1.0) ? 0.4 * R * (r * r - 1.0) : 0.0;
        z += 0.08 * R * std::exp(-((r - 1.0) / 0.25) * ((r - 1.0) / 0.25));
        z += 0.01 * R * (std::sin(9.1 * x / R + 0.3) * std::cos(7.7 * y / R - 1.1) +
                         0.5 * std::sin(23.3 * (x + y) / R));
        return z;
    };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const double x = -L + 2.0 * L * i / (n - 1), y = -L + 2.0 * L * j / (n - 1);
            s.V.push_back({f32(x), f32(y), f32(height(x, y))});
        }
    auto id = [n](int i, int j) { return static_cast<uint32_t>(j * n + i); };
    for (int j = 0; j + 1 < n; ++j)
        for (int i = 0; i + 1 < n; ++i) {
            const uint32_t a = id(i, j), b = id(i + 1, j), c = id(i + 1, j + 1), d = id(i, j + 1);
            if ((i + j) % 2 == 0) { s.F.push_back({a, b, c}); s.F.push_back({a, c, d}); }
            else                  { s.F.push_back({a, b, d}); s.F.push_back({b, c, d}); }
        }
    for (auto& t : s.F) {
        Vec3 nrm = cross(sub(s.V[t[1]], s.V[t[0]]), sub(s.V[t[2]], s.V[t[0]]));
        if (nrm[2] < 0) { std::swap(t[1], t[2]); nrm = {-nrm[0], -nrm[1], -nrm[2]}; }
        const double l = norm(nrm);
        s.N.push_back({nrm[0] / l, nrm[1] / l, nrm[2] / l});
        Vec3 c{};
        for (int k = 0; k < 3; ++k) c[k] = (s.V[t[0]][k] + s.V[t[1]][k] + s.V[t[2]][k]) / 3.0;
        s.C.push_back(c);
    }
    return s;
}

// Shadow rays in the layout traceShadowRaysEx takes.
struct Rays {
    std::string cls;
    std::vector<float> ox, oy, oz, dx, dy, dz, tmin, tmax, dn, heps;
    std::vector<int32_t> skip;
    size_t size() const { return ox.size(); }
};

// The production ray classes. 0 "lifted": origin a quarter facet above the
// face, plain any-hit occlusion. 1 "at-face": origin on the face, its own
// facet skipped, hits rising less than 0.2 facet out of its plane ignored
// (PTM_ILLUM_RISE_EPS with PTM_ILLUM_RISE_KEPS=0.2). 2 "at-face0": the same
// with h_eps = 0, i.e. only hits behind the face plane are ignored
// (PTM_ILLUM_RISE_HABS=0, the validated production recipe).
Rays makeRays(const Scene& s, int cls, int per_dir, uint32_t seed) {
    const bool at_face = cls != 0;
    Rays r;
    r.cls = cls == 0 ? "lifted" : (cls == 1 ? "at-face" : "at-face0");
    const double elev[8] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 60.0};
    const double deg = std::acos(-1.0) / 180.0;
    for (int k = 0; k < 8; ++k) {
        const double e = elev[k] * deg, a = (37.0 + 97.0 * k) * deg;
        const Vec3 d = {std::cos(e) * std::cos(a), std::cos(e) * std::sin(a), std::sin(e)};
        std::vector<size_t> facing;
        for (size_t f = 0; f < s.F.size(); ++f)
            if (dot(s.N[f], d) > 0.02) facing.push_back(f);
        if (facing.empty()) continue;
        std::mt19937 g(seed + 7919u * static_cast<uint32_t>(k));
        std::uniform_int_distribution<size_t> pick(0, facing.size() - 1);
        for (int q = 0; q < per_dir; ++q) {
            const size_t f = facing[pick(g)];
            const double lift = at_face ? 0.0 : 0.25 * s.edge;
            Vec3 o;
            for (int c = 0; c < 3; ++c) o[c] = s.C[f][c] + s.N[f][c] * lift;
            r.ox.push_back(static_cast<float>(o[0]));
            r.oy.push_back(static_cast<float>(o[1]));
            r.oz.push_back(static_cast<float>(o[2]));
            r.dx.push_back(static_cast<float>(d[0]));
            r.dy.push_back(static_cast<float>(d[1]));
            r.dz.push_back(static_cast<float>(d[2]));
            r.tmin.push_back(0.0f);
            r.tmax.push_back(1.0e30f);
            r.skip.push_back(at_face ? static_cast<int32_t>(f) : -1);
            r.dn.push_back(at_face ? static_cast<float>(dot(d, s.N[f])) : 0.0f);
            r.heps.push_back(cls == 1 ? static_cast<float>(0.2 * s.edge) : (cls == 2 ? 0.0f : -1.0f));
        }
    }
    return r;
}

// Depth of a built BVH (root = level 1).
int bvhDepth(const gpu::opencl::BVHBuilder& b) {
    const auto& N = b.nodes();
    if (N.empty()) return 0;
    int best = 0;
    std::vector<std::pair<uint32_t, int>> st = {{0u, 1}};
    while (!st.empty()) {
        auto [i, d] = st.back(); st.pop_back();
        best = std::max(best, d);
        if (N[i].tri_count == 0) { st.push_back({N[i].left_child, d + 1}); st.push_back({N[i].right_child, d + 1}); }
    }
    return best;
}

// float64 brute force on the float32 inputs. 0 clear, 1 blocked, 2 near-tie
// (the answer flips within the tolerances; either answer is acceptable).
std::vector<uint8_t> groundTruth(const Scene& s, const Rays& r) {
    const double eb = 1e-4;               // barycentric tolerance
    const double et = 1e-4 * s.edge;      // distance tolerances, m
    std::vector<uint8_t> out(r.size(), 0);
    const unsigned nth = std::max(1u, std::thread::hardware_concurrency());
    std::vector<std::thread> pool;
    for (unsigned th = 0; th < nth; ++th) {
        pool.emplace_back([&, th]() {
            for (size_t i = th; i < r.size(); i += nth) {
                const Vec3 o = {r.ox[i], r.oy[i], r.oz[i]};
                const Vec3 d = {r.dx[i], r.dy[i], r.dz[i]};
                const double tmin = r.tmin[i], dn = r.dn[i], heps = r.heps[i];
                bool loose = false, strict = false;
                for (size_t j = 0; j < s.F.size() && !strict; ++j) {
                    if (static_cast<int32_t>(j) == r.skip[i]) continue;
                    const Vec3& A = s.V[s.F[j][0]];
                    const Vec3 e1 = sub(s.V[s.F[j][1]], A), e2 = sub(s.V[s.F[j][2]], A);
                    const Vec3 p = cross(d, e2);
                    const double det = dot(e1, p);
                    if (std::fabs(det) < 1e-12 * norm(cross(e1, e2))) continue;  // parallel
                    const Vec3 sv = sub(o, A);
                    const double u = dot(sv, p) / det;
                    const Vec3 q = cross(sv, e1);
                    const double v = dot(d, q) / det;
                    const double t = dot(e2, q) / det;
                    const double m = std::min({u, v, 1.0 - u - v});
                    if (m < -eb || !std::isfinite(t)) continue;
                    const bool filt = heps >= 0.0;
                    const bool hit_l = t >= tmin - et && (!filt || t * dn >= heps - et);
                    const bool hit_s = m > eb && t >= tmin + et && (!filt || t * dn >= heps + et);
                    loose = loose || hit_l;
                    strict = strict || hit_s;
                }
                out[i] = (loose == strict) ? static_cast<uint8_t>(loose) : 2;
            }
        });
    }
    for (auto& t : pool) t.join();
    return out;
}

struct Tally {
    size_t rays = 0, decided = 0, blocked = 0, ties = 0, wrong = 0, other_than_ref = 0;
};

Tally score(const std::vector<uint8_t>& gt, const std::vector<uint8_t>& got,
            const std::vector<uint8_t>* ref) {
    Tally t;
    t.rays = gt.size();
    for (size_t i = 0; i < gt.size(); ++i) {
        if (ref && got[i] != (*ref)[i]) ++t.other_than_ref;
        if (gt[i] == 2) { ++t.ties; continue; }
        ++t.decided;
        t.blocked += gt[i];
        if (got[i] != gt[i]) ++t.wrong;
    }
    return t;
}

int usage(const char* argv0) {
    std::fprintf(stderr,
                 "usage: %s [--device all|N] [--pause-ms N] [--rays-per-direction N]\n"
                 "  --device              which OpenCL GPU to check (default: all, one at a time)\n"
                 "  --pause-ms            host sleep between dispatches on a device (default 500,\n"
                 "                        the illumination path's validated setting)\n"
                 "  --rays-per-direction  rays per Sun direction per scene and class (default 800)\n"
                 "exit: 0 PASS, 1 FAIL, 2 SKIP (no GPU support or device), 3 INVALID\n",
                 argv0);
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    int device_sel = -1;  // -1 = all
    int pause_ms = 500;
    int per_dir = 800;
    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        const bool has_val = a + 1 < argc;
        if (arg == "--device" && has_val) {
            const std::string v = argv[++a];
            device_sel = (v == "all") ? -1 : std::atoi(v.c_str());
            if (v != "all" && (device_sel < 0 || v.find_first_not_of("0123456789") != std::string::npos))
                return usage(argv[0]);
        } else if (arg == "--pause-ms" && has_val) {
            pause_ms = std::atoi(argv[++a]);
            if (pause_ms < 0) return usage(argv[0]);
        } else if (arg == "--rays-per-direction" && has_val) {
            per_dir = std::atoi(argv[++a]);
            if (per_dir < 50) return usage(argv[0]);
        } else {
            return usage(argv[0]);
        }
    }

#ifndef PTM_GPU_OPENCL
    (void)device_sel; (void)pause_ms; (void)per_dir;
    std::printf("[gpu-selfcheck] SKIP: this build has no OpenCL GPU support "
                "(configure with -DENABLE_GPU=ON). Nothing was checked.\n");
    return 2;
#else
    using namespace gpu::opencl;
    std::printf("[gpu-selfcheck] GPU shadow-ray tracer vs float64 ground truth "
                "(synthetic meshes, no data files)\n");

    std::vector<Scene> scenes;
    scenes.push_back(makeCrater("crater_20m_facets", 500.0, 81));
    scenes.push_back(makeCrater("crater_2mm_facets", 0.05, 81));
    // A large mesh for a deep BVH (~19 levels): a ray's traversal stack then
    // holds more entries than the kernel's former 30-entry cap allowed for at
    // the smaller scales (PTM-795 dropped whole subtrees past that cap).
    scenes.push_back(makeCrater("crater_deep_bvh", 500.0, 361));
    const size_t kDeep = 2;

    struct Case { size_t scene; Rays rays; std::vector<uint8_t> gt; };
    std::vector<Case> cases;
    for (size_t si = 0; si < scenes.size(); ++si) {
        if (si == kDeep) {
            Case c{si, makeRays(scenes[si], 0, std::max(50, per_dir / 2), 20261104u), {}};
            c.gt = groundTruth(scenes[si], c.rays);
            cases.push_back(std::move(c));
            continue;
        }
        for (int cls = 0; cls < 3; ++cls) {
            Case c{si, makeRays(scenes[si], cls, per_dir, 20261103u + 31u * si + cls), {}};
            c.gt = groundTruth(scenes[si], c.rays);
            cases.push_back(std::move(c));
        }
    }

    // The check must be able to fail: every case needs both answers, and few ties.
    bool invalid = false;
    for (const auto& c : cases) {
        Tally t = score(c.gt, c.gt, nullptr);
        const double fb = t.decided ? double(t.blocked) / t.decided : 0.0;
        const double ft = t.rays ? double(t.ties) / t.rays : 1.0;
        if (fb < 0.05 || fb > 0.95 || ft > 0.02) {
            std::printf("[gpu-selfcheck] INVALID: %s %s is not a usable test (%.1f%% blocked, "
                        "%.2f%% near-ties)\n", scenes[c.scene].name.c_str(), c.rays.cls.c_str(),
                        100.0 * fb, 100.0 * ft);
            invalid = true;
        }
    }
    if (invalid) return 3;

    const auto devices = OpenCLContext::enumerateGPUs(true);
    if (devices.empty()) {
        std::printf("[gpu-selfcheck] SKIP: no OpenCL GPU device found. Nothing was checked; "
                    "GPU illumination cannot run here.\n");
        return 2;
    }
    if (device_sel >= static_cast<int>(devices.size())) {
        std::printf("[gpu-selfcheck] SKIP: --device %d, but only %zu device(s) found.\n",
                    device_sel, devices.size());
        return 2;
    }

    BVHConfig bcfg;  // the production settings (illumination_model.cpp)
    bcfg.max_triangles_per_leaf = 4;
    bcfg.use_sah = true;
    bcfg.verbose = false;
    std::vector<BVHBuilder> bvhs(scenes.size());
    std::vector<std::vector<std::array<float, 3>>> verts(scenes.size());
    for (size_t si = 0; si < scenes.size(); ++si) {
        for (const auto& v : scenes[si].V)
            verts[si].push_back({static_cast<float>(v[0]), static_cast<float>(v[1]),
                                 static_cast<float>(v[2])});
        bvhs[si].build(verts[si], scenes[si].F, bcfg);
    }
    {
        const int depth = bvhDepth(bvhs[kDeep]);
        std::printf("[gpu-selfcheck] BVH depth: %d / %d / %d levels\n", bvhDepth(bvhs[0]),
                    bvhDepth(bvhs[1]), depth);
        if (depth < 17) {
            std::printf("[gpu-selfcheck] INVALID: the deep-BVH scene built only %d levels; it "
                        "does not exercise the traversal stack.\n", depth);
            return 3;
        }
    }

    auto trace = [&](BVHRayTracer& rt, const Rays& r, std::vector<uint8_t>& occ) {
        rt.traceShadowRaysEx(r.ox, r.oy, r.oz, r.dx, r.dy, r.dz, r.tmin, r.tmax,
                             r.skip, r.dn, r.heps, occ);
    };

    // Reference: the same kernel text compiled for the host (the engine's
    // PTM_BVH_HOST_KERNEL path), on the first device's context. If it is wrong,
    // the harness is wrong, and no device verdict can be trusted.
    // The variable is unset afterwards whatever it was: inherited from the
    // caller's environment it would turn the device check into a second host
    // check and a PASS that tested nothing.
    std::vector<std::vector<uint8_t>> ref(cases.size());
    {
        if (std::getenv("PTM_BVH_HOST_KERNEL"))
            std::printf("[gpu-selfcheck] NOTE: PTM_BVH_HOST_KERNEL is set in the environment; "
                        "it is ignored here (the device check always runs on the device).\n");
        setenv("PTM_BVH_HOST_KERNEL", "1", 1);
        bool ok = true;
        for (size_t si = 0; si < scenes.size() && ok; ++si) {
            BVHRayTracer rt;
            ok = rt.initialize(devices[0]);
            if (!ok) break;
            rt.uploadMesh(verts[si], scenes[si].F);
            rt.uploadBVH(bvhs[si]);
            for (size_t ci = 0; ci < cases.size(); ++ci)
                if (cases[ci].scene == si) {
                    try { trace(rt, cases[ci].rays, ref[ci]); }
                    catch (const std::exception& e) {
                        std::printf("[gpu-selfcheck] INVALID: host reference failed: %s\n", e.what());
                        ok = false;
                    }
                }
        }
        unsetenv("PTM_BVH_HOST_KERNEL");
        if (!ok) {
            std::printf("[gpu-selfcheck] INVALID: the host reference could not be set up.\n");
            return 3;
        }
        for (size_t ci = 0; ci < cases.size(); ++ci) {
            const Tally t = score(cases[ci].gt, ref[ci], nullptr);
            if (t.wrong) {
                std::printf("[gpu-selfcheck] INVALID: the host-compiled kernel disagrees with the "
                            "float64 ground truth on %zu rays (%s %s); the check itself is "
                            "broken.\n", t.wrong, scenes[cases[ci].scene].name.c_str(),
                            cases[ci].rays.cls.c_str());
                return 3;
            }
        }
    }

    std::printf("[gpu-selfcheck] %zu device(s); pause between dispatches %d ms\n",
                devices.size(), pause_ms);
    int n_checked = 0, n_pass = 0;
    for (size_t di = 0; di < devices.size(); ++di) {
        if (device_sel >= 0 && static_cast<int>(di) != device_sel) continue;
        ++n_checked;
        std::printf("\n[gpu-selfcheck] device %zu: %s\n", di, devices[di].name.c_str());
        bool dev_ok = true, aborted = false;  // wrong answers: keep reporting; errors: stop
        double worst_ms = 0.0;
        bool first_dispatch = true;
        for (size_t si = 0; si < scenes.size() && !aborted; ++si) {
            BVHRayTracer rt;
            if (!rt.initialize(devices[di])) {
                std::printf("  FAIL: could not initialize the ray tracer on this device\n");
                dev_ok = false;
                aborted = true;
                break;
            }
            rt.uploadMesh(verts[si], scenes[si].F);
            rt.uploadBVH(bvhs[si]);
            for (size_t ci = 0; ci < cases.size(); ++ci) {
                if (cases[ci].scene != si) continue;
                if (!first_dispatch && pause_ms > 0)
                    std::this_thread::sleep_for(std::chrono::milliseconds(pause_ms));
                first_dispatch = false;
                std::vector<uint8_t> occ;
                const auto t0 = std::chrono::steady_clock::now();
                try { trace(rt, cases[ci].rays, occ); }
                catch (const std::exception& e) {
                    std::printf("  FAIL: %s %s: %s\n", scenes[si].name.c_str(),
                                cases[ci].rays.cls.c_str(), e.what());
                    dev_ok = false;
                    aborted = true;
                    break;
                }
                const double ms = std::chrono::duration<double, std::milli>(
                                      std::chrono::steady_clock::now() - t0).count();
                worst_ms = std::max(worst_ms, ms);
                const Tally t = score(cases[ci].gt, occ, &ref[ci]);
                std::printf("  %-18s %-8s %6zu rays: %6zu decided (%5.1f%% blocked), %3zu near-ties "
                            "skipped | WRONG %zu | differs from host kernel on %zu | %.1f ms\n",
                            scenes[si].name.c_str(), cases[ci].rays.cls.c_str(), t.rays, t.decided,
                            100.0 * t.blocked / std::max<size_t>(1, t.decided), t.ties, t.wrong,
                            t.other_than_ref, ms);
                if (t.wrong) dev_ok = false;
            }
        }
        std::printf("[gpu-selfcheck] device %zu: %s (longest dispatch %.1f ms)\n", di,
                    dev_ok ? "PASS" : "FAIL", worst_ms);
        n_pass += dev_ok;
    }

    const bool all_ok = (n_pass == n_checked);
    std::printf("\n[gpu-selfcheck] RESULT: %s (%d of %d device(s) passed)\n",
                all_ok ? "PASS" : "FAIL", n_pass, n_checked);
    if (!all_ok)
        std::printf("[gpu-selfcheck] Do not use GPU illumination on the failing device(s); "
                    "keep illumination.use_gpu: false (CPU Embree).\n");
    return all_ok ? 0 : 1;
#endif
}
