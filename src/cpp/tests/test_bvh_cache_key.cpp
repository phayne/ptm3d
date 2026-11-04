// test_bvh_cache_key -- the GPU BVH disk cache must be keyed on content (PTM-800).
//
// The cache (illumination.bvh_cache_file, ON by default at
// <output>/bvh_cache.bvh1) used to be accepted when the triangle COUNT
// matched. A different mesh with the same face count then reused a BVH
// bounding other triangles. This test drives BVHBuilder::contentKey / write /
// read, the calls IlluminationModel makes, on synthetic meshes. No device.
//
// C1 control: the perturbed mesh has the same face count, so the old rule
//    would have accepted the stale cache.
// T1 same content -> hit, identical nodes.   T2 one vertex moved 1 mm -> miss.
// T3 one triangle's indices rotated -> miss. T4 scene offset changed -> miss.
// T5 leaf size changed -> miss.              T6 an old BVH1 file -> miss.
// T7 a truncated file -> miss.               T8 a miss leaves the builder empty.

#include "bvh_opencl.h"

#include <array>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace gpu::opencl;
using V3 = std::array<float, 3>;
using T3 = std::array<uint32_t, 3>;

namespace {
int fails = 0;
void check(const char* name, bool ok, const std::string& detail = "") {
    std::printf("  [%s] %s%s%s\n", name, ok ? "PASS" : "FAIL", ok ? "" : "  ", ok ? "" : detail.c_str());
    fails += !ok;
}
void grid(int n, std::vector<V3>& v, std::vector<T3>& f) {
    v.clear(); f.clear();
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            v.push_back({float(i), float(j), 0.1f * float((i * 7 + j * 3) % 5)});
    for (int j = 0; j + 1 < n; ++j)
        for (int i = 0; i + 1 < n; ++i) {
            uint32_t a = j * n + i, b = a + 1, c = a + n + 1, d = a + n;
            f.push_back({a, b, c}); f.push_back({a, c, d});
        }
}
bool sameNodes(const BVHBuilder& x, const BVHBuilder& y) {
    if (x.numNodes() != y.numNodes() || x.triIndices() != y.triIndices()) return false;
    return std::memcmp(x.nodes().data(), y.nodes().data(), x.numNodes() * sizeof(BVHNode)) == 0;
}
}  // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : ".";
    const std::string path = dir + "/test_bvh_cache_key.bvh";
    std::vector<V3> v; std::vector<T3> f; grid(40, v, f);
    const std::array<double, 3> off = {1000.0, -2000.0, 1737400.0};
    BVHConfig cfg; cfg.max_triangles_per_leaf = 4; cfg.use_sah = true;

    BVHBuilder built; built.build(v, f, cfg);
    const uint64_t key = BVHBuilder::contentKey(v, f, off, cfg);
    if (!built.write(path, key)) { std::printf("cannot write %s\n", path.c_str()); return 2; }
    std::string why;

    BVHBuilder hit;
    check("T1 same content hits", hit.read(path, key, f.size(), &why) && sameNodes(hit, built), why);

    std::vector<V3> moved = v; moved[25 * 40 + 17][2] += 0.001f;  // same face count
    check("C1 control: perturbed mesh has the same face count (old rule accepts)",
          moved.size() == v.size() && f.size() == built.numTriangles());
    const uint64_t key_moved = BVHBuilder::contentKey(moved, f, off, cfg);
    BVHBuilder m;
    const bool m_hit = m.read(path, key_moved, f.size(), &why);
    check("T2 one vertex moved 1 mm misses", !m_hit && why.find("content key") != std::string::npos, why);
    check("T8 a miss leaves the builder empty", m.numNodes() == 0 && m.numTriangles() == 0);

    std::vector<T3> rot = f; rot[100] = {rot[100][1], rot[100][2], rot[100][0]};
    BVHBuilder r;
    check("T3 rotated indices miss", !r.read(path, BVHBuilder::contentKey(v, rot, off, cfg), f.size(), &why));

    std::array<double, 3> off2 = off; off2[0] += 1.0;
    BVHBuilder o;
    check("T4 scene offset change misses", !o.read(path, BVHBuilder::contentKey(v, f, off2, cfg), f.size(), &why));

    BVHConfig cfg2 = cfg; cfg2.max_triangles_per_leaf = 8;
    BVHBuilder c;
    check("T5 leaf size change misses", !c.read(path, BVHBuilder::contentKey(v, f, off, cfg2), f.size(), &why));

    {   // the pre-PTM-800 format: magic "BVH1", no key, same counts
        const std::string p1 = dir + "/test_bvh_cache_key_v1.bvh";
        std::ofstream o1(p1, std::ios::binary);
        const uint32_t magic = 0x31485642, nn = uint32_t(built.numNodes()),
                       nt = uint32_t(built.numTriangles()), ni = uint32_t(built.triIndices().size());
        o1.write((const char*)&magic, 4); o1.write((const char*)&nn, 4);
        o1.write((const char*)&nt, 4); o1.write((const char*)&ni, 4);
        o1.write((const char*)built.nodes().data(), nn * sizeof(BVHNode));
        o1.write((const char*)built.triIndices().data(), ni * 4);
        o1.close();
        BVHBuilder l;
        check("T6 old BVH1 cache misses", !l.read(p1, key, f.size(), &why) &&
              why.find("BVH1") != std::string::npos, why);
        std::remove(p1.c_str());
    }
    {
        const std::string pt = dir + "/test_bvh_cache_key_trunc.bvh";
        std::ifstream in(path, std::ios::binary);
        std::vector<char> buf((std::istreambuf_iterator<char>(in)), {});
        std::ofstream(pt, std::ios::binary).write(buf.data(), std::streamsize(buf.size() / 2));
        BVHBuilder t;
        check("T7 truncated cache misses", !t.read(pt, key, f.size(), &why) && t.numNodes() == 0, why);
        std::remove(pt.c_str());
    }
    std::remove(path.c_str());
    std::printf("[test_bvh_cache_key] %s\n", fails ? "FAIL" : "PASS");
    return fails ? 1 : 0;
}
