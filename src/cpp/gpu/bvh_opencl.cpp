/**
 * @file bvh_opencl.cpp
 * @brief Implementation of OpenCL-accelerated BVH ray tracer
 */

#ifdef PTM_GPU_OPENCL

#include "bvh_opencl.h"
#include "gpu_dispatch_log.h"
#include <fstream>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <iostream>
#include <thread>
#include <cmath>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace gpu {
namespace opencl {

//=============================================================================
// Embedded OpenCL kernel source
//=============================================================================

#define PTM_BVH_KERNEL(...) #__VA_ARGS__
static const char* kBVHKernelText =
#include "bvh_shadow_kernel.inc"
;
#undef PTM_BVH_KERNEL

std::string getBVHKernelSource() {
    // No FP contraction: the watertight test needs the shared-edge products
    // rounded separately so neighbors compute exactly negated edge functions.
    return std::string("#pragma OPENCL FP_CONTRACT OFF\n") + kBVHKernelText;
}

// Host execution of the SAME kernel text (test-only reference path).
namespace host_kernel {
typedef unsigned int uint;
typedef unsigned char uchar;
static thread_local int g_global_id = 0;
static inline int get_global_id(int) { return g_global_id; }
static inline float max(float a, float b) { return std::fmax(a, b); }
static inline float min(float a, float b) { return std::fmin(a, b); }
#define __kernel static
#define __global
#define PTM_BVH_KERNEL(...) __VA_ARGS__
#include "bvh_shadow_kernel.inc"
#undef PTM_BVH_KERNEL
#undef __kernel
#undef __global
}  // namespace host_kernel

//=============================================================================
// BVHBuilder implementation
//=============================================================================

void BVHBuilder::build(const std::vector<std::array<float, 3>>& vertices,
                       const std::vector<std::array<uint32_t, 3>>& indices,
                       const BVHConfig& config) {
    nodes_.clear();
    tri_indices_.clear();
    num_triangles_ = indices.size();

    if (num_triangles_ == 0) {
        return;
    }

    // Build primitives with centroids and bounds
    std::vector<BuildPrimitive> primitives(num_triangles_);

    for (size_t i = 0; i < num_triangles_; ++i) {
        const auto& tri = indices[i];
        const auto& v0 = vertices[tri[0]];
        const auto& v1 = vertices[tri[1]];
        const auto& v2 = vertices[tri[2]];

        primitives[i].index = static_cast<uint32_t>(i);

        // Compute bounding box
        primitives[i].min_bound = {{
            std::min({v0[0], v1[0], v2[0]}),
            std::min({v0[1], v1[1], v2[1]}),
            std::min({v0[2], v1[2], v2[2]})
        }};
        primitives[i].max_bound = {{
            std::max({v0[0], v1[0], v2[0]}),
            std::max({v0[1], v1[1], v2[1]}),
            std::max({v0[2], v1[2], v2[2]})
        }};

        // Compute centroid
        primitives[i].centroid = {{
            (v0[0] + v1[0] + v2[0]) / 3.0f,
            (v0[1] + v1[1] + v2[1]) / 3.0f,
            (v0[2] + v1[2] + v2[2]) / 3.0f
        }};
    }

    // Reserve space for nodes (roughly 2N-1 for N triangles)
    nodes_.reserve(2 * num_triangles_);
    tri_indices_.reserve(num_triangles_);

    // Build recursively
    buildRecursive(primitives, 0, num_triangles_, 0, config);

    if (config.verbose) {
        std::cout << "BVH built: " << nodes_.size() << " nodes, "
                  << num_triangles_ << " triangles" << std::endl;
    }
}

void BVHBuilder::build(const std::vector<std::array<double, 3>>& vertices,
                       const std::vector<std::array<int, 3>>& indices,
                       const BVHConfig& config) {
    // Convert double to float
    std::vector<std::array<float, 3>> float_vertices(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        float_vertices[i] = {{
            static_cast<float>(vertices[i][0]),
            static_cast<float>(vertices[i][1]),
            static_cast<float>(vertices[i][2])
        }};
    }

    // Convert int to uint32_t
    std::vector<std::array<uint32_t, 3>> uint_indices(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        uint_indices[i] = {{
            static_cast<uint32_t>(indices[i][0]),
            static_cast<uint32_t>(indices[i][1]),
            static_cast<uint32_t>(indices[i][2])
        }};
    }

    build(float_vertices, uint_indices, config);
}

uint32_t BVHBuilder::buildRecursive(std::vector<BuildPrimitive>& primitives,
                                     size_t start, size_t end,
                                     int depth, const BVHConfig& config) {
    uint32_t node_idx = static_cast<uint32_t>(nodes_.size());
    nodes_.emplace_back();
    BVHNode& node = nodes_.back();

    size_t count = end - start;

    // Compute bounds of all primitives in this node
    node.min_bound = {{1e30f, 1e30f, 1e30f}};
    node.max_bound = {{-1e30f, -1e30f, -1e30f}};

    for (size_t i = start; i < end; ++i) {
        for (int j = 0; j < 3; ++j) {
            node.min_bound[j] = std::min(node.min_bound[j], primitives[i].min_bound[j]);
            node.max_bound[j] = std::max(node.max_bound[j], primitives[i].max_bound[j]);
        }
    }

    // Make leaf if few primitives or max depth
    if (count <= static_cast<size_t>(config.max_triangles_per_leaf) ||
        depth >= config.max_depth) {
        // Leaf node
        node.tri_start = static_cast<uint32_t>(tri_indices_.size());
        node.tri_count = static_cast<uint32_t>(count);
        node.left_child = 0;
        node.right_child = 0;

        for (size_t i = start; i < end; ++i) {
            tri_indices_.push_back(primitives[i].index);
        }
        return node_idx;
    }

    // Find best split using SAH
    int best_axis = 0;
    float best_pos = 0;
    float best_cost = 1e30f;

    if (config.use_sah) {
        // Compute centroid bounds
        std::array<float, 3> centroid_min = {{1e30f, 1e30f, 1e30f}};
        std::array<float, 3> centroid_max = {{-1e30f, -1e30f, -1e30f}};

        for (size_t i = start; i < end; ++i) {
            for (int j = 0; j < 3; ++j) {
                centroid_min[j] = std::min(centroid_min[j], primitives[i].centroid[j]);
                centroid_max[j] = std::max(centroid_max[j], primitives[i].centroid[j]);
            }
        }

        // Try each axis
        for (int axis = 0; axis < 3; ++axis) {
            float extent = centroid_max[axis] - centroid_min[axis];
            if (extent < 1e-6f) continue;

            // Initialize bins
            std::vector<SAHBucket> buckets(config.sah_bins);

            // Fill bins
            for (size_t i = start; i < end; ++i) {
                float rel_pos = (primitives[i].centroid[axis] - centroid_min[axis]) / extent;
                int bucket_idx = std::min(static_cast<int>(rel_pos * config.sah_bins),
                                         config.sah_bins - 1);
                buckets[bucket_idx].count++;
                for (int j = 0; j < 3; ++j) {
                    buckets[bucket_idx].min_bound[j] = std::min(buckets[bucket_idx].min_bound[j],
                                                                primitives[i].min_bound[j]);
                    buckets[bucket_idx].max_bound[j] = std::max(buckets[bucket_idx].max_bound[j],
                                                                primitives[i].max_bound[j]);
                }
            }

            // Evaluate each split position
            for (int split = 1; split < config.sah_bins; ++split) {
                // Left side
                int left_count = 0;
                std::array<float, 3> left_min = {{1e30f, 1e30f, 1e30f}};
                std::array<float, 3> left_max = {{-1e30f, -1e30f, -1e30f}};

                for (int b = 0; b < split; ++b) {
                    left_count += buckets[b].count;
                    for (int j = 0; j < 3; ++j) {
                        left_min[j] = std::min(left_min[j], buckets[b].min_bound[j]);
                        left_max[j] = std::max(left_max[j], buckets[b].max_bound[j]);
                    }
                }

                // Right side
                int right_count = 0;
                std::array<float, 3> right_min = {{1e30f, 1e30f, 1e30f}};
                std::array<float, 3> right_max = {{-1e30f, -1e30f, -1e30f}};

                for (int b = split; b < config.sah_bins; ++b) {
                    right_count += buckets[b].count;
                    for (int j = 0; j < 3; ++j) {
                        right_min[j] = std::min(right_min[j], buckets[b].min_bound[j]);
                        right_max[j] = std::max(right_max[j], buckets[b].max_bound[j]);
                    }
                }

                if (left_count == 0 || right_count == 0) continue;

                // Compute SAH cost
                float parent_area = surfaceArea(node.min_bound, node.max_bound);
                float left_area = surfaceArea(left_min, left_max);
                float right_area = surfaceArea(right_min, right_max);

                float cost = config.traversal_cost +
                            config.intersection_cost * (
                                left_count * left_area / parent_area +
                                right_count * right_area / parent_area
                            );

                if (cost < best_cost) {
                    best_cost = cost;
                    best_axis = axis;
                    best_pos = centroid_min[axis] + extent * split / config.sah_bins;
                }
            }
        }
    }

    // Fall back to median split if SAH didn't find good split
    if (best_cost >= static_cast<float>(count) * config.intersection_cost) {
        // Find longest axis
        float max_extent = 0;
        for (int axis = 0; axis < 3; ++axis) {
            float extent = node.max_bound[axis] - node.min_bound[axis];
            if (extent > max_extent) {
                max_extent = extent;
                best_axis = axis;
            }
        }

        // Use median position
        std::nth_element(primitives.begin() + start,
                        primitives.begin() + start + count / 2,
                        primitives.begin() + end,
                        [best_axis](const BuildPrimitive& a, const BuildPrimitive& b) {
                            return a.centroid[best_axis] < b.centroid[best_axis];
                        });
        best_pos = primitives[start + count / 2].centroid[best_axis];
    }

    // Partition primitives
    auto mid = std::partition(primitives.begin() + start, primitives.begin() + end,
                              [best_axis, best_pos](const BuildPrimitive& p) {
                                  return p.centroid[best_axis] < best_pos;
                              });
    size_t mid_idx = mid - primitives.begin();

    // Handle degenerate cases
    if (mid_idx == start || mid_idx == end) {
        mid_idx = start + count / 2;
    }

    // Mark as internal node
    node.tri_count = 0;
    node.tri_start = 0;

    // Build children (update node reference after recursive calls since vector may reallocate)
    uint32_t left_idx = buildRecursive(primitives, start, mid_idx, depth + 1, config);
    uint32_t right_idx = buildRecursive(primitives, mid_idx, end, depth + 1, config);

    nodes_[node_idx].left_child = left_idx;
    nodes_[node_idx].right_child = right_idx;

    return node_idx;
}

float BVHBuilder::surfaceArea(const std::array<float, 3>& min_b,
                              const std::array<float, 3>& max_b) {
    float dx = max_b[0] - min_b[0];
    float dy = max_b[1] - min_b[1];
    float dz = max_b[2] - min_b[2];
    return 2.0f * (dx * dy + dy * dz + dz * dx);
}

namespace {
// FNV-1a 64, byte-wise. Not cryptographic: the cache is the user's own file
// and the threat is an accidental reuse, which any 64-bit hash of the full
// content rules out.
struct Fnv64 {
    uint64_t h = 1469598103934665603ULL;
    void bytes(const void* p, size_t n) {
        const unsigned char* c = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) { h ^= c[i]; h *= 1099511628211ULL; }
    }
    template <class T> void pod(const T& v) { bytes(&v, sizeof(v)); }
};
}  // namespace

uint64_t BVHBuilder::contentKey(const std::vector<std::array<float, 3>>& vertices,
                                const std::vector<std::array<uint32_t, 3>>& indices,
                                const std::array<double, 3>& scene_offset,
                                const BVHConfig& config) {
    Fnv64 f;
    f.pod(BVH_MAGIC);
    f.pod(static_cast<uint32_t>(sizeof(BVHNode)));
    f.pod(static_cast<uint64_t>(vertices.size()));
    f.pod(static_cast<uint64_t>(indices.size()));
    if (!vertices.empty()) f.bytes(vertices.data(), vertices.size() * sizeof(vertices[0]));
    if (!indices.empty()) f.bytes(indices.data(), indices.size() * sizeof(indices[0]));
    f.pod(scene_offset);
    f.pod(config.max_triangles_per_leaf);
    f.pod(config.max_depth);
    f.pod(config.use_sah);
    f.pod(config.sah_bins);
    f.pod(config.traversal_cost);
    f.pod(config.intersection_cost);
    return f.h;
}

bool BVHBuilder::write(const std::string& filename, uint64_t content_key) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    file.write(reinterpret_cast<const char*>(&BVH_MAGIC), sizeof(BVH_MAGIC));
    file.write(reinterpret_cast<const char*>(&content_key), sizeof(content_key));

    uint32_t num_nodes = static_cast<uint32_t>(nodes_.size());
    uint32_t num_tris = static_cast<uint32_t>(num_triangles_);
    uint32_t num_indices = static_cast<uint32_t>(tri_indices_.size());

    file.write(reinterpret_cast<const char*>(&num_nodes), sizeof(num_nodes));
    file.write(reinterpret_cast<const char*>(&num_tris), sizeof(num_tris));
    file.write(reinterpret_cast<const char*>(&num_indices), sizeof(num_indices));

    file.write(reinterpret_cast<const char*>(nodes_.data()),
               nodes_.size() * sizeof(BVHNode));
    file.write(reinterpret_cast<const char*>(tri_indices_.data()),
               tri_indices_.size() * sizeof(uint32_t));

    return file.good();
}

bool BVHBuilder::read(const std::string& filename, uint64_t expected_key,
                      size_t expected_triangles, std::string* why) {
    auto reject = [&](const std::string& reason) {
        nodes_.clear();
        tri_indices_.clear();
        num_triangles_ = 0;
        if (why) *why = reason;
        return false;
    };
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return reject("no cache file");
    }

    uint32_t magic = 0;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    if (!file.good()) return reject("truncated header");
    if (magic == BVH_MAGIC_V1) {
        return reject("old unkeyed BVH1 cache (accepted on triangle count alone; PTM-800)");
    }
    if (magic != BVH_MAGIC) return reject("not a BVH cache file");

    uint64_t key = 0;
    file.read(reinterpret_cast<char*>(&key), sizeof(key));
    if (!file.good()) return reject("truncated header");
    if (key != expected_key) {
        return reject("built from different geometry or build settings (content key mismatch)");
    }

    uint32_t num_nodes = 0, num_tris = 0, num_indices = 0;
    file.read(reinterpret_cast<char*>(&num_nodes), sizeof(num_nodes));
    file.read(reinterpret_cast<char*>(&num_tris), sizeof(num_tris));
    file.read(reinterpret_cast<char*>(&num_indices), sizeof(num_indices));
    if (!file.good()) return reject("truncated header");
    if (num_tris != expected_triangles) return reject("triangle count mismatch");

    nodes_.resize(num_nodes);
    file.read(reinterpret_cast<char*>(nodes_.data()), num_nodes * sizeof(BVHNode));
    tri_indices_.resize(num_indices);
    file.read(reinterpret_cast<char*>(tri_indices_.data()), num_indices * sizeof(uint32_t));
    if (!file.good()) return reject("truncated body");

    num_triangles_ = num_tris;
    if (why) why->clear();
    return true;
}

// ---------------------------------------------------------------------------
// TEST-ONLY ray dump (PTM-794 forensics). PTM_BVH_DUMP_DIR=<dir> with
// PTM_BVH_DUMP_PRIMS=<i,j,...> writes, from every tracer in the process:
//   <dir>/mesh_f32.bin   uint64 nv, uint64 nf, float32 xyz[nv], uint32 idx[nf][3]
//                        (the exact geometry uploaded, in the tracer's frame)
//   <dir>/rays.bin       per ray whose skip_prim is listed: uint64 call,
//                        uint64 ray, float32 org[3], dir[3], tmin, tmax, int32 skip,
//                        float32 dn, heps, uint8 occluded, uint8 host_kernel
// Never on by default; announces itself once.
namespace {
struct RayDump {
    std::string dir;
    std::vector<int32_t> prims;
    bool on = false;
    uint64_t calls = 0;
    std::mutex mu;
    RayDump() {
        const char* d = std::getenv("PTM_BVH_DUMP_DIR");
        const char* p = std::getenv("PTM_BVH_DUMP_PRIMS");
        if (!d || !*d) return;
        dir = d;
        if (p) {
            std::string tok; std::string all = p;
            for (size_t i = 0; i <= all.size(); ++i) {
                if (i == all.size() || all[i] == ',') {
                    if (!tok.empty()) prims.push_back(static_cast<int32_t>(std::atol(tok.c_str())));
                    tok.clear();
                } else tok += all[i];
            }
        }
        std::sort(prims.begin(), prims.end());
        on = true;
        std::cerr << "[BVHRayTracer] WARNING: PTM_BVH_DUMP_DIR=" << dir << " -- dumping the "
                     "uploaded mesh and the rays of " << prims.size() << " source facet(s). "
                     "TEST ONLY." << std::endl;
    }
    bool wants(int32_t prim) const { return std::binary_search(prims.begin(), prims.end(), prim); }
};
RayDump& rayDump() { static RayDump d; return d; }
}  // namespace

//=============================================================================
// BVHRayTracer implementation
//=============================================================================

BVHRayTracer::BVHRayTracer() = default;

BVHRayTracer::~BVHRayTracer() {
    release();
}

bool BVHRayTracer::initialize() {
    if (initialized_) return true;

    try {
        context_ = std::make_unique<OpenCLContext>();
        if (!context_->initialize(true)) {
            std::cerr << "Warning: Failed to initialize OpenCL for BVH ray tracer" << std::endl;
            return false;
        }

        buildKernels();

        // Get optimal work group size
        work_group_size_ = kernelWorkGroupSize();
        if (const char* hk = std::getenv("PTM_BVH_HOST_KERNEL")) {
            if (hk[0] == '1') {
                host_kernel_ = true;
                std::cerr << "[BVHRayTracer] WARNING: PTM_BVH_HOST_KERNEL=1 -- shadow rays are "
                             "traced by the kernel text compiled for the HOST, not by the device. "
                             "TEST ONLY (reference for the device result)." << std::endl;
            }
        }

        initialized_ = true;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error initializing BVH ray tracer: " << e.what() << std::endl;
        return false;
    }
}

bool BVHRayTracer::initialize(const OpenCLDeviceInfo& device_handle) {
    if (initialized_) return true;

    try {
        context_ = std::make_unique<OpenCLContext>();
        if (!context_->initialize(device_handle)) {
            std::cerr << "Warning: Failed to initialize OpenCL for BVH ray tracer on device "
                      << device_handle.name << std::endl;
            return false;
        }

        buildKernels();

        work_group_size_ = kernelWorkGroupSize();
        if (const char* hk = std::getenv("PTM_BVH_HOST_KERNEL")) {
            if (hk[0] == '1') {
                host_kernel_ = true;
                std::cerr << "[BVHRayTracer] WARNING: PTM_BVH_HOST_KERNEL=1 -- shadow rays are "
                             "traced by the kernel text compiled for the HOST, not by the device. "
                             "TEST ONLY (reference for the device result)." << std::endl;
            }
        }

        initialized_ = true;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error initializing BVH ray tracer on " << device_handle.name
                  << ": " << e.what() << std::endl;
        return false;
    }
}

void BVHRayTracer::buildKernels() {
    std::string source = getBVHKernelSource();

    // Build program manually (not using context_->buildProgram to avoid exception issues
    // with CL_CHECK macro across shared library boundaries on macOS)
    const char* src = source.c_str();
    size_t len = source.size();

    cl_int err;
    program_ = clCreateProgramWithSource(context_->context(), 1, &src, &len, &err);
    if (err != CL_SUCCESS) {
        throw std::runtime_error(std::string("clCreateProgramWithSource failed: ") +
                                 clErrorString(err));
    }

    // Use simpler build options to avoid driver issues on some platforms
    // NO -cl-fast-relaxed-math. It was here (PTM-794) and licensed approximate
    // reciprocals and finite-math assumptions in the slab and Moller-Trumbore
    // tests, which is incompatible with matching the CPU (Embree) shadow
    // answer. PTM_BVH_CL_OPTIONS overrides, for experiments only.
    std::string options = "";
    if (const char* e = std::getenv("PTM_BVH_CL_OPTIONS")) {
        options = e;
        std::cerr << "[BVHRayTracer] WARNING: kernel build options overridden by "
                     "PTM_BVH_CL_OPTIONS='" << options << "' (experiment only)" << std::endl;
    }
    cl_device_id device = context_->device();
    err = clBuildProgram(program_, 1, &device, options.c_str(), nullptr, nullptr);
    if (err != CL_SUCCESS) {
        // Get build log
        size_t log_size;
        clGetProgramBuildInfo(program_, device, CL_PROGRAM_BUILD_LOG,
                             0, nullptr, &log_size);
        std::vector<char> log(log_size);
        clGetProgramBuildInfo(program_, device, CL_PROGRAM_BUILD_LOG,
                             log_size, log.data(), nullptr);

        std::string msg = std::string("clBuildProgram failed: ") + clErrorString(err) +
                         "\nBuild log:\n" + log.data();
        clReleaseProgram(program_);
        program_ = nullptr;
        throw std::runtime_error(msg);
    }

    // Create kernel
    shadow_kernel_ = clCreateKernel(program_, "shadow_ray_kernel", &err);
    if (err != CL_SUCCESS) {
        throw std::runtime_error(std::string("Failed to create shadow_ray_kernel: ") +
                                 clErrorString(err));
    }
}

void BVHRayTracer::uploadMesh(const std::vector<std::array<float, 3>>& vertices,
                              const std::vector<std::array<uint32_t, 3>>& indices) {
    if (!initialized_) {
        std::cerr << "BVHRayTracer::uploadMesh: Not initialized" << std::endl;
        return;
    }

    num_triangles_ = indices.size();

    if (rayDump().on) {
        std::lock_guard<std::mutex> lk(rayDump().mu);
        std::ofstream mf(rayDump().dir + "/mesh_f32.bin", std::ios::binary);
        const uint64_t nv = vertices.size(), nf = indices.size();
        mf.write(reinterpret_cast<const char*>(&nv), 8);
        mf.write(reinterpret_cast<const char*>(&nf), 8);
        mf.write(reinterpret_cast<const char*>(vertices.data()), nv * sizeof(vertices[0]));
        mf.write(reinterpret_cast<const char*>(indices.data()), nf * sizeof(indices[0]));
    }

    std::cerr << "[BVHRayTracer] Uploading mesh: " << num_triangles_ << " triangles, "
              << vertices.size() << " vertices" << std::endl;

    // Convert to SoA format
    std::vector<float> v0_x(num_triangles_), v0_y(num_triangles_), v0_z(num_triangles_);
    std::vector<float> v1_x(num_triangles_), v1_y(num_triangles_), v1_z(num_triangles_);
    std::vector<float> v2_x(num_triangles_), v2_y(num_triangles_), v2_z(num_triangles_);

    for (size_t i = 0; i < num_triangles_; ++i) {
        const auto& tri = indices[i];
        v0_x[i] = vertices[tri[0]][0];
        v0_y[i] = vertices[tri[0]][1];
        v0_z[i] = vertices[tri[0]][2];
        v1_x[i] = vertices[tri[1]][0];
        v1_y[i] = vertices[tri[1]][1];
        v1_z[i] = vertices[tri[1]][2];
        v2_x[i] = vertices[tri[2]][0];
        v2_y[i] = vertices[tri[2]][1];
        v2_z[i] = vertices[tri[2]][2];
    }

    std::cerr << "[BVHRayTracer] SoA conversion complete" << std::endl;
    if (host_kernel_) {
        h_tri_[0] = v0_x; h_tri_[1] = v0_y; h_tri_[2] = v0_z;
        h_tri_[3] = v1_x; h_tri_[4] = v1_y; h_tri_[5] = v1_z;
        h_tri_[6] = v2_x; h_tri_[7] = v2_y; h_tri_[8] = v2_z;
    }

    // Create and upload buffers
    cl_context ctx = context_->context();
    cl_command_queue queue = context_->queue();
    cl_int err;

    size_t buf_size = num_triangles_ * sizeof(float);

    std::cerr << "[BVHRayTracer] Buffer size: " << (buf_size / (1024*1024.0)) << " MB per array" << std::endl;

    // Release existing buffers
    if (buf_tri_v0_x_) clReleaseMemObject(buf_tri_v0_x_);
    if (buf_tri_v0_y_) clReleaseMemObject(buf_tri_v0_y_);
    if (buf_tri_v0_z_) clReleaseMemObject(buf_tri_v0_z_);
    if (buf_tri_v1_x_) clReleaseMemObject(buf_tri_v1_x_);
    if (buf_tri_v1_y_) clReleaseMemObject(buf_tri_v1_y_);
    if (buf_tri_v1_z_) clReleaseMemObject(buf_tri_v1_z_);
    if (buf_tri_v2_x_) clReleaseMemObject(buf_tri_v2_x_);
    if (buf_tri_v2_y_) clReleaseMemObject(buf_tri_v2_y_);
    if (buf_tri_v2_z_) clReleaseMemObject(buf_tri_v2_z_);

    // Use non-throwing error checks (macOS shared library exception issue)
    #define CL_CHECK_RETURN(call, msg) \
        do { \
            cl_int _err = (call); \
            if (_err != CL_SUCCESS) { \
                std::cerr << "OpenCL error in " << msg << ": " << clErrorString(_err) << std::endl; \
                return; \
            } \
        } while (0)

    std::cerr << "[BVHRayTracer] Creating triangle buffers..." << std::endl;

    buf_tri_v0_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v0_x");
    buf_tri_v0_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v0_y");
    buf_tri_v0_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v0_z");
    buf_tri_v1_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v1_x");
    buf_tri_v1_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v1_y");
    buf_tri_v1_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v1_z");
    buf_tri_v2_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v2_x");
    buf_tri_v2_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v2_y");
    buf_tri_v2_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer v2_z");

    std::cerr << "[BVHRayTracer] Uploading triangle data..." << std::endl;

    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v0_x_, CL_FALSE, 0, buf_size, v0_x.data(), 0, nullptr, nullptr), "write v0_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v0_y_, CL_FALSE, 0, buf_size, v0_y.data(), 0, nullptr, nullptr), "write v0_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v0_z_, CL_FALSE, 0, buf_size, v0_z.data(), 0, nullptr, nullptr), "write v0_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v1_x_, CL_FALSE, 0, buf_size, v1_x.data(), 0, nullptr, nullptr), "write v1_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v1_y_, CL_FALSE, 0, buf_size, v1_y.data(), 0, nullptr, nullptr), "write v1_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v1_z_, CL_FALSE, 0, buf_size, v1_z.data(), 0, nullptr, nullptr), "write v1_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v2_x_, CL_FALSE, 0, buf_size, v2_x.data(), 0, nullptr, nullptr), "write v2_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v2_y_, CL_FALSE, 0, buf_size, v2_y.data(), 0, nullptr, nullptr), "write v2_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_v2_z_, CL_FALSE, 0, buf_size, v2_z.data(), 0, nullptr, nullptr), "write v2_z");

    std::cerr << "[BVHRayTracer] Finishing queue..." << std::endl;

    err = clFinish(queue);
    CL_CHECK_RETURN(err, "clFinish");

    mesh_uploaded_ = true;
    std::cerr << "[BVHRayTracer] Mesh upload complete" << std::endl;

    #undef CL_CHECK_RETURN
}

void BVHRayTracer::uploadMesh(const std::vector<std::array<double, 3>>& vertices,
                              const std::vector<std::array<int, 3>>& indices) {
    // Convert to float
    std::vector<std::array<float, 3>> float_vertices(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        float_vertices[i] = {{
            static_cast<float>(vertices[i][0]),
            static_cast<float>(vertices[i][1]),
            static_cast<float>(vertices[i][2])
        }};
    }

    // Convert to uint32_t
    std::vector<std::array<uint32_t, 3>> uint_indices(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        uint_indices[i] = {{
            static_cast<uint32_t>(indices[i][0]),
            static_cast<uint32_t>(indices[i][1]),
            static_cast<uint32_t>(indices[i][2])
        }};
    }

    uploadMesh(float_vertices, uint_indices);
}

void BVHRayTracer::uploadBVH(const BVHBuilder& bvh) {
    if (!initialized_) {
        throw std::runtime_error("BVHRayTracer not initialized");
    }

    const auto& nodes = bvh.nodes();
    const auto& tri_indices = bvh.triIndices();
    num_nodes_ = nodes.size();

    if (num_nodes_ == 0) {
        throw std::runtime_error("Cannot upload empty BVH");
    }

    // Convert to SoA format
    std::vector<float> min_x(num_nodes_), min_y(num_nodes_), min_z(num_nodes_);
    std::vector<float> max_x(num_nodes_), max_y(num_nodes_), max_z(num_nodes_);
    std::vector<uint32_t> left(num_nodes_), right(num_nodes_);
    std::vector<uint32_t> tri_start(num_nodes_), tri_count(num_nodes_);

    for (size_t i = 0; i < num_nodes_; ++i) {
        const auto& node = nodes[i];
        min_x[i] = node.min_bound[0];
        min_y[i] = node.min_bound[1];
        min_z[i] = node.min_bound[2];
        max_x[i] = node.max_bound[0];
        max_y[i] = node.max_bound[1];
        max_z[i] = node.max_bound[2];
        left[i] = node.left_child;
        right[i] = node.right_child;
        tri_start[i] = node.tri_start;
        tri_count[i] = node.tri_count;
    }

    if (host_kernel_) {
        h_node_[0] = min_x; h_node_[1] = min_y; h_node_[2] = min_z;
        h_node_[3] = max_x; h_node_[4] = max_y; h_node_[5] = max_z;
        h_left_ = left; h_right_ = right; h_tstart_ = tri_start; h_tcount_ = tri_count;
        h_tidx_.assign(tri_indices.begin(), tri_indices.end());
    }

    cl_context ctx = context_->context();
    cl_command_queue queue = context_->queue();
    cl_int err;

    // Release existing BVH buffers
    if (buf_node_min_x_) clReleaseMemObject(buf_node_min_x_);
    if (buf_node_min_y_) clReleaseMemObject(buf_node_min_y_);
    if (buf_node_min_z_) clReleaseMemObject(buf_node_min_z_);
    if (buf_node_max_x_) clReleaseMemObject(buf_node_max_x_);
    if (buf_node_max_y_) clReleaseMemObject(buf_node_max_y_);
    if (buf_node_max_z_) clReleaseMemObject(buf_node_max_z_);
    if (buf_node_left_) clReleaseMemObject(buf_node_left_);
    if (buf_node_right_) clReleaseMemObject(buf_node_right_);
    if (buf_node_tri_start_) clReleaseMemObject(buf_node_tri_start_);
    if (buf_node_tri_count_) clReleaseMemObject(buf_node_tri_count_);
    if (buf_tri_indices_) clReleaseMemObject(buf_tri_indices_);

    size_t node_size_float = num_nodes_ * sizeof(float);
    size_t node_size_uint = num_nodes_ * sizeof(uint32_t);
    size_t tri_idx_size = tri_indices.size() * sizeof(uint32_t);

    std::cerr << "[BVHRayTracer] BVH buffer sizes: nodes=" << (node_size_float / (1024.0)) << " KB, "
              << "indices=" << (tri_idx_size / (1024.0)) << " KB" << std::endl;

    // Use non-throwing error checks (macOS shared library exception issue)
    #define CL_CHECK_RETURN(call, msg) \
        do { \
            cl_int _err = (call); \
            if (_err != CL_SUCCESS) { \
                std::cerr << "OpenCL error in " << msg << ": " << clErrorString(_err) << std::endl; \
                return; \
            } \
        } while (0)

    std::cerr << "[BVHRayTracer] Creating BVH buffers..." << std::endl;

    // Create buffers
    buf_node_min_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_float, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_min_x");
    buf_node_min_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_float, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_min_y");
    buf_node_min_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_float, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_min_z");
    buf_node_max_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_float, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_max_x");
    buf_node_max_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_float, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_max_y");
    buf_node_max_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_float, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_max_z");
    buf_node_left_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_uint, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_left");
    buf_node_right_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_uint, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_right");
    buf_node_tri_start_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_uint, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_tri_start");
    buf_node_tri_count_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, node_size_uint, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer node_tri_count");
    buf_tri_indices_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, tri_idx_size, nullptr, &err);
    CL_CHECK_RETURN(err, "clCreateBuffer tri_indices");

    std::cerr << "[BVHRayTracer] Uploading BVH data..." << std::endl;

    // Upload data
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_min_x_, CL_FALSE, 0, node_size_float, min_x.data(), 0, nullptr, nullptr), "write node_min_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_min_y_, CL_FALSE, 0, node_size_float, min_y.data(), 0, nullptr, nullptr), "write node_min_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_min_z_, CL_FALSE, 0, node_size_float, min_z.data(), 0, nullptr, nullptr), "write node_min_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_max_x_, CL_FALSE, 0, node_size_float, max_x.data(), 0, nullptr, nullptr), "write node_max_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_max_y_, CL_FALSE, 0, node_size_float, max_y.data(), 0, nullptr, nullptr), "write node_max_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_max_z_, CL_FALSE, 0, node_size_float, max_z.data(), 0, nullptr, nullptr), "write node_max_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_left_, CL_FALSE, 0, node_size_uint, left.data(), 0, nullptr, nullptr), "write node_left");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_right_, CL_FALSE, 0, node_size_uint, right.data(), 0, nullptr, nullptr), "write node_right");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_tri_start_, CL_FALSE, 0, node_size_uint, tri_start.data(), 0, nullptr, nullptr), "write node_tri_start");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_node_tri_count_, CL_FALSE, 0, node_size_uint, tri_count.data(), 0, nullptr, nullptr), "write node_tri_count");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_tri_indices_, CL_FALSE, 0, tri_idx_size, tri_indices.data(), 0, nullptr, nullptr), "write tri_indices");

    std::cerr << "[BVHRayTracer] Finishing BVH upload..." << std::endl;

    err = clFinish(queue);
    CL_CHECK_RETURN(err, "clFinish");

    bvh_uploaded_ = true;
    std::cerr << "[BVHRayTracer] BVH upload complete" << std::endl;

    #undef CL_CHECK_RETURN
}

bool BVHRayTracer::ensureRayBuffers(size_t num_rays) {
    if (num_rays <= max_rays_) return true;

    cl_context ctx = context_->context();
    cl_int err;

    // Release existing ray buffers
    if (buf_ray_org_x_) clReleaseMemObject(buf_ray_org_x_);
    if (buf_ray_org_y_) clReleaseMemObject(buf_ray_org_y_);
    if (buf_ray_org_z_) clReleaseMemObject(buf_ray_org_z_);
    if (buf_ray_dir_x_) clReleaseMemObject(buf_ray_dir_x_);
    if (buf_ray_dir_y_) clReleaseMemObject(buf_ray_dir_y_);
    if (buf_ray_dir_z_) clReleaseMemObject(buf_ray_dir_z_);
    if (buf_ray_inv_dir_x_) clReleaseMemObject(buf_ray_inv_dir_x_);
    if (buf_ray_inv_dir_y_) clReleaseMemObject(buf_ray_inv_dir_y_);
    if (buf_ray_inv_dir_z_) clReleaseMemObject(buf_ray_inv_dir_z_);
    if (buf_ray_tmax_) clReleaseMemObject(buf_ray_tmax_);
    if (buf_ray_tmin_) clReleaseMemObject(buf_ray_tmin_);
    if (buf_ray_skip_) clReleaseMemObject(buf_ray_skip_);
    if (buf_ray_dn_) clReleaseMemObject(buf_ray_dn_);
    if (buf_ray_heps_) clReleaseMemObject(buf_ray_heps_);
    if (buf_occluded_) clReleaseMemObject(buf_occluded_);

    // Allocate with some headroom
    max_rays_ = num_rays + num_rays / 4;

    size_t buf_size_float = max_rays_ * sizeof(float);
    size_t buf_size_uchar = max_rays_ * sizeof(uint8_t);

    // Use non-throwing error check macro
    #define CL_CHECK_FALSE(call, msg) \
        do { \
            cl_int _err = (call); \
            if (_err != CL_SUCCESS) { \
                std::cerr << "OpenCL error in " << msg << ": " << clErrorString(_err) << std::endl; \
                return false; \
            } \
        } while (0)

    buf_ray_org_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_org_x");
    buf_ray_org_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_org_y");
    buf_ray_org_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_org_z");
    buf_ray_dir_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_dir_x");
    buf_ray_dir_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_dir_y");
    buf_ray_dir_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_dir_z");
    buf_ray_inv_dir_x_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_inv_dir_x");
    buf_ray_inv_dir_y_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_inv_dir_y");
    buf_ray_inv_dir_z_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_inv_dir_z");
    buf_ray_tmin_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_tmin");
    buf_ray_skip_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, max_rays_ * sizeof(cl_int), nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_skip_prim");
    buf_ray_dn_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_dir_dot_n");
    buf_ray_heps_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_h_eps");
    buf_ray_tmax_ = clCreateBuffer(ctx, CL_MEM_READ_ONLY, buf_size_float, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer ray_tmax");
    buf_occluded_ = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY, buf_size_uchar, nullptr, &err);
    CL_CHECK_FALSE(err, "clCreateBuffer occluded");

    #undef CL_CHECK_FALSE
    return true;
}

void BVHRayTracer::traceShadowRays(const std::vector<float>& ray_org_x,
                                   const std::vector<float>& ray_org_y,
                                   const std::vector<float>& ray_org_z,
                                   const std::vector<float>& ray_dir_x,
                                   const std::vector<float>& ray_dir_y,
                                   const std::vector<float>& ray_dir_z,
                                   std::vector<uint8_t>& occluded_out) {
    // Backward-compatible: infinite tmax (original shadow-to-sun semantics).
    std::vector<float> inf_tmax(ray_org_x.size(), 1e10f);
    traceShadowRays(ray_org_x, ray_org_y, ray_org_z,
                    ray_dir_x, ray_dir_y, ray_dir_z,
                    inf_tmax, occluded_out);
}

void BVHRayTracer::traceShadowRays(const std::vector<float>& ray_org_x,
                                   const std::vector<float>& ray_org_y,
                                   const std::vector<float>& ray_org_z,
                                   const std::vector<float>& ray_dir_x,
                                   const std::vector<float>& ray_dir_y,
                                   const std::vector<float>& ray_dir_z,
                                   const std::vector<float>& ray_tmax,
                                   std::vector<uint8_t>& occluded_out) {
    const size_t n = ray_org_x.size();
    traceShadowRaysEx(ray_org_x, ray_org_y, ray_org_z, ray_dir_x, ray_dir_y, ray_dir_z,
                      std::vector<float>(n, kLegacyTmin), ray_tmax,
                      std::vector<int32_t>(n, -1), std::vector<float>(n, 0.0f),
                      std::vector<float>(n, -1.0f), occluded_out);
}

size_t BVHRayTracer::kernelWorkGroupSize() const {
    // The kernel's own limit can be far below the device's (the OpenCL CPU
    // device reports 1 for this kernel); launching above it fails with
    // CL_INVALID_WORK_GROUP_SIZE on every batch.
    size_t wg = std::min(context_->deviceInfo().max_work_group_size, size_t(256));
    size_t kwg = 0;
    if (shadow_kernel_ &&
        clGetKernelWorkGroupInfo(shadow_kernel_, context_->device(), CL_KERNEL_WORK_GROUP_SIZE,
                                 sizeof(kwg), &kwg, nullptr) == CL_SUCCESS && kwg > 0) {
        wg = std::min(wg, kwg);
    }
    return std::max<size_t>(wg, 1);
}

namespace {
void dumpRays(const std::vector<float>& ox, const std::vector<float>& oy, const std::vector<float>& oz,
              const std::vector<float>& dx, const std::vector<float>& dy, const std::vector<float>& dz,
              const std::vector<float>& tmin, const std::vector<float>& tmax,
              const std::vector<int32_t>& skip, const std::vector<float>& dn,
              const std::vector<float>& heps, const std::vector<uint8_t>& occ, bool host) {
    RayDump& D = rayDump();
    std::lock_guard<std::mutex> lk(D.mu);
    const uint64_t call = D.calls++;
    std::ofstream rf(D.dir + "/rays.bin", std::ios::binary | std::ios::app);
    for (size_t i = 0; i < ox.size(); ++i) {
        if (!D.wants(skip[i])) continue;
        const uint64_t ri = i;
        const float f[8] = {ox[i], oy[i], oz[i], dx[i], dy[i], dz[i], tmin[i], tmax[i]};
        const float g[2] = {dn[i], heps[i]};
        const uint8_t b[2] = {occ[i], static_cast<uint8_t>(host)};
        rf.write(reinterpret_cast<const char*>(&call), 8);
        rf.write(reinterpret_cast<const char*>(&ri), 8);
        rf.write(reinterpret_cast<const char*>(f), sizeof(f));
        rf.write(reinterpret_cast<const char*>(&skip[i]), 4);
        rf.write(reinterpret_cast<const char*>(g), sizeof(g));
        rf.write(reinterpret_cast<const char*>(b), 2);
    }
}
}  // namespace

void BVHRayTracer::traceShadowRaysEx(const std::vector<float>& ray_org_x,
                                     const std::vector<float>& ray_org_y,
                                     const std::vector<float>& ray_org_z,
                                     const std::vector<float>& ray_dir_x,
                                     const std::vector<float>& ray_dir_y,
                                     const std::vector<float>& ray_dir_z,
                                     const std::vector<float>& ray_tmin,
                                     const std::vector<float>& ray_tmax,
                                     const std::vector<int32_t>& ray_skip_prim,
                                     const std::vector<float>& ray_dir_dot_n,
                                     const std::vector<float>& ray_h_eps,
                                     std::vector<uint8_t>& occluded_out) {
    if (!initialized_ || !mesh_uploaded_ || !bvh_uploaded_) {
        // Formerly returned every ray OCCLUDED and let the run finish with a
        // dark Sun. Fail instead; the caller decides (and says) what to do.
        throw std::runtime_error("BVHRayTracer not ready (init=" + std::to_string(initialized_) +
                                 ", mesh=" + std::to_string(mesh_uploaded_) +
                                 ", bvh=" + std::to_string(bvh_uploaded_) + ")");
    }
    const size_t n_in = ray_org_x.size();
    if (ray_tmin.size() != n_in || ray_tmax.size() != n_in || ray_skip_prim.size() != n_in ||
        ray_dir_dot_n.size() != n_in || ray_h_eps.size() != n_in) {
        throw std::runtime_error("BVHRayTracer::traceShadowRaysEx: per-ray array size mismatch");
    }
    if (host_kernel_) {
        const size_t nr = n_in;
        occluded_out.assign(nr, 0);
        if (nr == 0) return;
        auto inv = [](float d) { return 1.0f / (std::abs(d) < 1e-8f ? (d >= 0 ? 1e-8f : -1e-8f) : d); };
        std::vector<float> ix(nr), iy(nr), iz(nr);
        for (size_t r = 0; r < nr; ++r) { ix[r] = inv(ray_dir_x[r]); iy[r] = inv(ray_dir_y[r]); iz[r] = inv(ray_dir_z[r]); }
        const unsigned nth = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::thread> pool;
        const int nri = static_cast<int>(nr);
        for (unsigned t = 0; t < nth; ++t) {
            pool.emplace_back([&, t]() {
                for (int r = static_cast<int>(t); r < nri; r += static_cast<int>(nth)) {
                    host_kernel::g_global_id = r;
                    host_kernel::shadow_ray_kernel(
                        h_node_[0].data(), h_node_[1].data(), h_node_[2].data(),
                        h_node_[3].data(), h_node_[4].data(), h_node_[5].data(),
                        h_left_.data(), h_right_.data(), h_tstart_.data(), h_tcount_.data(),
                        h_tri_[0].data(), h_tri_[1].data(), h_tri_[2].data(),
                        h_tri_[3].data(), h_tri_[4].data(), h_tri_[5].data(),
                        h_tri_[6].data(), h_tri_[7].data(), h_tri_[8].data(), h_tidx_.data(),
                        ray_org_x.data(), ray_org_y.data(), ray_org_z.data(),
                        ray_dir_x.data(), ray_dir_y.data(), ray_dir_z.data(),
                        ix.data(), iy.data(), iz.data(), ray_tmax.data(), ray_tmin.data(),
                        ray_skip_prim.data(), ray_dir_dot_n.data(), ray_h_eps.data(),
                        occluded_out.data(), nri);
                }
            });
        }
        for (auto& th : pool) th.join();
        for (size_t r = 0; r < nr; ++r)
            if (occluded_out[r] > 1)
                throw std::runtime_error("BVHRayTracer: traversal stack overflow (host kernel)");
        stats_.rays_traced += nr;
        if (rayDump().on)
            dumpRays(ray_org_x, ray_org_y, ray_org_z, ray_dir_x, ray_dir_y, ray_dir_z, ray_tmin,
                     ray_tmax, ray_skip_prim, ray_dir_dot_n, ray_h_eps, occluded_out, true);
        return;
    }

    size_t num_rays = ray_org_x.size();
    if (num_rays == 0) {
        occluded_out.clear();
        return;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Ensure buffers are large enough
    if (!ensureRayBuffers(num_rays)) {
        throw std::runtime_error("BVHRayTracer: failed to allocate ray buffers");
    }

    // Compute inverse directions for slab intersection
    std::vector<float> inv_dir_x(num_rays), inv_dir_y(num_rays), inv_dir_z(num_rays);
    for (size_t i = 0; i < num_rays; ++i) {
        // Add small epsilon to avoid division by zero
        inv_dir_x[i] = 1.0f / (std::abs(ray_dir_x[i]) < 1e-8f ?
                               (ray_dir_x[i] >= 0 ? 1e-8f : -1e-8f) : ray_dir_x[i]);
        inv_dir_y[i] = 1.0f / (std::abs(ray_dir_y[i]) < 1e-8f ?
                               (ray_dir_y[i] >= 0 ? 1e-8f : -1e-8f) : ray_dir_y[i]);
        inv_dir_z[i] = 1.0f / (std::abs(ray_dir_z[i]) < 1e-8f ?
                               (ray_dir_z[i] >= 0 ? 1e-8f : -1e-8f) : ray_dir_z[i]);
    }

    cl_command_queue queue = context_->queue();
    size_t buf_size_float = num_rays * sizeof(float);

    // Use non-throwing error check macro
    #define CL_CHECK_RETURN(call, msg) \
        do { \
            cl_int _err = (call); \
            if (_err != CL_SUCCESS) { \
                throw std::runtime_error(std::string("BVHRayTracer: OpenCL error in ") + \
                                         msg + ": " + clErrorString(_err)); \
            } \
        } while (0)

    // Upload ray data
    auto upload_start = std::chrono::high_resolution_clock::now();

    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_org_x_, CL_FALSE, 0, buf_size_float, ray_org_x.data(), 0, nullptr, nullptr), "write ray_org_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_org_y_, CL_FALSE, 0, buf_size_float, ray_org_y.data(), 0, nullptr, nullptr), "write ray_org_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_org_z_, CL_FALSE, 0, buf_size_float, ray_org_z.data(), 0, nullptr, nullptr), "write ray_org_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_dir_x_, CL_FALSE, 0, buf_size_float, ray_dir_x.data(), 0, nullptr, nullptr), "write ray_dir_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_dir_y_, CL_FALSE, 0, buf_size_float, ray_dir_y.data(), 0, nullptr, nullptr), "write ray_dir_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_dir_z_, CL_FALSE, 0, buf_size_float, ray_dir_z.data(), 0, nullptr, nullptr), "write ray_dir_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_inv_dir_x_, CL_FALSE, 0, buf_size_float, inv_dir_x.data(), 0, nullptr, nullptr), "write ray_inv_dir_x");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_inv_dir_y_, CL_FALSE, 0, buf_size_float, inv_dir_y.data(), 0, nullptr, nullptr), "write ray_inv_dir_y");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_inv_dir_z_, CL_FALSE, 0, buf_size_float, inv_dir_z.data(), 0, nullptr, nullptr), "write ray_inv_dir_z");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_tmax_, CL_FALSE, 0, buf_size_float, ray_tmax.data(), 0, nullptr, nullptr), "write ray_tmax");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_tmin_, CL_FALSE, 0, buf_size_float, ray_tmin.data(), 0, nullptr, nullptr), "write ray_tmin");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_skip_, CL_FALSE, 0, num_rays * sizeof(cl_int), ray_skip_prim.data(), 0, nullptr, nullptr), "write ray_skip_prim");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_dn_, CL_FALSE, 0, buf_size_float, ray_dir_dot_n.data(), 0, nullptr, nullptr), "write ray_dir_dot_n");
    CL_CHECK_RETURN(clEnqueueWriteBuffer(queue, buf_ray_heps_, CL_FALSE, 0, buf_size_float, ray_h_eps.data(), 0, nullptr, nullptr), "write ray_h_eps");

    cl_int err = clFinish(queue);
    CL_CHECK_RETURN(err, "clFinish after ray upload");
    auto upload_end = std::chrono::high_resolution_clock::now();

    // Set kernel arguments
    int arg_idx = 0;
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_min_x_), "setarg node_min_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_min_y_), "setarg node_min_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_min_z_), "setarg node_min_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_max_x_), "setarg node_max_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_max_y_), "setarg node_max_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_max_z_), "setarg node_max_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_left_), "setarg node_left");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_right_), "setarg node_right");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_tri_start_), "setarg node_tri_start");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_node_tri_count_), "setarg node_tri_count");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v0_x_), "setarg tri_v0_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v0_y_), "setarg tri_v0_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v0_z_), "setarg tri_v0_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v1_x_), "setarg tri_v1_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v1_y_), "setarg tri_v1_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v1_z_), "setarg tri_v1_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v2_x_), "setarg tri_v2_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v2_y_), "setarg tri_v2_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_v2_z_), "setarg tri_v2_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_tri_indices_), "setarg tri_indices");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_org_x_), "setarg ray_org_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_org_y_), "setarg ray_org_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_org_z_), "setarg ray_org_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_dir_x_), "setarg ray_dir_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_dir_y_), "setarg ray_dir_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_dir_z_), "setarg ray_dir_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_inv_dir_x_), "setarg ray_inv_dir_x");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_inv_dir_y_), "setarg ray_inv_dir_y");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_inv_dir_z_), "setarg ray_inv_dir_z");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_tmax_), "setarg ray_tmax");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_tmin_), "setarg ray_tmin");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_skip_), "setarg ray_skip_prim");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_dn_), "setarg ray_dir_dot_n");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_ray_heps_), "setarg ray_h_eps");
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(cl_mem), &buf_occluded_), "setarg occluded");
    int num_rays_int = static_cast<int>(num_rays);
    CL_CHECK_RETURN(clSetKernelArg(shadow_kernel_, arg_idx++, sizeof(int), &num_rays_int), "setarg num_rays");

    // Launch kernel
    auto kernel_start = std::chrono::high_resolution_clock::now();

    size_t global_size = ((num_rays + work_group_size_ - 1) / work_group_size_) * work_group_size_;
    CL_CHECK_RETURN(clEnqueueNDRangeKernel(queue, shadow_kernel_, 1, nullptr,
                                    &global_size, &work_group_size_, 0, nullptr, nullptr), "enqueue kernel");
    err = clFinish(queue);
    CL_CHECK_RETURN(err, "clFinish after kernel");

    auto kernel_end = std::chrono::high_resolution_clock::now();

    // Read results
    occluded_out.resize(num_rays);
    CL_CHECK_RETURN(clEnqueueReadBuffer(queue, buf_occluded_, CL_TRUE, 0,
                                 num_rays * sizeof(uint8_t), occluded_out.data(),
                                 0, nullptr, nullptr), "read occluded");

    auto end_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < num_rays; ++i) {
        if (occluded_out[i] > 1) {
            throw std::runtime_error("BVHRayTracer: traversal stack overflow on ray " +
                                     std::to_string(i) + " (BVH deeper than the kernel stack)");
        }
    }

    if (rayDump().on)
        dumpRays(ray_org_x, ray_org_y, ray_org_z, ray_dir_x, ray_dir_y, ray_dir_z, ray_tmin,
                 ray_tmax, ray_skip_prim, ray_dir_dot_n, ray_h_eps, occluded_out, false);

    // Update stats
    stats_.rays_traced += num_rays;
    stats_.upload_time_ms += std::chrono::duration<double, std::milli>(upload_end - upload_start).count();
    stats_.kernel_time_ms += std::chrono::duration<double, std::milli>(kernel_end - kernel_start).count();
    stats_.total_time_ms += std::chrono::duration<double, std::milli>(end_time - start_time).count();

    if (gpu::gpuDispatchLogOn()) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "path=bvh device=\"%s\" rays=%zu upload_ms=%.3f kernel_ms=%.3f total_ms=%.3f",
                      getDeviceInfo().name.c_str(), num_rays,
                      std::chrono::duration<double, std::milli>(upload_end - upload_start).count(),
                      std::chrono::duration<double, std::milli>(kernel_end - kernel_start).count(),
                      std::chrono::duration<double, std::milli>(end_time - start_time).count());
        gpu::gpuDispatchLog(buf);
    }

    #undef CL_CHECK_RETURN
}

void BVHRayTracer::traceShadowRays(const std::vector<std::array<double, 3>>& origins,
                                   const std::vector<std::array<double, 3>>& directions,
                                   std::vector<uint8_t>& occluded_out) {
    size_t num_rays = origins.size();

    // Convert AoS to SoA
    std::vector<float> org_x(num_rays), org_y(num_rays), org_z(num_rays);
    std::vector<float> dir_x(num_rays), dir_y(num_rays), dir_z(num_rays);

    for (size_t i = 0; i < num_rays; ++i) {
        org_x[i] = static_cast<float>(origins[i][0]);
        org_y[i] = static_cast<float>(origins[i][1]);
        org_z[i] = static_cast<float>(origins[i][2]);
        dir_x[i] = static_cast<float>(directions[i][0]);
        dir_y[i] = static_cast<float>(directions[i][1]);
        dir_z[i] = static_cast<float>(directions[i][2]);
    }

    traceShadowRays(org_x, org_y, org_z, dir_x, dir_y, dir_z, occluded_out);
}

const OpenCLDeviceInfo& BVHRayTracer::getDeviceInfo() const {
    if (!context_) {
        static OpenCLDeviceInfo empty_info;
        return empty_info;
    }
    return context_->deviceInfo();
}

void BVHRayTracer::release() {
    if (shadow_kernel_) {
        clReleaseKernel(shadow_kernel_);
        shadow_kernel_ = nullptr;
    }
    if (program_) {
        clReleaseProgram(program_);
        program_ = nullptr;
    }

    // Release mesh buffers
    if (buf_tri_v0_x_) { clReleaseMemObject(buf_tri_v0_x_); buf_tri_v0_x_ = nullptr; }
    if (buf_tri_v0_y_) { clReleaseMemObject(buf_tri_v0_y_); buf_tri_v0_y_ = nullptr; }
    if (buf_tri_v0_z_) { clReleaseMemObject(buf_tri_v0_z_); buf_tri_v0_z_ = nullptr; }
    if (buf_tri_v1_x_) { clReleaseMemObject(buf_tri_v1_x_); buf_tri_v1_x_ = nullptr; }
    if (buf_tri_v1_y_) { clReleaseMemObject(buf_tri_v1_y_); buf_tri_v1_y_ = nullptr; }
    if (buf_tri_v1_z_) { clReleaseMemObject(buf_tri_v1_z_); buf_tri_v1_z_ = nullptr; }
    if (buf_tri_v2_x_) { clReleaseMemObject(buf_tri_v2_x_); buf_tri_v2_x_ = nullptr; }
    if (buf_tri_v2_y_) { clReleaseMemObject(buf_tri_v2_y_); buf_tri_v2_y_ = nullptr; }
    if (buf_tri_v2_z_) { clReleaseMemObject(buf_tri_v2_z_); buf_tri_v2_z_ = nullptr; }

    // Release BVH buffers
    if (buf_node_min_x_) { clReleaseMemObject(buf_node_min_x_); buf_node_min_x_ = nullptr; }
    if (buf_node_min_y_) { clReleaseMemObject(buf_node_min_y_); buf_node_min_y_ = nullptr; }
    if (buf_node_min_z_) { clReleaseMemObject(buf_node_min_z_); buf_node_min_z_ = nullptr; }
    if (buf_node_max_x_) { clReleaseMemObject(buf_node_max_x_); buf_node_max_x_ = nullptr; }
    if (buf_node_max_y_) { clReleaseMemObject(buf_node_max_y_); buf_node_max_y_ = nullptr; }
    if (buf_node_max_z_) { clReleaseMemObject(buf_node_max_z_); buf_node_max_z_ = nullptr; }
    if (buf_node_left_) { clReleaseMemObject(buf_node_left_); buf_node_left_ = nullptr; }
    if (buf_node_right_) { clReleaseMemObject(buf_node_right_); buf_node_right_ = nullptr; }
    if (buf_node_tri_start_) { clReleaseMemObject(buf_node_tri_start_); buf_node_tri_start_ = nullptr; }
    if (buf_node_tri_count_) { clReleaseMemObject(buf_node_tri_count_); buf_node_tri_count_ = nullptr; }
    if (buf_tri_indices_) { clReleaseMemObject(buf_tri_indices_); buf_tri_indices_ = nullptr; }

    // Release ray buffers
    if (buf_ray_org_x_) { clReleaseMemObject(buf_ray_org_x_); buf_ray_org_x_ = nullptr; }
    if (buf_ray_org_y_) { clReleaseMemObject(buf_ray_org_y_); buf_ray_org_y_ = nullptr; }
    if (buf_ray_org_z_) { clReleaseMemObject(buf_ray_org_z_); buf_ray_org_z_ = nullptr; }
    if (buf_ray_dir_x_) { clReleaseMemObject(buf_ray_dir_x_); buf_ray_dir_x_ = nullptr; }
    if (buf_ray_dir_y_) { clReleaseMemObject(buf_ray_dir_y_); buf_ray_dir_y_ = nullptr; }
    if (buf_ray_dir_z_) { clReleaseMemObject(buf_ray_dir_z_); buf_ray_dir_z_ = nullptr; }
    if (buf_ray_inv_dir_x_) { clReleaseMemObject(buf_ray_inv_dir_x_); buf_ray_inv_dir_x_ = nullptr; }
    if (buf_ray_inv_dir_y_) { clReleaseMemObject(buf_ray_inv_dir_y_); buf_ray_inv_dir_y_ = nullptr; }
    if (buf_ray_inv_dir_z_) { clReleaseMemObject(buf_ray_inv_dir_z_); buf_ray_inv_dir_z_ = nullptr; }
    if (buf_ray_tmax_) { clReleaseMemObject(buf_ray_tmax_); buf_ray_tmax_ = nullptr; }
    if (buf_ray_tmin_) { clReleaseMemObject(buf_ray_tmin_); buf_ray_tmin_ = nullptr; }
    if (buf_ray_skip_) { clReleaseMemObject(buf_ray_skip_); buf_ray_skip_ = nullptr; }
    if (buf_ray_dn_) { clReleaseMemObject(buf_ray_dn_); buf_ray_dn_ = nullptr; }
    if (buf_ray_heps_) { clReleaseMemObject(buf_ray_heps_); buf_ray_heps_ = nullptr; }
    if (buf_occluded_) { clReleaseMemObject(buf_occluded_); buf_occluded_ = nullptr; }

    context_.reset();
    initialized_ = false;
    mesh_uploaded_ = false;
    bvh_uploaded_ = false;
    max_rays_ = 0;
}

}  // namespace opencl
}  // namespace gpu

#endif  // PTM_GPU_OPENCL
