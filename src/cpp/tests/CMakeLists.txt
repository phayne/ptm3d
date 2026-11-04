# CMakeLists.txt for the standalone cross-module tests.
#
# ⛔ THIS DIRECTORY WAS NOT IN THE BUILD BEFORE 2026-09-14. Three sources sat
# here (test_gpu_guardrails.cpp, test_machine_profile.cpp,
# test_riseeps_occlusion.cpp) with no CMake target anywhere in the tree; the
# only thing that ever compiled one is a hand-written `c++` line in INSTALL.md
# section for test_machine_profile. They are deliberately LEFT as they are here
# -- adding targets for three tests nobody has compiled in months is a separate
# change with its own gating -- but the gap is real and is tracked
# as an open issue.

# The shadow-depth falsifier: synthetic meshes, closed-form umbral ceilings, no
# SPICE, no GPU, no Embree, no engine. Links only the pass, which is why the
# pass takes its ray casting as a callback.
add_executable(test_shadow_depth test_shadow_depth.cpp)
target_link_libraries(test_shadow_depth PRIVATE shadow_depth)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_shadow_depth PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(test_shadow_depth PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# PTM-780: the open-mesh body sphere R_occ must include the occluder ring.
# Synthetic high core on a low ring; drives the real IlluminationModel through
# both the batch (production) and single-face kernels from a precomputed
# positions file. No SPICE kernels, no GPU. Fails on the pre-fix engine.
add_executable(test_rocc_core_plus_ring test_rocc_core_plus_ring.cpp)
target_link_libraries(test_rocc_core_plus_ring PRIVATE
    illumination_model mesh_utils thermal constants)
set_target_properties(test_rocc_core_plus_ring PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# PTM-794: GPU shadow-ray tracer self-check against a float64 brute-force ground
# truth on synthetic meshes (no data files). Built in every GPU build so a user
# can check new GPU hardware before turning on illumination.use_gpu; run it
# directly or as `ptm gpu-selfcheck`. In a build without OpenCL it compiles to
# a stub that reports SKIP. Like bvh_opencl.cpp it opts out of the project-wide
# -ffast-math: the ground truth's tie tolerances assume IEEE arithmetic.
if(ENABLE_GPU AND TARGET ptm_gpu)
    add_executable(ptm_gpu_selfcheck ptm_gpu_selfcheck.cpp)
    target_link_libraries(ptm_gpu_selfcheck PRIVATE ptm_gpu)
    set_source_files_properties(ptm_gpu_selfcheck.cpp PROPERTIES
        COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off")
    set_target_properties(ptm_gpu_selfcheck PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
    install(TARGETS ptm_gpu_selfcheck RUNTIME DESTINATION bin)
endif()
