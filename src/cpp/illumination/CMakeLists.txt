# Illumination Model Library

# Add all source files
file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h *.hpp)

# Create the illumination model library
add_library(illumination_model
    illumination_model.cpp
)

# Include directories
target_include_directories(illumination_model
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CSPICE_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}/include
)

# Link dependencies
target_link_libraries(illumination_model
    PUBLIC
    constants
    mesh_utils
    orbital
    ${CSPICE_LIBRARY}
    ${EMBREE_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# GPU support for BVH ray tracing
if(ENABLE_GPU AND TARGET ptm_gpu)
    target_link_libraries(illumination_model PUBLIC ptm_gpu)
    target_compile_definitions(illumination_model PUBLIC PTM_GPU_OPENCL)
endif()
# ---------------------------------------------------------------------------
# Shadow depth: the umbral-ceiling post-solve pass.
#
# A SEPARATE, dependency-light library on purpose. Like flux_cubes it takes its
# ray casting as a callback and links no Embree, no CSPICE and no illumination
# model, so the standalone falsifier can exercise it against synthetic geometry
# with closed-form answers. Linking it into illumination_model would make the
# test drag in a scene lifetime that only exists inside a running simulation.
# ---------------------------------------------------------------------------
add_library(shadow_depth STATIC
    shadow_depth_pass.cpp
)
target_include_directories(shadow_depth PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/cpp
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(shadow_depth PUBLIC OpenMP::OpenMP_CXX)
endif()

# ---------------------------------------------------------------------------
# Horizon classifier (history marcher M0, 2026-11-09): per-face solar horizon
# envelope + classify/fall-through census + keyed .hzc file + audit selection,
# and the Embree scene/predicate it certifies against. NOT linked into
# ptm_thermal_sim (no default can move). IEEE arithmetic on purpose: the
# classification compares disk extents with bracket ends found by the same
# arithmetic, so -ffast-math / FMA contraction would flip boundary cases and
# break bit-parity with the measurement tool and across ISAs.
# ---------------------------------------------------------------------------
add_library(horizon_classifier STATIC
    horizon_classifier.cpp
    horizon_classifier_scene.cpp
)
target_include_directories(horizon_classifier PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/cpp
    ${EMBREE_INCLUDE_DIRS}
)
target_link_libraries(horizon_classifier PUBLIC ${EMBREE_LIBRARIES})
target_compile_options(horizon_classifier PRIVATE -fno-fast-math -ffp-contract=off)
