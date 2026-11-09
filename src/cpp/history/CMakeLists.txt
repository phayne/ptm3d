# History marcher (M1, 2026-11-09): pause / checkpoint / resume, and its byte-identical gate.
# Not linked into ptm_thermal_sim.
add_library(history_marcher STATIC
    history_checkpoint.cpp
)
target_include_directories(history_marcher PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/cpp
)

add_executable(test_history_checkpoint test_history_checkpoint.cpp)
target_include_directories(test_history_checkpoint PRIVATE ${CMAKE_SOURCE_DIR}/src/cpp)
target_link_libraries(test_history_checkpoint PRIVATE history_marcher radiation thermal mesh_utils)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_history_checkpoint PRIVATE OpenMP::OpenMP_CXX)
endif()
if(ENABLE_GPU AND TARGET ptm_gpu)
    target_link_libraries(test_history_checkpoint PRIVATE ptm_gpu)
    target_compile_definitions(test_history_checkpoint PRIVATE PTM_GPU_OPENCL)
endif()
set_target_properties(test_history_checkpoint PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
