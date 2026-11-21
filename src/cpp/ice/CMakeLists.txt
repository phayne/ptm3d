add_library(ice_table
    ice_table.cpp
    ice_retention.cpp
)

target_include_directories(ice_table
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(ice_table
    PUBLIC
    ${OpenMP_CXX_LIBRARIES}
)

set_target_properties(ice_table PROPERTIES
    OUTPUT_NAME ice_table
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Analytic gate for the ice-stability criterion (surface rule, S&T 2007 Eq. 11,
# threshold in physical units). Synthetic accumulators only: no mesh, no SPICE,
# no GPU; under a second.
add_executable(test_ice_stability_criterion test_ice_stability_criterion.cpp)
target_link_libraries(test_ice_stability_criterion PRIVATE ice_table)
set_target_properties(test_ice_stability_criterion PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Analytic gates for the ice retention criterion (isothermal and two-layer
# closed forms, underflow, grid convergence). Synthetic columns only.
add_executable(test_ice_retention test_ice_retention.cpp)
target_link_libraries(test_ice_retention PRIVATE ice_table)
set_target_properties(test_ice_retention PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
