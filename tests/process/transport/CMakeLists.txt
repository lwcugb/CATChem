# CMakeLists.txt for transport process tests
#
# Follows the same layout as the other process tests (tests/process/<name>/),
# auto-discovered by tests/CMakeLists.txt. The transport numerics library is
# built with `-fdefault-real-8`; the test itself is not, so the test source
# declares every array passed to a transport routine as `real(r8)` explicitly
# (see the header note in unit/test_transport_unit.F90).

# Stage the shared default configuration + species files the integration test
# loads at runtime (its WORKING_DIRECTORY is ${CMAKE_CURRENT_BINARY_DIR}). The
# test reads CATChem_new_config_standalone.yml, which in turn references
# CATChem_species.yml; copy every default *.yml so those (and any files they
# reference) are present next to the executable.
file(GLOB TEST_DATA_FILES "${CATCHEM_TESTS_CONFIG_DIR}/*.yml")
foreach(test_data_file ${TEST_DATA_FILES})
  get_filename_component(filename ${test_data_file} NAME)
  configure_file(
    "${test_data_file}"
    "${CMAKE_CURRENT_BINARY_DIR}/${filename}"
    COPYONLY
  )
endforeach()

# Test executables
set(
  TEST_SOURCES
  unit/test_transport_unit.F90
  integration/test_transport_integration.F90
)

foreach(test_source ${TEST_SOURCES})
  get_filename_component(test_name ${test_source} NAME_WE)

  add_executable(${test_name} ${test_source})

  # Link against the transport process library + core + test harness
  target_link_libraries(
    ${test_name}
    CATChem_process_transport
    CATChem_core
    testing
  )

  # Set Fortran module directory and force Fortran linker
  set_target_properties(
    ${test_name}
    PROPERTIES
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
      LINKER_LANGUAGE Fortran
  )

  # Transport library .mod files are exported to ${CMAKE_BINARY_DIR}/include
  target_include_directories(
    ${test_name}
    PRIVATE "${CMAKE_BINARY_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/modules"
  )

  add_dependencies(${test_name} CATChem_process_transport)

  # Add to test suite
  add_test(NAME ${test_name} COMMAND ${test_name})

  set_tests_properties(
    ${test_name}
    PROPERTIES TIMEOUT 60 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )
endforeach()

# Convenience target to run the transport unit tests
add_custom_target(
  run_transport_unit
  COMMAND test_transport_unit
  DEPENDS test_transport_unit
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Running unit tests for transport process"
)

# Convenience target to run the transport integration test
add_custom_target(
  run_transport_integration
  COMMAND test_transport_integration
  DEPENDS test_transport_integration
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Running the end-to-end integration test for the transport process"
)
