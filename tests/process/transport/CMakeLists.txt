# CMakeLists.txt for transport process tests
#
# Follows the same layout as the other process tests (tests/process/<name>/),
# auto-discovered by tests/CMakeLists.txt. The transport numerics library is
# built with `-fdefault-real-8`; the test itself is not, so the test source
# declares every array passed to a transport routine as `real(r8)` explicitly
# (see the header note in unit/test_transport_unit.F90).

# Test executables
set(
  TEST_SOURCES
  unit/test_transport_unit.F90
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
