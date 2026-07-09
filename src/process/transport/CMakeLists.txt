# Transport Process CMakeLists.txt
# Horizontal (and later vertical) tracer transport.
#
# Numerics: vendored GFDL FV3 tp_core kernel (kernels/fv3_tp_core_mod.F90, LGPL v3)
# built in-tree with a trimmed grid-metric types module (no FMS dependency).

# Kernel sources (vendored FV3 flux operator + trimmed metric types)
set(
  TRANSPORT_KERNEL_SOURCES
  kernels/fv3_grid_types_mod.F90
  kernels/fv3_tp_core_mod.F90
  kernels/fv3_vremap_mod.F90
  kernels/fv3_pfix_mod.F90
)

# Process sources (CATChem-owned driver + factory)
set(
  TRANSPORT_PROCESS_SOURCES
  TransportHalo_Mod.F90
  TransportGridMetrics_Mod.F90
  TransportMassFlux_Mod.F90
  ProcessTransportInterface_Mod.F90
  TransportProcessCreator_Mod.F90
)

set(TRANSPORT_ALL_SOURCES ${TRANSPORT_KERNEL_SOURCES} ${TRANSPORT_PROCESS_SOURCES})

# Create the transport process library
set(_lib CATChem_process_transport)
add_library(${_lib} ${TRANSPORT_ALL_SOURCES})

# Link with required libraries
target_link_libraries(${_lib} PUBLIC CATChem_core)

# Set module directory
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# Include directories
target_include_directories(${_lib} PRIVATE ${CMAKE_BINARY_DIR}/include)

# Compiler-specific flags.
# NOTE: -r8 / -fdefault-real-8 is REQUIRED here: the vendored FV3 kernel uses
# bare `real`, which must be 8-byte to interoperate with CATChem's real(fp)=8.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(${_lib} PRIVATE -fdefault-real-8 -fdefault-double-8)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(${_lib} PRIVATE -r8)
endif()

# Add preprocessor definitions
target_compile_definitions(
  ${_lib}
  PRIVATE PROCESS_NAME="transport" PROCESS_VERSION="0.1.0"
)

# Install targets
install(
  TARGETS ${_lib}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

# Install module files
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.mod"
)

message(STATUS "Configured Transport process:")
message(STATUS "  Process type: transport")
message(STATUS "  Kernel: vendored FV3 tp_core + fv_operators mappm (LGPL v3)")
