set(_srcs CCPr_Megan_Mod.F90 ccpr_megan_common_mod.F90 ccpr_scheme_megan_mod.F90)

set(_lib CATChem_process_megan)

add_library(${_lib} ${_srcs})
target_link_libraries(${_lib} PUBLIC CATChem_core)
set_target_properties(${_lib} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
