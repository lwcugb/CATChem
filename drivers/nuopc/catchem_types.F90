!> \file catchem_types.F90
!> \details
!> Defines the core data structures used by CATCHEM chemistry model including
!> arrays for meteorological states, chemical species concentrations, emissions,
!> and diagnostic outputs. Provides the fundamental data container types needed
!> for chemistry calculations.
!>
!> \author Barry Baker
!>
!> \date 11/2024
!>
!> \defgroup catchem_nuopc_group CATChem NUOPC Interface
!> \brief NUOPC interface drivers and data types for CATChem
!> \ingroup catchem
!>
!> This group contains all NUOPC-compliant interface modules and data types
!> for integrating CATChem with the NUOPC framework.
!!!>

module catchem_types

  use CATChem, only: GridStateType, MetStateType, ChemStateType, EmisStateType, DiagStateType

  implicit none
  private

  public :: catchem_container_type

  type :: catchem_container_type
    integer :: im
    type(GridStateType)  :: GridState
    type(MetStateType),  allocatable :: MetState(:)
    type(ChemStateType), allocatable :: ChemState(:)
    type(EmisStateType), allocatable :: EmisState(:)
    type(DiagStateType), allocatable :: DiagState(:)
  end type

end module catchem_types
