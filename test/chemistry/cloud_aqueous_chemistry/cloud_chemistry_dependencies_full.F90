! Mocked dependencies of cloud aqueous chemistry module
! modified to include randomly generated data for species
! not included in the configuration used to generate the
! snapshot file.
module chemistry_test_data
  implicit none
  public
  logical :: do_debug_logging = .false.
  integer :: debug_column = 13
  integer :: debug_level = 26
end module chemistry_test_data

module shr_kind_mod
  implicit none
  public :: shr_kind_r8
  integer, parameter :: shr_kind_r8 = kind(0.0d0)
end module shr_kind_mod

module spmd_utils
  implicit none
  public :: masterproc
  logical :: masterproc = .true.
end module spmd_utils

module cam_logfile
  implicit none
  public :: iulog
  integer, parameter :: iulog = 6
end module cam_logfile

module physics_buffer
  implicit none
  public :: physics_buffer_desc, pbuf_get_index, pbuf_add_field, dtype_r8
  integer, parameter :: dtype_r8 = 1
  integer, parameter :: pbuf_get_index = 1
  integer, parameter :: pbuf_add_field = 1
  type :: physics_buffer_desc
  end type physics_buffer_desc
end module physics_buffer

module physics_types
  implicit none
  public :: physics_state
  type :: physics_state
  end type physics_state
end module physics_types

module mo_chem_utls
  implicit none
  public :: get_spc_ndx, get_inv_ndx
contains
  integer function get_spc_ndx(spc_name)
    character(len=*), intent(in) :: spc_name
    select case (trim(adjustl(spc_name)))
    case ('SO2')
      get_spc_ndx = 18
    case ('H2O2')
      get_spc_ndx = 7
    case ('SO4')
      get_spc_ndx = 0
    case ('H2SO4')
      get_spc_ndx = 8
    case ('NH3')
      get_spc_ndx = 1
    case ('HNO3')
      get_spc_ndx = 2
    case ('MSA')
      get_spc_ndx = 3
    case default
      get_spc_ndx = -1
    end select
  end function get_spc_ndx
  integer function get_inv_ndx(inv_name)
    character(len=*), intent(in) :: inv_name
    select case (trim(adjustl(inv_name)))
    case ('O3')
      get_inv_ndx = 4
    case ('HO2')
      get_inv_ndx = 6
    case default
      get_inv_ndx = -1
    end select
  end function get_inv_ndx
end module mo_chem_utls

module carma_flags_mod
  implicit none
  public :: carma_do_cloudborne
  logical, parameter :: carma_do_cloudborne = .false.
end module carma_flags_mod

module phys_control
  implicit none
  public :: phys_getopts, cam_chempkg_is
contains
  subroutine phys_getopts(prog_modal_aero_out, history_aerosol_out)
    logical, optional, intent(out) :: prog_modal_aero_out
    logical, optional, intent(out) :: history_aerosol_out
    if (present(prog_modal_aero_out)) then
      prog_modal_aero_out = .true.
    end if
    if (present(history_aerosol_out)) then
      history_aerosol_out = .false.
    end if
  end subroutine phys_getopts
  logical function cam_chempkg_is(pkg)
    character(len=*), intent(in) :: pkg
    cam_chempkg_is = .false.
  end function cam_chempkg_is
end module phys_control

module ppgrid
  implicit none
  public :: pcols, pver
  integer, parameter :: pcols = 16
  integer, parameter :: pver = 32
end module ppgrid

module chem_mods
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  public :: gas_pcnst, nfs, adv_mass
  integer, parameter :: gas_pcnst = 26
  integer, parameter :: nfs = 1
  real(kind=r8), parameter :: adv_mass(*) = (/ &
          12.010999999999999_r8, 12.010999999999999_r8, 62.132399999999997_r8, &
          135.06403900000001_r8, 135.06403900000001_r8, 135.06403900000001_r8, 34.013599999999997_r8, &
          98.078400000000002_r8, 58.442467999999998_r8, 58.442467999999998_r8, 58.442467999999998_r8, &
          1.0074000000000001_r8, 1.0074000000000001_r8, 1.0074000000000001_r8, 1.0074000000000001_r8, &
          12.010999999999999_r8, 12.010999999999999_r8, 64.064800000000005_r8, 115.10733999999999_r8, &
          115.10733999999999_r8, 115.10733999999999_r8, 12.010999999999999_r8, 12.010999999999999_r8, &
          12.010999999999999_r8, 12.010999999999999_r8, 18.014199999999999_r8 /)
end module chem_mods

module physconst
  use shr_kind_mod, only : shr_kind_r8
  implicit none
  public :: mwdry, gravit
  real(kind=shr_kind_r8), parameter :: mwdry = -999.0_shr_kind_r8
  real(kind=shr_kind_r8), parameter :: gravit = 9.7976399999999995_shr_kind_r8
end module physconst

module mo_constants
  use shr_kind_mod, only : shr_kind_r8
  implicit none
  public :: pi
  real(kind=shr_kind_r8), parameter :: pi = 3.1415926535897931_shr_kind_r8
end module mo_constants

module cam_abortutils
  implicit none
  public :: endrun
contains
  subroutine endrun(msg)
    character(len=*), intent(in) :: msg
    write(*,*) msg
    stop 3
  end subroutine endrun
end module cam_abortutils

module modal_aero_data
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  public :: ntot_amode, modeptr_accum, lptr_so4_cw_amode, lptr_msa_cw_amode, &
            numptrcw_amode, lptr_nh4_cw_amode, cnst_name_cw, specmw_so4_amode
  integer, parameter :: ntot_amode = 4
  integer, parameter :: modeptr_accum = 1
  integer, parameter :: lptr_so4_cw_amode(*) = (/ 28, 29, 30, -999888777 /)
  integer, parameter :: lptr_msa_cw_amode(*) = (/ -999888777, -999888777, -999888777, -999888777 /)
  integer, parameter :: numptrcw_amode(*) = (/ 21, 22, 23, 24 /)
  integer, parameter :: lptr_nh4_cw_amode(*) = (/ -999888777, -999888777, -999888777, -999888777 /)
  character(len=6), parameter :: cnst_name_cw(*) = (/ &
    ' bc_c1', ' bc_c4', 'dst_c1', 'dst_c2', 'dst_c3', 'ncl_c1', 'ncl_c2', 'ncl_c3', 'num_c1', &
    'num_c2', 'num_c3', 'num_c4', 'pom_c1', 'pom_c4', 'so4_c1', 'so4_c2', 'so4_c3', 'soa_c1', &
    'soa_c2' /)
  real(kind=r8), parameter :: specmw_so4_amode = 115.10733999999999_r8
end module modal_aero_data

