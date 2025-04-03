! Mocked dependencies of cloud aqueous chemistry module
module chemistry_test_data
  implicit none
  public
  logical :: do_debug_logging = .false.
  integer :: debug_column = -1
  integer :: debug_level = -1
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
    case ('O3')
      get_spc_ndx = 1
    case ('SO2')
      get_spc_ndx = 83
    case ('NH3')
      get_spc_ndx = 86
    case ('HNO3')
      get_spc_ndx = 8
    case ('H2O2')
      get_spc_ndx = 14
    case ('HO2')
      get_spc_ndx = 13
    case ('SO4')
      get_spc_ndx = 85
    case default
      get_spc_ndx = -1
    end select
  end function get_spc_ndx
  integer function get_inv_ndx(inv_name)
    character(len=*), intent(in) :: inv_name
    select case (trim(adjustl(inv_name)))
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
      prog_modal_aero_out = .false.
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
  integer, parameter :: pver = 26
end module ppgrid

module chem_mods
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  public :: gas_pcnst, nfs, adv_mass
  integer, parameter :: gas_pcnst = 103
  integer, parameter :: nfs = 4
  real(kind=r8), parameter :: adv_mass(*) = (/ &
    47.998199999999997, 15.999400000000000, 15.999400000000000, 44.012880000000003, &
    30.006139999999998, 46.005540000000003, 62.004939999999998, 63.012340000000002, &
    79.011740000000003, 108.01048000000000, 2.0148000000000001, 17.006799999999998, &
    33.006200000000000, 34.013599999999997, 16.040600000000001, 28.010400000000001, &
    47.031999999999996, 48.039400000000001, 30.025200000000002, 32.039999999999999, &
    46.065800000000003, 28.051600000000001, 61.057800000000000, 77.057199999999995, &
    60.050400000000003, 60.050400000000003, 30.066400000000002, 61.057800000000000, &
    62.065199999999997, 44.051000000000002, 75.042400000000001, 76.049800000000005, &
    42.077399999999997, 44.092199999999998, 75.083600000000004, 76.090999999999994, &
    91.082999999999998, 92.090400000000002, 58.076799999999999, 89.068200000000004, &
    90.075599999999994, 56.103200000000001, 105.10880000000000, 72.102599999999995, &
    103.09399999999999, 104.10140000000000, 72.143799999999999, 103.13520000000000, &
    104.14260000000000, 68.114199999999997, 117.11980000000000, 118.12720000000000, &
    70.087800000000001, 70.087800000000001, 119.09340000000000, 120.10080000000001, &
    101.07920000000000, 100.11300000000000, 74.076200000000000, 72.061400000000006, &
    149.11859999999999, 150.12600000000000, 136.22839999999999, 185.23400000000001, &
    186.24140000000000, 92.136200000000002, 108.13560000000000, 173.14060000000001, &
    174.14800000000000, 190.14740000000000, 98.098200000000006, 58.035600000000002, &
    121.04794000000000, 119.07434000000001, 147.08474000000001, 162.11794000000000, &
    147.12594000000001, 12.010999999999999, 12.010999999999999, 12.010999999999999, &
    12.010999999999999, 144.13200000000001, 64.064800000000005, 62.132399999999997, &
    96.063599999999994, 17.028939999999999, 18.036339999999999, 80.041280000000000, &
    58.442467999999998, 58.442467999999998, 58.442467999999998, 58.442467999999998, &
    135.06403900000001, 135.06403900000001, 135.06403900000001, 135.06403900000001, &
    222.00000000000000, 207.19999999999999, 27.025140000000000, 41.050939999999997, &
    26.036799999999999, 46.024600000000000, 63.031399999999998 /)
end module chem_mods

module physconst
  use shr_kind_mod, only : shr_kind_r8
  implicit none
  public :: mwdry, gravit
  real(kind=shr_kind_r8), parameter :: mwdry = 28.966000000000001_shr_kind_r8
  real(kind=shr_kind_r8), parameter :: gravit = 9.8061600000000002_shr_kind_r8
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

module rad_constituents
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  public :: rad_cnst_get_info, rad_cnst_get_info_by_bin, rad_cnst_get_bin_props_by_idx
contains
  subroutine rad_cnst_get_info(list_idx, nbins)
    integer, intent(in) :: list_idx
    integer, optional, intent(out) :: nbins
    if(present(nbins)) then
      nbins = 0
    end if
  end subroutine rad_cnst_get_info
  subroutine rad_cnst_get_info_by_bin(list_idx, bin_idx, nspec, bin_name)
    integer, intent(in) :: list_idx
    integer, intent(in) :: bin_idx
    integer, optional, intent(out) :: nspec
    character(len=*), optional, intent(out) :: bin_name
    if(present(bin_name)) then
      bin_name = 'something'
    end if
    if(present(nspec)) then
      nspec = 0
    end if
  end subroutine rad_cnst_get_info_by_bin
  subroutine rad_cnst_get_bin_props_by_idx(list_idx, bin_idx, spec_idx, spectype)
    integer, intent(in) :: list_idx
    integer, intent(in) :: bin_idx
    integer, intent(in) :: spec_idx
    character(len=*), optional, intent(out) :: spectype
    if(present(spectype)) then
      spectype = 'sulfate'
    end if
  end subroutine rad_cnst_get_bin_props_by_idx
end module rad_constituents

module aerosol_properties_mod
  implicit none
  public :: aero_name_len
  integer, parameter :: aero_name_len = 32
end module aerosol_properties_mod

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

module carma_intr
  implicit none
  public :: carma_get_group_by_name, carma_get_dry_radius
contains
  subroutine carma_get_group_by_name(short_name, igroup, rc)
    character(len=*), intent(in) :: short_name
    integer, intent(out) :: igroup
    integer, intent(out) :: rc
    igroup = 1
    rc = 0
  end subroutine carma_get_group_by_name
  subroutine carma_get_dry_radius(state, igroup, ibin, dryr, rho, rc)
    use physics_types, only: physics_state
    use shr_kind_mod, only: r8 => shr_kind_r8
    type(physics_state), intent(in) :: state
    integer, intent(in) :: igroup
    integer, intent(in) :: ibin
    real(kind=r8), intent(out) :: dryr(:,:)
    real(kind=r8), intent(out) :: rho(:,:)
    integer, intent(out) :: rc
    dryr(:,:) = 0.0_r8 ! m
    rho(:,:) = 0.0_r8 ! kg m-3
    rc = 0
  end subroutine carma_get_dry_radius
end module carma_intr


