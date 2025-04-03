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
      get_spc_ndx = 122
    case ('SO2')
      get_spc_ndx = 137
    case ('NH3')
      get_spc_ndx = 112
    case ('HNO3')
      get_spc_ndx = 84
    case ('H2O2')
      get_spc_ndx = 74
    case ('HO2')
      get_spc_ndx = 176
    case ('H2SO4')
      get_spc_ndx = 75
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
  logical, parameter :: carma_do_cloudborne = .true.
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
  integer, parameter :: pver = 32
end module ppgrid

module chem_mods
  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none
  public :: gas_pcnst, nfs, adv_mass
  integer, parameter :: gas_pcnst = 202
  integer, parameter :: nfs = 3
  real(kind=r8), parameter :: adv_mass(*) = (/ &
    133.14134000000001_r8, 104.14260000000000_r8, 28.010400000000001_r8, 204.34260000000000_r8, &
    78.110399999999998_r8, 160.12219999999999_r8, 126.10860000000000_r8, 98.098200000000006_r8, &
    84.072400000000002_r8, 98.098200000000006_r8, 98.098200000000006_r8, 112.12400000000000_r8, &
    72.143799999999999_r8, 56.103200000000001_r8, 79.903999999999996_r8, 115.35670000000000_r8, &
    95.903400000000005_r8, 141.90894000000000_r8, 99.716849999999994_r8, 106.12080000000000_r8, &
    124.13500000000001_r8, 26.036799999999999_r8, 28.051600000000001_r8, 46.065800000000003_r8, &
    62.065199999999997_r8, 30.066400000000002_r8, 42.077399999999997_r8, 76.090999999999994_r8, &
    44.092199999999998_r8, 110.10920000000000_r8, 153.82180000000000_r8, 165.36450600000001_r8, &
    148.91021000000001_r8, 137.36750300000000_r8, 187.37531000000001_r8, 170.92101299999999_r8, &
    154.46671599999999_r8, 120.91320600000000_r8, 173.83380000000000_r8, 30.025200000000002_r8, &
    94.937200000000004_r8, 133.40230000000000_r8, 44.051000000000002_r8, 50.485900000000001_r8, &
    41.050939999999997_r8, 58.076799999999999_r8, 72.061400000000006_r8, 60.050400000000003_r8, &
    76.049800000000005_r8, 32.039999999999999_r8, 48.039400000000001_r8, 16.040600000000001_r8, &
    252.73040000000000_r8, 35.452700000000000_r8, 70.905400000000000_r8, 102.90420000000000_r8, &
    51.452100000000002_r8, 97.457639999999998_r8, 100.91685000000000_r8, 28.010400000000001_r8, &
    44.009799999999998_r8, 66.007205999999996_r8, 82.461502999999993_r8, 108.13560000000000_r8, &
    62.132399999999997_r8, 28.010400000000001_r8, 78.064599999999999_r8, 18.998403000000000_r8, &
    60.050400000000003_r8, 58.035600000000002_r8, 1.0074000000000001_r8, 2.0148000000000001_r8, &
    259.82361300000002_r8, 34.013599999999997_r8, 98.078400000000002_r8, 80.911400000000000_r8, &
    116.94800300000000_r8, 100.49370600000000_r8, 86.467905999999999_r8, 36.460099999999997_r8, &
    27.025140000000000_r8, 46.024600000000000_r8, 20.005803000000000_r8, 63.012340000000002_r8, &
    79.011740000000003_r8, 96.910799999999995_r8, 52.459499999999998_r8, 135.11493999999999_r8, &
    116.11239999999999_r8, 74.076200000000000_r8, 100.11300000000000_r8, 118.12720000000000_r8, &
    68.114199999999997_r8, 147.12594000000001_r8, 147.12594000000001_r8, 162.11794000000000_r8, &
    163.12533999999999_r8, 118.12720000000000_r8, 184.35020000000000_r8, 70.087800000000001_r8, &
    120.10080000000001_r8, 72.102599999999995_r8, 104.10140000000000_r8, 147.08474000000001_r8, &
    136.22839999999999_r8, 70.087800000000001_r8, 14.006740000000001_r8, 44.012880000000003_r8, &
    108.01048000000000_r8, 147.12594000000001_r8, 145.11114000000001_r8, 17.028939999999999_r8, &
    18.036339999999999_r8, 28.010400000000001_r8, 28.010400000000001_r8, 30.006139999999998_r8, &
    46.005540000000003_r8, 62.004939999999998_r8, 119.07434000000001_r8, 231.23954000000001_r8, &
    15.999400000000000_r8, 47.998199999999997_r8, 47.998199999999997_r8, 67.451499999999996_r8, &
    60.076400000000000_r8, 133.10014000000001_r8, 121.04794000000000_r8, 183.11774000000000_r8, &
    93.102400000000003_r8, 94.109800000000007_r8, 176.12160000000000_r8, 92.090400000000002_r8, &
    90.075599999999994_r8, 32.066000000000003_r8, 146.05641900000001_r8, 48.065399999999997_r8, &
    64.064800000000005_r8, 80.064200000000000_r8, 250.44499999999999_r8, 250.44499999999999_r8, &
    250.44499999999999_r8, 250.44499999999999_r8, 250.44499999999999_r8, 28.010400000000001_r8, &
    310.58240000000001_r8, 140.13440000000000_r8, 200.22600000000000_r8, 215.24014000000000_r8, &
    186.24140000000000_r8, 168.22720000000001_r8, 154.20140000000001_r8, 174.14800000000000_r8, &
    92.136200000000002_r8, 150.12600000000000_r8, 106.16200000000001_r8, 188.17380000000000_r8, &
    122.16140000000000_r8, 204.17320000000001_r8, 14.006740000000001_r8, 14.006740000000001_r8, &
    137.11220000000000_r8, 103.13520000000000_r8, 253.34819999999999_r8, 159.11480000000000_r8, &
    159.11480000000000_r8, 123.12760000000000_r8, 61.057800000000000_r8, 75.083600000000004_r8, &
    109.10180000000000_r8, 75.042400000000001_r8, 47.031999999999996_r8, 129.08959999999999_r8, &
    105.10880000000000_r8, 61.057800000000000_r8, 77.057199999999995_r8, 33.006200000000000_r8, &
    63.031399999999998_r8, 117.11980000000000_r8, 117.11980000000000_r8, 117.11980000000000_r8, &
    233.35579999999999_r8, 119.09340000000000_r8, 115.06380000000000_r8, 101.07920000000000_r8, &
    117.07859999999999_r8, 103.09399999999999_r8, 185.23400000000001_r8, 230.23213999999999_r8, &
    15.999400000000000_r8, 17.006799999999998_r8, 175.11420000000001_r8, 91.082999999999998_r8, &
    89.068200000000004_r8, 199.21860000000001_r8, 185.23400000000001_r8, 173.14060000000001_r8, &
    173.14060000000001_r8, 149.11859999999999_r8, 187.16640000000001_r8, 187.16640000000001_r8, &
    203.16579999999999_r8, 18.014199999999999_r8 /)
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
      nbins = 40
    end if
  end subroutine rad_cnst_get_info
  subroutine rad_cnst_get_info_by_bin(list_idx, bin_idx, nspec, bin_name)
    integer, intent(in) :: list_idx
    integer, intent(in) :: bin_idx
    integer, optional, intent(out) :: nspec
    character(len=*), optional, intent(out) :: bin_name
    if(present(bin_name)) then
      if (bin_idx > 20) then
        write(bin_name, '(A,I2.2)') 'PRSUL', bin_idx
      else
        write(bin_name, '(A,I2.2)') 'MXAER', bin_idx
      end if
    end if
    if(present(nspec)) then
      if (bin_idx > 20) then
        nspec = 1
      else
        nspec = 10
      end if
    end if
  end subroutine rad_cnst_get_info_by_bin
  subroutine rad_cnst_get_bin_props_by_idx(list_idx, bin_idx, spec_idx, spectype)
    integer, intent(in) :: list_idx
    integer, intent(in) :: bin_idx
    integer, intent(in) :: spec_idx
    character(len=*), optional, intent(out) :: spectype
    if(present(spectype)) then
      if (spec_idx == 1) then
        spectype = 'sulfate'
      else
        spectype = 'something else'
      end if
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
    if (trim(adjustl(short_name)) == 'PRSUL') then
      igroup = 1
    else if (trim(adjustl(short_name)) == 'MXAER') then
      igroup = 2
    else
      igroup = 0
    end if
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
    dryr(:,:) = 1.0e-7_r8 ! m
    rho(:,:) = 1.0e-11_r8 ! kg m-3
    rc = 0
  end subroutine carma_get_dry_radius
end module carma_intr


