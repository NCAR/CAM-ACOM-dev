module test_cloud_aqueous_chemistry_mod

  use shr_kind_mod, only: r8 => shr_kind_r8

! No invariants are needed for CARMA
#define NO_INVARIANTS

  implicit none

  integer, parameter :: n_ranks = 29
  integer, parameter :: n_timesteps = 1
  integer, parameter :: n_lat = 19
  integer, parameter :: n_lon = 24
  integer, parameter :: n_columns_per_rank = 16
  integer, parameter :: n_layers = 32
  integer, parameter :: n_invariants = 3
  integer, parameter :: n_species = 202
  integer, parameter :: n_modes = 40

#ifndef NO_INVARIANTS
  integer, parameter :: needed_invariant_indices(0) = (/  /)
#endif

  integer, parameter :: needed_species_indices(*) = (/ 74, 75, 84, 112, 122, 137, 176 /)
  integer, parameter :: extra_qcw_indices(*) = (/ 1, 11, 21, 31, 41, 51, 61, &
                                                 71, 81, 91, 101, 111, 121, 131, &
                                                 141, 151, 161, 171, 181, 191 /)

  character(len=*), parameter :: data_file = &
      '../test/chemistry/data/QPCARMATS-f10_f10_mg37-cloud-chemistry-mod.cam.h1i.0001-01-01-05400.nc'

  type :: chemistry_args
    integer :: ncol
    integer :: lchnk
    integer :: loffset
    real(r8) :: dtime
    real(r8) :: time
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: lon(:)
    real(r8), allocatable :: pres(:,:)
    real(r8), allocatable :: pdel(:,:)
    real(r8), allocatable :: tfld(:,:)
    real(r8), allocatable :: mbar(:,:)
    real(r8), allocatable :: lwc(:,:)
    real(r8), allocatable :: cldfrc(:,:)
    real(r8), allocatable :: cldnum(:,:)
    real(r8), allocatable :: xhnm(:,:)
    real(r8), allocatable :: invariants(:,:,:)
    real(r8), allocatable :: qcw(:,:,:)
    real(r8), allocatable :: qin(:,:,:)
    real(r8), allocatable :: xphlwc(:,:)
    real(r8), allocatable :: aqso4(:,:)
    real(r8), allocatable :: aqh2so4(:,:)
    real(r8), allocatable :: aqso4_h2o2(:)
    real(r8), allocatable :: aqso4_o3(:)
  end type chemistry_args

  interface read_and_condense_columns
    module procedure read_and_condense_columns_2D
    module procedure read_and_condense_columns_1D
  end interface read_and_condense_columns

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_as_primary_rank_against_snapshot()
   
    use spmd_utils, only: masterproc
    use physics_buffer, only: physics_buffer_desc
    use physics_types, only: physics_state
    use cloud_aqueous_chemistry, only: sox_inti => initialize, &
                                       setsox => calculate

    type(chemistry_args), allocatable :: chem_args(:), expected_outputs(:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state) :: state
    integer :: i

    allocate(pbuf(n_columns_per_rank))
    masterproc = .true.
    chem_args = get_inputs()
    expected_outputs = get_expected_outputs()
    call sox_inti()
    do i = 1, size(chem_args)
      call setsox( &
          state, &
          pbuf, &
          chem_args(i)%ncol, &
          chem_args(i)%lchnk, &
          chem_args(i)%loffset, &
          chem_args(i)%dtime, &
          chem_args(i)%pres, &
          chem_args(i)%pdel, &
          chem_args(i)%tfld, &
          chem_args(i)%mbar, &
          chem_args(i)%lwc, &
          chem_args(i)%cldfrc, &
          chem_args(i)%cldnum, &
          chem_args(i)%xhnm, &
          chem_args(i)%invariants, &
          chem_args(i)%qcw, &
          chem_args(i)%qin, &
          chem_args(i)%xphlwc, &
          chem_args(i)%aqso4, &
          chem_args(i)%aqh2so4, &
          chem_args(i)%aqso4_h2o2, &
          chem_args(i)%aqso4_o3 &
        )
    end do
    deallocate(pbuf)
    if (compare_outputs(chem_args, expected_outputs, &
                        1.0e-6_r8, 1.0e-40_r8)) then
      write(*,*) 'Primary rank test against snapshot passed'
    else
      write(*,*) 'Primary rank test against snapshot failed'
      stop 3
    end if

  end subroutine test_as_primary_rank_against_snapshot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_as_other_rank_against_snapshot()

    use spmd_utils, only: masterproc
    use physics_buffer, only: physics_buffer_desc
    use physics_types, only: physics_state
    use cloud_aqueous_chemistry, only: sox_inti => initialize, &
                                       setsox => calculate

    type(chemistry_args), allocatable :: chem_args(:), expected_outputs(:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state) :: state
    integer :: i

    allocate(pbuf(n_columns_per_rank))
    masterproc = .false.
    chem_args = get_inputs()
    expected_outputs = get_expected_outputs()
    call sox_inti()
    do i = 1, size(chem_args)
      call setsox( &
          state, &
          pbuf, &
          chem_args(i)%ncol, &
          chem_args(i)%lchnk, &
          chem_args(i)%loffset, &
          chem_args(i)%dtime, &
          chem_args(i)%pres, &
          chem_args(i)%pdel, &
          chem_args(i)%tfld, &
          chem_args(i)%mbar, &
          chem_args(i)%lwc, &
          chem_args(i)%cldfrc, &
          chem_args(i)%cldnum, &
          chem_args(i)%xhnm, &
          chem_args(i)%invariants, &
          chem_args(i)%qcw, &
          chem_args(i)%qin, &
          chem_args(i)%xphlwc, &
          chem_args(i)%aqso4, &
          chem_args(i)%aqh2so4, &
          chem_args(i)%aqso4_h2o2, &
          chem_args(i)%aqso4_o3 &
        )
      end do
    deallocate(pbuf)
    if (compare_outputs(chem_args, expected_outputs, &
                        1.0e-6_r8, 1.0e-40_r8)) then
      write(*,*) 'Other rank test against snapshot passed'
    else
      write(*,*) 'Other rank test against snapshot failed'
      stop 3
    end if

  end subroutine test_as_other_rank_against_snapshot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_as_primary_against_original_module()
  
    use spmd_utils, only: masterproc
    use physics_buffer, only: physics_buffer_desc
    use physics_types, only: physics_state
    use cloud_aqueous_chemistry, only: new_sox_inti => initialize, &
                                       new_setsox => calculate
    use mo_setsox, only: old_sox_inti => sox_inti, old_setsox => setsox

    type(chemistry_args), allocatable :: new_args(:), old_args(:)
    type(chemistry_args), allocatable :: expected_outputs(:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state) :: state
    integer :: i

    allocate(pbuf(n_columns_per_rank))
    new_args = get_inputs()
    old_args = get_inputs()
    expected_outputs = get_expected_outputs()
    call new_sox_inti()
    call old_sox_inti()
    do i = 1, size(new_args)
      call new_setsox( &
          state, &
          pbuf, &
          new_args(i)%ncol, &
          new_args(i)%lchnk, &
          new_args(i)%loffset, &
          new_args(i)%dtime, &
          new_args(i)%pres, &
          new_args(i)%pdel, &
          new_args(i)%tfld, &
          new_args(i)%mbar, &
          new_args(i)%lwc, &
          new_args(i)%cldfrc, &
          new_args(i)%cldnum, &
          new_args(i)%xhnm, &
          new_args(i)%invariants, &
          new_args(i)%qcw, &
          new_args(i)%qin, &
          new_args(i)%xphlwc, &
          new_args(i)%aqso4, &
          new_args(i)%aqh2so4, &
          new_args(i)%aqso4_h2o2, &
          new_args(i)%aqso4_o3 &
        )
    end do
    do i = 1, size(old_args)
      call old_setsox( &
          state, &
          pbuf, &
          old_args(i)%ncol, &
          old_args(i)%lchnk, &
          old_args(i)%loffset, &
          old_args(i)%dtime, &
          old_args(i)%pres, &
          old_args(i)%pdel, &
          old_args(i)%tfld, &
          old_args(i)%mbar, &
          old_args(i)%lwc, &
          old_args(i)%cldfrc, &
          old_args(i)%cldnum, &
          old_args(i)%xhnm, &
          old_args(i)%invariants, &
          old_args(i)%qcw, &
          old_args(i)%qin, &
          old_args(i)%xphlwc, &
          old_args(i)%aqso4, &
          old_args(i)%aqh2so4, &
          old_args(i)%aqso4_h2o2, &
          old_args(i)%aqso4_o3 &
        )
    end do
    deallocate(pbuf)
    if (compare_outputs(new_args, old_args)) then
      write(*,*) 'Primary rank test against original module passed'
    else
      write(*,*) 'Primary rank test against original module failed'
      stop 3
    end if

  end subroutine test_as_primary_against_original_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_as_other_rank_against_original_module()
  
    use spmd_utils, only: masterproc
    use physics_buffer, only: physics_buffer_desc
    use physics_types, only: physics_state
    use cloud_aqueous_chemistry, only: new_sox_inti => initialize, &
                                       new_setsox => calculate
    use mo_setsox, only: old_sox_inti => sox_inti, old_setsox => setsox

    type(chemistry_args), allocatable :: new_args(:), old_args(:)
    type(chemistry_args), allocatable :: expected_outputs(:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state) :: state
    integer :: i

    allocate(pbuf(n_columns_per_rank))
    masterproc = .false.
    new_args = get_inputs()
    old_args = get_inputs()
    expected_outputs = get_expected_outputs()
    call new_sox_inti()
    call old_sox_inti()
    do i = 1, size(new_args)
      call new_setsox( &
          state, &
          pbuf, &
          new_args(i)%ncol, &
          new_args(i)%lchnk, &
          new_args(i)%loffset, &
          new_args(i)%dtime, &
          new_args(i)%pres, &
          new_args(i)%pdel, &
          new_args(i)%tfld, &
          new_args(i)%mbar, &
          new_args(i)%lwc, &
          new_args(i)%cldfrc, &
          new_args(i)%cldnum, &
          new_args(i)%xhnm, &
          new_args(i)%invariants, &
          new_args(i)%qcw, &
          new_args(i)%qin, &
          new_args(i)%xphlwc, &
          new_args(i)%aqso4, &
          new_args(i)%aqh2so4, &
          new_args(i)%aqso4_h2o2, &
          new_args(i)%aqso4_o3 &
        )
    end do
    do i = 1, size(old_args)
      call old_setsox( &
          state, &
          pbuf, &
          old_args(i)%ncol, &
          old_args(i)%lchnk, &
          old_args(i)%loffset, &
          old_args(i)%dtime, &
          old_args(i)%pres, &
          old_args(i)%pdel, &
          old_args(i)%tfld, &
          old_args(i)%mbar, &
          old_args(i)%lwc, &
          old_args(i)%cldfrc, &
          old_args(i)%cldnum, &
          old_args(i)%xhnm, &
          old_args(i)%invariants, &
          old_args(i)%qcw, &
          old_args(i)%qin, &
          old_args(i)%xphlwc, &
          old_args(i)%aqso4, &
          old_args(i)%aqh2so4, &
          old_args(i)%aqso4_h2o2, &
          old_args(i)%aqso4_o3 &
        )
    end do
    deallocate(pbuf)
    if (compare_outputs(new_args, old_args)) then
      write(*,*) 'Other rank test against original module passed'
    else
      write(*,*) 'Other rank test against original module failed'
      stop 3
    end if

  end subroutine test_as_other_rank_against_original_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_inputs() result(chem_args)

    use file_io, only: file_io_t

    type(chemistry_args), allocatable :: chem_args(:)
    type(file_io_t), pointer :: file

    integer :: i, j, i_time, i_lat, i_lon
    character(len=10) :: index_string
    real(r8), allocatable :: temp1d(:), temp2d(:,:), lats(:), lons(:)

    call initialize_args(chem_args)
    file => file_io_t(data_file)

    allocate(temp1d(1))
    call file%read('time', temp1d)
    do i = 1, size(chem_args)
      chem_args(i)%time = temp1d(1)
    end do
    allocate(lats(n_lat))
    allocate(lons(n_lon))
    call file%read('lat', lats)
    call file%read('lon', lons)
    i = 1
    j = 0
    do i_time = 1, n_timesteps
      do i_lat = 1, n_lat
        do i_lon = 1, n_lon
          j = j+1
          if (j>n_columns_per_rank) then
            i = i + 1
            j = 1
          end if
          chem_args(i)%lat(j) = lats(i_lat)
          chem_args(i)%lon(j) = lons(i_lon)
        end do
      end do
    end do
    call read_and_condense_columns(file, 'cloud_press_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%pres = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_pdel_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%pdel = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_tfld_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%tfld = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_mbar_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%mbar = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_lwc_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%lwc = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_cldfrc_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%cldfrc = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_cldnum_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%cldnum = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    call read_and_condense_columns(file, 'cloud_xhnm_in', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%xhnm = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    ! No invariants in BAM configuration
    do i = 1, size(chem_args)
      chem_args(i)%invariants = 0.0_r8
    end do
#ifndef NO_INVARIANTS
    do j = 1, size(needed_invariant_indices)
      write(index_string, '(I3)') needed_invariant_indices(j)
      call read_and_condense_columns(file, 'cloud_invariants_'// &
              trim(adjustl(index_string))//'_in', temp2d, n_layers)
      do i = 1, size(chem_args)
        chem_args(i)%invariants(:,:,needed_invariant_indices(j)) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
      end do
    end do
#endif
    do i = 1, size(chem_args)
      chem_args(i)%qcw = 0.0_r8
      chem_args(i)%qin = 0.0_r8
    end do
    do j = 1, size(needed_species_indices)
      write(index_string, '(I3)') needed_species_indices(j)
      call read_and_condense_columns(file, &
          'cloud_qcw_'//trim(adjustl(index_string)) //'_in', temp2d, n_layers)
      do i = 1, size(chem_args)
        chem_args(i)%qcw(:,:,needed_species_indices(j)) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
      end do
      call read_and_condense_columns(file, &
          'cloud_qin_'//trim(adjustl(index_string)) //'_in', temp2d, n_layers)
      do i = 1, size(chem_args)
        chem_args(i)%qin(:,:,needed_species_indices(j)) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
      end do
    end do
    do i = 1, size(chem_args)
      chem_args(i)%xphlwc = -1.0e300_r8
      chem_args(i)%aqso4 = -1.0e300_r8
      chem_args(i)%aqh2so4 = -1.0e300_r8
      chem_args(i)%aqso4_h2o2 = -1.0e300_r8
      chem_args(i)%aqso4_o3 = -1.0e300_r8
    end do
    deallocate(file)

  end function get_inputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_expected_outputs() result(chem_args)

    use file_io, only: file_io_t

    type(chemistry_args), allocatable :: chem_args(:)

    type(file_io_t), pointer :: file
    integer :: i, j
    character(len=10) :: index_string
    real(r8), allocatable :: temp1d(:), temp2d(:,:)

    chem_args = get_inputs()
    file => file_io_t(data_file)
    do i = 1, size(chem_args)
      chem_args(i)%qcw = 0.0_r8
      chem_args(i)%qin = 0.0_r8
    end do
    do j = 1, size(needed_species_indices)
      write(index_string, '(I3)') needed_species_indices(j)
      call read_and_condense_columns(file, &
          'cloud_qcw_'//trim(adjustl(index_string))//'_out', temp2d, n_layers)
      do i = 1, size(chem_args)
        chem_args(i)%qcw(:,:,needed_species_indices(j)) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
      end do
      call read_and_condense_columns(file, &
          'cloud_qin_'//trim(adjustl(index_string))//'_out', temp2d, n_layers)
      do i = 1, size(chem_args)
        chem_args(i)%qin(:,:,needed_species_indices(j)) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
      end do
    end do
    do j = 1, size(extra_qcw_indices)
      write(index_string, '(I3)') extra_qcw_indices(j)
      call read_and_condense_columns(file, &
          'cloud_qcw_'//trim(adjustl(index_string))//'_out', temp2d, n_layers)
      do i = 1, size(chem_args)
        chem_args(i)%qcw(:,:,extra_qcw_indices(j)) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
      end do
    end do
    call read_and_condense_columns(file, 'cloud_xphlwc_out', temp2d, n_layers)
    do i = 1, size(chem_args)
      chem_args(i)%xphlwc(:,:) = temp2d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp2d, dim=1)),:)
    end do
    do j = 1, n_modes
      write(index_string, '(I2)') j
      call read_and_condense_columns(file, &
          'cloud_aqso4_'//trim(adjustl(index_string))//'_out', temp1d)
      do i = 1, size(chem_args)
        chem_args(i)%aqso4(:,j) = temp1d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp1d)))
      end do
      call read_and_condense_columns(file, &
          'cloud_aqh2so4_'//trim(adjustl(index_string))//'_out', temp1d)
      do i = 1, size(chem_args)
        chem_args(i)%aqh2so4(:,j) = temp1d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp1d)))
      end do
    end do
    call read_and_condense_columns(file, 'cloud_aqso4_h2o2_out', temp1d)
    do i = 1, size(chem_args)
      chem_args(i)%aqso4_h2o2(:) = temp1d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp1d)))
    end do
    call read_and_condense_columns(file, 'cloud_aqso4_o3_out', temp1d)
    do i = 1, size(chem_args)
      chem_args(i)%aqso4_o3(:) = temp1d((i-1)*n_columns_per_rank+1 &
                            : min(i*n_columns_per_rank,size(temp1d)))
    end do
    deallocate(file)

  end function get_expected_outputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_args(chem_args)

    type(chemistry_args), allocatable, intent(inout) :: chem_args(:)

    integer :: i

    if (allocated(chem_args)) deallocate(chem_args)
    allocate(chem_args(n_ranks))
    do i = 1, size(chem_args)
      chem_args(i)%ncol = min(n_columns_per_rank, n_lat*n_lon - &
                                (i-1)*n_columns_per_rank)
      chem_args(i)%lchnk = 37
      chem_args(i)%loffset = 9
      chem_args(i)%dtime = 1800.0_r8
      allocate(chem_args(i)%lat(chem_args(i)%ncol))
      allocate(chem_args(i)%lon(chem_args(i)%ncol))
      allocate(chem_args(i)%pres(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%pdel(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%tfld(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%mbar(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%lwc(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%cldfrc(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%cldnum(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%xhnm(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%invariants(chem_args(i)%ncol, n_layers, &
                                       n_invariants))
      ! CARMA has an extra 18 elements in the third dimension of qcw for some reason
      allocate(chem_args(i)%qcw(chem_args(i)%ncol, n_layers, n_species+18))
      allocate(chem_args(i)%qin(chem_args(i)%ncol, n_layers, n_species))
      allocate(chem_args(i)%xphlwc(chem_args(i)%ncol, n_layers))
      allocate(chem_args(i)%aqso4(chem_args(i)%ncol, n_modes))
      allocate(chem_args(i)%aqh2so4(chem_args(i)%ncol, n_modes))
      allocate(chem_args(i)%aqso4_h2o2(chem_args(i)%ncol))
      allocate(chem_args(i)%aqso4_o3(chem_args(i)%ncol))
    end do

  end subroutine initialize_args
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_and_condense_columns_2D(file, variable_name, data, dim2)

    use file_io, only: file_io_t

    type(file_io_t), intent(in) :: file
    character(len=*), intent(in) :: variable_name
    real(r8), dimension(:,:), allocatable, intent(inout) :: data
    integer, intent(in) :: dim2

    integer :: i_time, i_lat, i_lon, i_col
    real(r8), allocatable :: temp(:,:,:,:)

    if (allocated(data)) deallocate(data)
    allocate(data(n_lat*n_lon*n_timesteps, dim2))
    allocate(temp(n_lon, n_lat, dim2, n_timesteps))
    call file%read(variable_name, temp)
    i_col = 1
    do i_time = 1, n_timesteps
      do i_lat = 1, n_lat
        do i_lon = 1, n_lon
          data(i_col,:) = temp(i_lon,i_lat,:,i_time)
          i_col = i_col + 1
        end do
      end do
    end do

  end subroutine read_and_condense_columns_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_and_condense_columns_1D(file, variable_name, data)

    use file_io, only: file_io_t

    type(file_io_t), intent(in) :: file
    character(len=*), intent(in) :: variable_name
    real(r8), dimension(:), allocatable, intent(inout) :: data

    integer :: i_time, i_lat, i_lon, i_col
    real(r8), allocatable :: temp(:,:,:)

    if (allocated(data)) deallocate(data)
    allocate(data(n_lat*n_lon*n_timesteps))
    allocate(temp(n_lon, n_lat, n_timesteps))
    call file%read(variable_name, temp)
    i_col = 1
    do i_time = 1, n_timesteps
      do i_lat = 1, n_lat
        do i_lon = 1, n_lon
          data(i_col) = temp(i_lon,i_lat,i_time)
          i_col = i_col + 1
        end do
      end do
    end do

  end subroutine read_and_condense_columns_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function compare_outputs(calculated, expected, relative_tolerance, &
      absolute_tolerance) result(passed)
  
    type(chemistry_args), intent(in) :: calculated(:)
    type(chemistry_args), intent(in) :: expected(:)
    real(r8), optional, intent(in) :: relative_tolerance
    real(r8), optional, intent(in) :: absolute_tolerance

    integer :: i, j, k, l
    passed = .true.

    do i = 1, size(calculated)
      do j = 1, calculated(i)%ncol
        do k = 1, n_layers
          if (.not. check_close(calculated(i)%pres(j,k), &
              expected(i)%pres(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'pres mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%pres(j,k), ' expected: ', &
              expected(i)%pres(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%pdel(j,k), &
              expected(i)%pdel(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'pdel mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%pdel(j,k), ' expected: ', &
              expected(i)%pdel(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%tfld(j,k), &
              expected(i)%tfld(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'tfld mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%tfld(j,k), ' expected: ', &
              expected(i)%tfld(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%mbar(j,k), &
              expected(i)%mbar(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'mbar mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%mbar(j,k), ' expected: ', &
              expected(i)%mbar(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%lwc(j,k), &
              expected(i)%lwc(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'lwc mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%lwc(j,k), ' expected: ', &
              expected(i)%lwc(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%cldfrc(j,k), &
              expected(i)%cldfrc(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'cldfrc mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%cldfrc(j,k), ' expected: ', &
              expected(i)%cldfrc(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%cldnum(j,k), &
              expected(i)%cldnum(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'cldnum mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%cldnum(j,k), ' expected: ', &
              expected(i)%cldnum(j,k)
            passed = .false.
          end if
          if (.not. check_close(calculated(i)%xhnm(j,k), &
              expected(i)%xhnm(j,k), relative_tolerance, absolute_tolerance)) then
            print *, 'xhnm mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%xhnm(j,k), ' expected: ', &
              expected(i)%xhnm(j,k)
            passed = .false.
          end if
          do l = 1, n_invariants
            if (.not. check_close(calculated(i)%invariants(j,k,l), &
                expected(i)%invariants(j,k,l), relative_tolerance, &
                absolute_tolerance)) then
              print *, 'invariants mismatch at column ', (i-1)*j, ' and layer ', &
                k, ' and invariant ', l, ' calculated: ', &
                calculated(i)%invariants(j,k,l), ' expected: ', &
                expected(i)%invariants(j,k,l)
              passed = .false.
            end if
          end do
          do l = 1, n_species
            if (.not. check_close(calculated(i)%qcw(j,k,l), &
                expected(i)%qcw(j,k,l), relative_tolerance, &
                absolute_tolerance)) then
              print *, 'qcw mismatch at column ', (i-1)*j, ' and layer ', k, &
                ' and species ', l, ' calculated: ', calculated(i)%qcw(j,k,l), &
                ' expected: ', expected(i)%qcw(j,k,l)
              passed = .false.
            end if
            if (.not. check_close(calculated(i)%qin(j,k,l), &
                expected(i)%qin(j,k,l), relative_tolerance, &
                absolute_tolerance)) then
              print *, 'qin mismatch at column ', (i-1)*j, ' and layer ', k, &
                ' and species ', l, ' calculated: ', calculated(i)%qin(j,k,l), &
                ' expected: ', expected(i)%qin(j,k,l)
              passed = .false.
            end if
          end do
          if (.not. check_close(calculated(i)%xphlwc(j,k), &
              expected(i)%xphlwc(j,k), relative_tolerance, &
              absolute_tolerance)) then
            print *, 'xphlwc mismatch at column ', (i-1)*j, ' and layer ', k, &
              ' calculated: ', calculated(i)%xphlwc(j,k), ' expected: ', &
              expected(i)%xphlwc(j,k)
            passed = .false.
          end if
          do l = 1, n_modes-1 ! the last mode is not actually set in the origical code
            if (.not. check_close(calculated(i)%aqso4(j,l), &
                expected(i)%aqso4(j,l), relative_tolerance, &
                absolute_tolerance)) then
              print *, 'aqso4 mismatch at column ', (i-1)*j, ' and mode ', l, &
                ' calculated: ', calculated(i)%aqso4(j,l), ' expected: ', &
                expected(i)%aqso4(j,l)
              passed = .false.
            end if
            if (.not. check_close(calculated(i)%aqh2so4(j,l), &
                expected(i)%aqh2so4(j,l), relative_tolerance, &
                absolute_tolerance)) then
              print *, 'aqh2so4 mismatch at column ', (i-1)*j, ' and mode ', &
                l, ' calculated: ', calculated(i)%aqh2so4(j,l), ' expected: ',&
                expected(i)%aqh2so4(j,l)
              passed = .false.
            end if
          end do
        end do
        if (.not. check_close(calculated(i)%aqso4_h2o2(j), &
            expected(i)%aqso4_h2o2(j), relative_tolerance, &
            absolute_tolerance)) then
          print *, 'aqso4_h2o2 mismatch at column ', (i-1)*j, &
            ' calculated: ', calculated(i)%aqso4_h2o2(j), ' expected: ', &
            expected(i)%aqso4_h2o2(j)
          passed = .false.
        end if
        if (.not. check_close(calculated(i)%aqso4_o3(j), &
            expected(i)%aqso4_o3(j), relative_tolerance, &
            absolute_tolerance)) then
          print *, 'aqso4_o3 mismatch at column ', (i-1)*j, &
            ' calculated: ', calculated(i)%aqso4_o3(j), ' expected: ', &
            expected(i)%aqso4_o3(j)
          passed = .false.
        end if
      end do
    end do

  end function compare_outputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function check_close(a, b, relative_tolerance, absolute_tolerance)

    real(r8), intent(in) :: a, b
    real(r8), optional, intent(in) :: relative_tolerance
    real(r8), optional, intent(in) :: absolute_tolerance

    real(r8) :: l_rel_tol = 1.0e-13_r8
    real(r8) :: l_abs_tol = 0.0_r8

    if (present(relative_tolerance)) l_rel_tol = relative_tolerance
    if (present(absolute_tolerance)) l_abs_tol = absolute_tolerance
    check_close = abs(a - b) <= &
                  (abs(a) + abs(b))/2.0_r8 * l_rel_tol + l_abs_tol

  end function check_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module test_cloud_aqueous_chemistry_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_cloud_aqueous_chemistry

  use test_cloud_aqueous_chemistry_mod, only: &
      test_as_primary_rank_against_snapshot, &
      test_as_other_rank_against_snapshot, &
      test_as_primary_against_original_module, &
      test_as_other_rank_against_original_module

  !call test_as_primary_rank_against_snapshot()
  !call test_as_other_rank_against_snapshot()
  call test_as_primary_against_original_module()
  call test_as_other_rank_against_original_module()

end program test_cloud_aqueous_chemistry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!