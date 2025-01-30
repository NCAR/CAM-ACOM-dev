! Captures inputs and output to cloud aqueous chemistry routines
module cloud_aqueous_chemistry_snapshot

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: cloud_snapshot_init, cloud_snapshot_capture_input, cloud_snapshot_capture_output

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine cloud_snapshot_init()

    use cam_history, only : addfld, horiz_only

    integer :: i_elem
    character(len=10) :: index_string

    call addfld( 'cloud_press_in',       (/ 'lev' /), 'I', 'Pa',       'mid-point pressure' )
    call addfld( 'cloud_pdel_in',        (/ 'lev' /), 'I', 'Pa',       'pressure thickness of levels' )
    call addfld( 'cloud_tfld_in',        (/ 'lev' /), 'I', 'K',        'temperature' )
    call addfld( 'cloud_mbar_in',        (/ 'lev' /), 'I', 'AMU',      'mean wet atmopheric mass' )
    call addfld( 'cloud_lwc_in',         (/ 'lev' /), 'I', 'kg kg-1',  'cloud liquid water content' )
    call addfld( 'cloud_cldfrc_in',      (/ 'lev' /), 'I', 'unitless', 'cloud fraction' )
    call addfld( 'cloud_cldnum_in',      (/ 'lev' /), 'I', 'kg-1',     'droplet number concentration' )
    call addfld( 'cloud_xhnm_in',        (/ 'lev' /), 'I', 'cm-3',     'total atmospheric density' )
    do i_elem = 1, 7
      write(index_string, '(I10)') i_elem
      call addfld( 'cloud_invariants_'//trim(adjustl(index_string))//'_in', (/ 'lev' /), 'I', 'unknown', 'invariant' )
    end do
    do i_elem = 1, 26
      write(index_string, '(I10)') i_elem
      call addfld( 'cloud_qcw_'//trim(adjustl(index_string))//'_in', (/ 'lev' /), 'I', 'vmr', 'cloud-borne aerosol' )
      call addfld( 'cloud_qcw_'//trim(adjustl(index_string))//'_out', (/ 'lev' /), 'I', 'vmr', 'transported species' )
      call addfld( 'cloud_qin_'//trim(adjustl(index_string))//'_in', (/ 'lev' /), 'I', 'vmr', 'transported species' )
      call addfld( 'cloud_qin_'//trim(adjustl(index_string))//'_out', (/ 'lev' /), 'I', 'vmr', 'transported species' )
    end do
    call addfld( 'cloud_xphlwc_out',     (/ 'lev' /), 'I', 'unitless', 'ph value multiplied by cloud liquid water content' )
    call addfld( 'cloud_aqso4_out',      (/ 'lev' /), 'I', 'unknown',  'aqueous phase SO4 chemistry?' )
    call addfld( 'cloud_aqh2so4_out',    (/ 'lev' /), 'I', 'unknown',  'aqueous phase H2SO4 chemistry?' )
    call addfld( 'cloud_aqso4_h2o2_out', horiz_only,  'I', 'kg m-2',   'SO4 aqueous phase chemistry due to H2O2?' )
    call addfld( 'cloud_aqso4_o3_out',   horiz_only,  'I', 'kg m-2',   'SO4 aqueous phase chemistry due to O3?' )
    
  end subroutine cloud_snapshot_init

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine cloud_snapshot_capture_input(ncol, lchnk, loffset, dtime, press, &
    pdel, tfld, mbar, lwc, cldfrc, cldnum, xhnm, invariants, qcw, qin)

    use cam_history, only : outfld
    use cam_logfile, only : iulog
    use spmd_utils,  only : is_main_process => masterproc

    integer,          intent(in)    :: ncol              ! num of columns in chunk
    integer,          intent(in)    :: lchnk             ! chunk id
    integer,          intent(in)    :: loffset           ! offset of chem tracers in the advected tracers array
    real(r8),         intent(in)    :: dtime             ! time step (sec)
    real(r8),         intent(in)    :: press(:,:)        ! midpoint pressure ( Pa )
    real(r8),         intent(in)    :: pdel(:,:)         ! pressure thickness of levels (Pa)
    real(r8),         intent(in)    :: tfld(:,:)         ! temperature
    real(r8),         intent(in)    :: mbar(:,:)         ! mean wet atmospheric mass ( amu )
    real(r8), target, intent(in)    :: lwc(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), target, intent(in)    :: cldfrc(:,:)       ! cloud fraction
    real(r8),         intent(in)    :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8),         intent(in)    :: xhnm(:,:)         ! total atms density ( /cm**3)
    real(r8),         intent(in)    :: invariants(:,:,:)
    real(r8), target, intent(inout) :: qcw(:,:,:)        ! cloud-borne aerosol (vmr)
    real(r8),         intent(inout) :: qin(:,:,:)        ! transported species ( vmr )

    integer :: i_elem
    character(len=10) :: index_string

    if (is_main_process) then
      write(iulog,*) "*****************************"
      write(iulog,*) "Cloud Chemistry scalar inputs"
      write(iulog,*) "*****************************"
      write(iulog,*) "ncol: ", ncol
      write(iulog,*) "lchnk: ", lchnk
      write(iulog,*) "loffset: ", loffset
      write(iulog,*) "dtime: ", dtime
      write(iulog,*) "press dims: ", size(press, dim=1), size(press, dim=2)
      write(iulog,*) "invariants dims: ", size(invariants, dim=1), size(invariants, dim=2), size(invariants, dim=3)
      write(iulog,*) "qcw dims: ", size(qcw, dim=1), size(qcw, dim=2), size(qcw, dim=3)
      write(iulog,*) "qin dims: ", size(qin, dim=1), size(qin, dim=2), size(qin, dim=3)
    end if
    call outfld( 'cloud_press_in',  press,  ncol, lchnk )
    call outfld( 'cloud_pdel_in',   pdel,   ncol, lchnk )
    call outfld( 'cloud_tfld_in',   tfld,   ncol, lchnk )
    call outfld( 'cloud_mbar_in',   mbar,   ncol, lchnk )
    call outfld( 'cloud_lwc_in',    lwc,    ncol, lchnk )
    call outfld( 'cloud_cldfrc_in', cldfrc, ncol, lchnk )
    call outfld( 'cloud_cldnum_in', cldnum, ncol, lchnk )
    call outfld( 'cloud_xhnm_in',   xhnm,   ncol, lchnk )
    do i_elem = 1, 7
      write(index_string, '(I10)') i_elem
      call outfld( 'cloud_invariants_'//trim(adjustl(index_string))//'_in', invariants(:,:,i_elem), ncol, lchnk )
    end do
    do i_elem = 1, 26
      write(index_string, '(I10)') i_elem
      call outfld( 'cloud_qcw_'//trim(adjustl(index_string))//'_in', qcw(:,:,i_elem), ncol, lchnk )
      call outfld( 'cloud_qin_'//trim(adjustl(index_string))//'_in', qin(:,:,i_elem), ncol, lchnk )
    end do

  end subroutine cloud_snapshot_capture_input

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine cloud_snapshot_capture_output(ncol, lchnk, qcw, qin, xphlwc, &
    aqso4, aqh2so4, aqso4_h2o2, aqso4_o3)

    use cam_history, only : outfld

    integer,          intent(in)    :: ncol              ! num of columns in chunk
    integer,          intent(in)    :: lchnk             ! chunk id
    real(r8), target, intent(inout) :: qcw(:,:,:)        ! cloud-borne aerosol (vmr)
    real(r8),         intent(inout) :: qin(:,:,:)        ! transported species ( vmr )
    real(r8),         intent(out)   :: xphlwc(:,:)       ! pH value multiplied by lwc

    real(r8),         intent(out)   :: aqso4(:,:)        ! aqueous phase chemistry
    real(r8),         intent(out)   :: aqh2so4(:,:)      ! aqueous phase chemistry
    real(r8),         intent(out)   :: aqso4_h2o2(:)     ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8),         intent(out)   :: aqso4_o3(:)       ! SO4 aqueous phase chemistry due to O3 (kg/m2)

    integer :: i_elem
    character(len=10) :: index_string

    call outfld( 'cloud_xphlwc_out',     xphlwc,     ncol, lchnk )
    call outfld( 'cloud_aqso4_out',      aqso4,      ncol, lchnk )
    call outfld( 'cloud_aqh2so4_out',    aqh2so4,    ncol, lchnk )
    call outfld( 'cloud_aqso4_h2o2_out', aqso4_h2o2, ncol, lchnk )
    call outfld( 'cloud_aqso4_o3_out',   aqso4_o3,   ncol, lchnk )
    do i_elem = 1, 26
      write(index_string, '(I10)') i_elem
      call outfld( 'cloud_qcw_'//trim(adjustl(index_string))//'_out', qcw(:,:,i_elem), ncol, lchnk )
      call outfld( 'cloud_qin_'//trim(adjustl(index_string))//'_out', qin(:,:,i_elem), ncol, lchnk )
    end do

  end subroutine cloud_snapshot_capture_output

end module cloud_aqueous_chemistry_snapshot