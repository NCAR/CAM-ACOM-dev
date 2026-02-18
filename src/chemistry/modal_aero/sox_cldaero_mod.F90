!----------------------------------------------------------------------------------
! Modal aerosol implementation
!----------------------------------------------------------------------------------
module sox_cldaero_mod

  use shr_kind_mod,    only : r8 => shr_kind_r8
  use cam_abortutils,  only : endrun
  use ppgrid,          only : pcols, pver
  use mo_chem_utls,    only : get_spc_ndx
  use cldaero_mod,     only : cldaero_conc_t, cldaero_allocate, cldaero_deallocate
  use physconst,       only : gravit
  use phys_control,    only : phys_getopts, cam_chempkg_is
  use cldaero_mod,     only : cldaero_uptakerate
  use aerosol_properties_mod, only: aerosol_properties

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_destroy_obj

  integer :: id_msa, id_h2so4, id_so2, id_h2o2, id_nh3

  real(r8), parameter :: small_value = 1.e-20_r8

  integer :: ncnst_tot                  ! total number of mode number conc + mode species

  class(aerosol_properties), pointer :: aero_props => null()

  ! Module-level helper variables
  integer :: modeptr_accum_local = -1  ! Index of accumulation mode
  real(r8) :: specmw_so4 = -1     ! SO4 molecular weight (will be queried from properties)

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

  subroutine sox_cldaero_init(aero_props_in)

    class(aerosol_properties), target, intent(in) :: aero_props_in

    integer :: l, m
    logical :: history_aerosol      ! Output the MAM aerosol tendencies
    character(len=32) :: spectype

    id_msa = get_spc_ndx( 'MSA' )
    id_h2so4 = get_spc_ndx( 'H2SO4' )
    id_so2 = get_spc_ndx( 'SO2' )
    id_h2o2 = get_spc_ndx( 'H2O2' )
    id_nh3 = get_spc_ndx( 'NH3' )

    if (id_h2so4<1 .or. id_so2<1 .or. id_h2o2<1) then
      call endrun('sox_cldaero_init:MAM mech does not include necessary species' &
                  //' -- should not invoke sox_cldaero_mod ')
    endif

    call phys_getopts( history_aerosol_out        = history_aerosol   )
    !
    !   add to history
    !

    aero_props => aero_props_in

    ncnst_tot = aero_props%ncnst_tot()

    ! Find accumulation mode - look for first soluble mode with a sulfate species
    do m = 1, aero_props%nbins()
       if (aero_props%soluble(m)) then
          do l = 1, aero_props%nspecies(m)
             call aero_props%get(m, l, spectype=spectype)
             if (trim(spectype) == 'sulfate') then
                modeptr_accum_local = m
                exit
             endif
          enddo
       endif
       if (modeptr_accum_local > 0) exit
    enddo

    if (modeptr_accum_local < 0) modeptr_accum_local = 1  ! default to first mode

    specmw_so4 = 96.0_r8  ! g/mol for SO4

  end subroutine sox_cldaero_init

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset) result( conc_obj )

    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: qcw(:,:,:)
    real(r8), intent(in) :: lwc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: loffset

    type(cldaero_conc_t), pointer :: conc_obj

    integer :: l, m, mm
    integer :: i, k
    character(len=32) :: spectype
    logical :: has_ammonium

    conc_obj => cldaero_allocate()

    do k = 1,pver
       do i = 1,ncol
          if( cldfrc(i,k) >0._r8) then
             conc_obj%xlwc(i,k) = lwc(i,k) *cfact(i,k) ! cloud water L(water)/L(air)
             conc_obj%xlwc(i,k) = conc_obj%xlwc(i,k) / cldfrc(i,k) ! liquid water in the cloudy fraction of cell
          else
             conc_obj%xlwc(i,k) = 0._r8
          endif
       enddo
    enddo

    conc_obj%no3c(:,:) = 0._r8
    conc_obj%nh4c(:,:) = 0._r8
    conc_obj%so4c(:,:) = 0._r8

    ! Determine whether ammonium is tracked as a separate aerosol species.
    ! When it is not (e.g., MAM3), SO4 is treated as NH4HSO4 (charge -1),
    ! so so4_fact = 1.  The default from cldaero_allocate is 2 (SO4^2-, charge -2).
    has_ammonium = .false.
    do m = 1, aero_props%nbins()
       do l = 1, aero_props%nspecies(m)
          call aero_props%get(m, l, spectype=spectype)
          if (trim(spectype) == 'ammonium') then
             has_ammonium = .true.
             exit
          endif
       enddo
       if (has_ammonium) exit
    enddo
    if (.not. has_ammonium) conc_obj%so4_fact = 1._r8

    do k = 1,pver
       do i = 1,ncol
          do m = 1, aero_props%nbins()
            do l = 1, aero_props%nspecies(m)
               mm = aero_props%indexer(m,l)
               call  aero_props%get(m,l, spectype=spectype)
               if (trim(spectype) == 'sulfate') then
                  conc_obj%so4c(i,k) = conc_obj%so4c(i,k) + qcw(i,k,mm)
               end if
               if (trim(spectype) == 'ammonium') then
                  conc_obj%nh4c(i,k) = conc_obj%nh4c(i,k) + qcw(i,k,mm)
               end if
            end do
          end do
       end do
    end do

  end function sox_cldaero_create_obj

!----------------------------------------------------------------------------------
! Update the mixing ratios
!----------------------------------------------------------------------------------
  subroutine sox_cldaero_update( &
       state, ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
       delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin, &
       aqso4, aqh2so4, aqso4_h2o2, aqso4_o3, aqso4_h2o2_3d, aqso4_o3_3d)

    use physics_types, only: physics_state

    ! args

    type(physics_state), intent(in) :: state     ! Physics state variables

    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk ! chunk id
    integer,  intent(in) :: loffset

    real(r8), intent(in) :: dtime ! time step (sec)

    real(r8), intent(in) :: mbar(:,:) ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: pdel(:,:)
    real(r8), intent(in) :: press(:,:)
    real(r8), intent(in) :: tfld(:,:)

    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    real(r8), intent(in) :: xlwc(:,:)

    real(r8), intent(in) :: delso4_hprxn(:,:)
    real(r8), intent(in) :: xh2so4(:,:)
    real(r8), intent(in) :: xso4(:,:)
    real(r8), intent(in) :: xso4_init(:,:)
    real(r8), intent(in) :: nh3g(:,:)
    real(r8), intent(in) :: hno3g(:,:)
    real(r8), intent(in) :: xnh3(:,:)
    real(r8), intent(in) :: xhno3(:,:)
    real(r8), intent(in) :: xnh4c(:,:)
    real(r8), intent(in) :: xmsa(:,:)
    real(r8), intent(in) :: xso2(:,:)
    real(r8), intent(in) :: xh2o2(:,:)
    real(r8), intent(in) :: xno3c(:,:)

    real(r8), intent(inout) :: qcw(:,:,:) ! cloud-borne aerosol (vmr)
    real(r8), intent(inout) :: qin(:,:,:) ! xported species ( vmr )

    real(r8), intent(out) :: aqso4(:,:)                   ! aqueous phase chemistry
    real(r8), intent(out) :: aqh2so4(:,:)                 ! aqueous phase chemistry
    real(r8), intent(out) :: aqso4_h2o2(:)                ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8), intent(out) :: aqso4_o3(:)                  ! SO4 aqueous phase chemistry due to O3 (kg/m2)
    real(r8), intent(out), optional :: aqso4_h2o2_3d(:,:) ! 3D SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8), intent(out), optional :: aqso4_o3_3d(:,:)   ! 3D SO4 aqueous phase chemistry due to O3 (kg/m2)


    ! local vars ...

    ! Tendency arrays are indexed by the aerosol element index (1..ncnst_tot)
    real(r8) :: dqdt_aqso4(ncol,pver,ncnst_tot), &
         dqdt_aqh2so4(ncol,pver,ncnst_tot), &
         dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)

    real(r8), allocatable :: faqgain_msa(:), faqgain_so4(:), qnum_c(:)

    real(r8) :: delso4_o3rxn, &
         dso4dt_aqrxn, dso4dt_hprxn, &
         dso4dt_gasuptk, dmsadt_gasuptk, &
         dmsadt_gasuptk_tomsa, dmsadt_gasuptk_toso4, &
         dqdt_aq, dqdt_wr, dqdt

    real(r8) :: fwetrem, sumf, uptkrate
    real(r8) :: delnh3, delnh4

    integer :: l, n, m, mm
    integer :: ntot_msa_c

    integer :: i,k
    real(r8) :: xl
    character(len=32) :: spectype

    allocate(faqgain_msa(aero_props%nbins()))
    allocate(faqgain_so4(aero_props%nbins()))
    allocate(qnum_c(aero_props%nbins()))

    ! make sure dqdt is zero initially, for budgets
    dqdt_aqso4(:,:,:) = 0.0_r8
    dqdt_aqh2so4(:,:,:) = 0.0_r8
    dqdt_aqhprxn(:,:) = 0.0_r8
    dqdt_aqo3rxn(:,:) = 0.0_r8

    ! Avoid double counting in-cloud sulfur oxidation when running with
    ! GEOS-Chem. If running with GEOS-Chem then sulfur oxidation
    ! is performed internally to GEOS-Chem. Here, we just return to the
    ! parent routine and thus we do not apply tendencies calculated by MAM.
    if ( cam_chempkg_is('geoschem_mam4') ) return

    lev_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          cloud: if (cldfrc(i,k) >= 1.0e-5_r8) then
             xl = xlwc(i,k) ! / cldfrc(i,k)

             IF (XL .ge. 1.e-8_r8) THEN !! WHEN CLOUD IS PRESENTED

                delso4_o3rxn = xso4(i,k) - xso4_init(i,k)

                if (id_nh3>0) then
                   delnh3 = nh3g(i,k) - xnh3(i,k)
                   delnh4 = - delnh3
                endif

                !-------------------------------------------------------------------------
                ! compute factors for partitioning aerosol mass gains among modes
                ! the factors are proportional to the activated particle MR for each
                ! mode, which is the MR of cloud drops "associated with" the mode
                ! thus we are assuming the cloud drop size is independent of the
                ! associated aerosol mode properties (i.e., drops associated with
                ! Aitken and coarse sea-salt particles are same size)
                !
                ! qnum_c(n) = activated particle number MR for mode n (these are just
                ! used for partitioning among modes, so don't need to divide by cldfrc)

                do m = 1, aero_props%nbins()
                   mm = aero_props%indexer(m,0)
                   qnum_c(m) = max( 0.0_r8, qcw(i,k,mm) )
                end do

                ! force qnum_c(n) to be positive for n=modeptr_accum or n=1
                n = modeptr_accum_local
                if (n <= 0) n = 1
                qnum_c(n) = max( 1.0e-10_r8, qnum_c(n) )

                ! faqgain_so4(n) = fraction of total so4_c gain going to mode n
                ! these are proportional to the activated particle MR for each mode
                sumf = 0.0_r8
                do n = 1, aero_props%nbins()
                   faqgain_so4(n) = 0.0_r8
                   do l = 1, aero_props%nspecies(n)
                      call aero_props%get(n, l, spectype=spectype)
                      if (trim(spectype) == 'sulfate') then
                         faqgain_so4(n) = qnum_c(n)
                         sumf = sumf + faqgain_so4(n)
                         exit
                      endif
                   enddo
                end do

                if (sumf > 0.0_r8) then
                   do n = 1, aero_props%nbins()
                      faqgain_so4(n) = faqgain_so4(n) / sumf
                   end do
                end if
                ! at this point (sumf <= 0.0) only when all the faqgain_so4 are zero

                ! faqgain_msa(n) = fraction of total msa_c gain going to mode n
                ntot_msa_c = 0
                sumf = 0.0_r8
                do n = 1, aero_props%nbins()
                   faqgain_msa(n) = 0.0_r8
                   do l = 1, aero_props%nspecies(n)
                      call aero_props%get(n, l, spectype=spectype)
                      if (trim(spectype) == 'msa') then
                         faqgain_msa(n) = qnum_c(n)
                         ntot_msa_c = ntot_msa_c + 1
                         exit
                      endif
                   enddo
                   sumf = sumf + faqgain_msa(n)
                end do

                if (sumf > 0.0_r8) then
                   do n = 1, aero_props%nbins()
                      faqgain_msa(n) = faqgain_msa(n) / sumf
                   end do
                end if
                ! at this point (sumf <= 0.0) only when all the faqgain_msa are zero

                uptkrate = cldaero_uptakerate( xl, cldnum(i,k), cfact(i,k), cldfrc(i,k), tfld(i,k), press(i,k) )
                ! average uptake rate over dtime
                uptkrate = (1.0_r8 - exp(-min(100._r8,dtime*uptkrate))) / dtime

                ! dso4dt_gasuptk = so4_c tendency from h2so4 gas uptake (mol/mol/s)
                ! dmsadt_gasuptk = msa_c tendency from msa gas uptake (mol/mol/s)
                dso4dt_gasuptk = xh2so4(i,k) * uptkrate
                if (id_msa > 0) then
                   dmsadt_gasuptk = xmsa(i,k) * uptkrate
                else
                   dmsadt_gasuptk = 0.0_r8
                end if

                ! if no modes have msa aerosol, then "rename" scavenged msa gas to so4
                dmsadt_gasuptk_toso4 = 0.0_r8
                dmsadt_gasuptk_tomsa = dmsadt_gasuptk
                if (ntot_msa_c == 0) then
                   dmsadt_gasuptk_tomsa = 0.0_r8
                   dmsadt_gasuptk_toso4 = dmsadt_gasuptk
                end if

                !-----------------------------------------------------------------------
                ! now compute TMR tendencies
                ! this includes the above aqueous so2 chemistry AND
                ! the uptake of highly soluble aerosol precursor gases (h2so4, msa, ...)
                ! AND the wetremoval of dissolved, unreacted so2 and h2o2

                dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn(i,k)) / dtime
                dso4dt_hprxn = delso4_hprxn(i,k) / dtime

                ! fwetrem = fraction of in-cloud-water material that is wet removed
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't have so4 & msa wet removal here

                ! compute TMR tendencies for so4 and msa aerosol-in-cloud-water
                do m = 1, aero_props%nbins()
                   do l = 1, aero_props%nspecies(m)
                      mm = aero_props%indexer(m,l)
                      call  aero_props%get(m,l, spectype=spectype)
                      if (trim(spectype) == 'sulfate') then

                         dqdt_aqso4(i,k,mm) = faqgain_so4(m)*dso4dt_aqrxn*cldfrc(i,k)

                         dqdt_aqh2so4(i,k,mm) = faqgain_so4(m)* &
                              (dso4dt_gasuptk + dmsadt_gasuptk_toso4)*cldfrc(i,k)
                         dqdt_aq = dqdt_aqso4(i,k,mm) + dqdt_aqh2so4(i,k,mm)
                         dqdt_wr = -fwetrem*dqdt_aq
                         dqdt= dqdt_aq + dqdt_wr
                         qcw(i,k,mm) = qcw(i,k,mm) + dqdt*dtime

                      end if
                      if (trim(spectype) == 'ammonium') then
                         if (delnh4 > 0.0_r8) then
                            dqdt_aq = faqgain_so4(m)*delnh4/dtime*cldfrc(i,k)
                            dqdt = dqdt_aq
                            qcw(i,k,mm) = qcw(i,k,mm) + dqdt*dtime
                         else
                            dqdt = (qcw(i,k,mm)/max(xnh4c(i,k),1.0e-35_r8)) &
                                 *delnh4/dtime*cldfrc(i,k)
                            qcw(i,k,mm) = qcw(i,k,mm) + dqdt*dtime
                         endif
                      end if
                   end do
                end do


                ! For gas species, tendency includes
                ! reactive uptake to cloud water that essentially transforms the gas to
                ! a different species. Wet removal associated with this is applied
                ! to the "new" species (e.g., so4_c) rather than to the gas.
                ! wet removal of the unreacted gas that is dissolved in cloud water.
                ! Need to multiply both these parts by cldfrc

                ! h2so4 (g) & msa (g)
                qin(i,k,id_h2so4) = qin(i,k,id_h2so4) - dso4dt_gasuptk * dtime * cldfrc(i,k)
                if (id_msa > 0) qin(i,k,id_msa) = qin(i,k,id_msa) - dmsadt_gasuptk * dtime * cldfrc(i,k)

                ! so2 -- the first order loss rate for so2 is frso2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frso2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't include so2 wet removal here

                dqdt_wr = -fwetrem*xso2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_aqrxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                qin(i,k,id_so2) = qin(i,k,id_so2) + dqdt * dtime

                ! h2o2 -- the first order loss rate for h2o2 is frh2o2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frh2o2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't include h2o2 wet removal here

                dqdt_wr = -fwetrem*xh2o2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_hprxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                qin(i,k,id_h2o2) = qin(i,k,id_h2o2) + dqdt * dtime

                ! NH3
                if (id_nh3>0) then
                   dqdt_aq = delnh3/dtime*cldfrc(i,k)
                   dqdt = dqdt_aq
                   qin(i,k,id_nh3) = qin(i,k,id_nh3) + dqdt * dtime
                endif

                ! for SO4 from H2O2/O3 budgets
                dqdt_aqhprxn(i,k) = dso4dt_hprxn*cldfrc(i,k)
                dqdt_aqo3rxn(i,k) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(i,k)

             ENDIF !! WHEN CLOUD IS PRESENTED
          endif cloud
       enddo col_loop
    enddo lev_loop

    !==============================================================
    ! ... Update the mixing ratios
    !==============================================================
    do k = 1,pver

       do n = 1, aero_props%nbins()
          do l = 1, aero_props%nspecies(n)
             mm = aero_props%indexer(n,l)
             call  aero_props%get(n,l, spectype=spectype)
             if (trim(spectype) == 'sulfate') then
                qcw(:,k,mm) = MAX(qcw(:,k,mm), small_value )
             end if
             if (trim(spectype) == 'msa') then
                qcw(:,k,mm) = MAX(qcw(:,k,mm), small_value )
             end if
             if (trim(spectype) == 'ammonium') then
                qcw(:,k,mm) = MAX(qcw(:,k,mm), small_value )
             end if
          end do
       end do

       qin(:,k,id_so2)   = MAX( qin(:,k,id_so2),   small_value )
       qin(:,k,id_h2o2)  = MAX( qin(:,k,id_h2o2),  small_value )
       qin(:,k,id_h2so4) = MAX( qin(:,k,id_h2so4), small_value )
       if ( id_msa > 0 ) qin(:,k,id_msa) = MAX( qin(:,k,id_msa), small_value )
       if ( id_nh3 > 0 ) qin(:,k,id_nh3) = MAX( qin(:,k,id_nh3), small_value )

    end do

    ! diagnostics

    ! Sum up diagnostics for each mode containing sulfate.
    ! Tendencies are in mol_SO4/mol_air/s; convert to kg_SO4/m2/s using
    ! specmw_so4 (g/mol) and mbar (g/mol) from the aerosol_properties interface.
    do n = 1, aero_props%nbins()
       aqso4(:,n) = 0._r8
       aqh2so4(:,n) = 0._r8

       do l = 1, aero_props%nspecies(n)
          call aero_props%get(n, l, spectype=spectype)
          if (trim(spectype) == 'sulfate') then
             mm = aero_props%indexer(n,l)
             do k=1,pver
                do i=1,ncol
                   aqso4(i,n)=aqso4(i,n)+dqdt_aqso4(i,k,mm)*specmw_so4/mbar(i,k) &
                        *pdel(i,k)/gravit ! kg/m2/s
                   aqh2so4(i,n)=aqh2so4(i,n)+dqdt_aqh2so4(i,k,mm)*specmw_so4/mbar(i,k) &
                        *pdel(i,k)/gravit ! kg/m2/s
                enddo
             enddo
             exit  ! Only process first sulfate species in mode
          endif
       enddo
    end do

    aqso4_h2o2(:) = 0._r8
    do k=1,pver
       do i=1,ncol
          aqso4_h2o2(i)=aqso4_h2o2(i)+dqdt_aqhprxn(i,k)*specmw_so4/mbar(i,k) &
                  *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo

    if (present(aqso4_h2o2_3d)) then
       aqso4_h2o2_3d(:,:) = 0._r8
       do k=1,pver
          do i=1,ncol
             aqso4_h2o2_3d(i,k)=dqdt_aqhprxn(i,k)*specmw_so4/mbar(i,k) &
                                *pdel(i,k)/gravit ! kg SO4 /m2/s
          enddo
       enddo
    end if

    aqso4_o3(:)=0._r8
    do k=1,pver
       do i=1,ncol
          aqso4_o3(i)=aqso4_o3(i)+dqdt_aqo3rxn(i,k)*specmw_so4/mbar(i,k) &
                  *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo

    if (present(aqso4_o3_3d)) then
       aqso4_o3_3d(:,:)=0._r8
       do k=1,pver
          do i=1,ncol
             aqso4_o3_3d(i,k)=dqdt_aqo3rxn(i,k)*specmw_so4/mbar(i,k) &
                              *pdel(i,k)/gravit ! kg SO4 /m2/s
          enddo
       enddo
    end if

    deallocate(faqgain_msa)
    deallocate(faqgain_so4)
    deallocate(qnum_c)

  end subroutine sox_cldaero_update

  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
  subroutine sox_cldaero_destroy_obj( conc_obj )
    type(cldaero_conc_t), pointer :: conc_obj

    call cldaero_deallocate( conc_obj )

  end subroutine sox_cldaero_destroy_obj

end module sox_cldaero_mod
