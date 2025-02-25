!----------------------------------------------------------------------------------
! Cloud aqueous chemistry
!
! The purpose of this module is to calculate the sulfate formation due to various
! oxidation/condensation pathways of cloud chemistry. These pathways include:
! - SO2 oxidation by O3
! - SO2 oxidation by H2O2
! - H2SO4 condensation
! Updated gas and aqueous species concentrations along with cloud pH changes are
! calculated.
!----------------------------------------------------------------------------------
module cloud_aqueous_chemistry

#define USE_MAM
#undef USE_CARMA

  use shr_kind_mod,      only : r8 => shr_kind_r8
  use cam_logfile,       only : iulog
  use physics_buffer,    only: physics_buffer_desc
  use physics_types,     only: physics_state

  implicit none

  private
  public :: initialize, calculate
  public :: do_cloud_aqueous_chemistry

  logical :: do_cloud_aqueous_chemistry = .false.

  integer, parameter :: CLOUD_INDEX_UNDEFINED = -1

  ! Cloud chemistry species information
  type :: cloud_species_t
    character(len=:), allocatable :: name_
    integer :: state_index_ = CLOUD_INDEX_UNDEFINED ! index in the state vector or the fixed concentrations array
    logical :: is_constant_ = .false.
  contains
    procedure :: exists => cloud_species_exists
    procedure :: mixing_ratio => cloud_species_get_mixing_ratio
  end type cloud_species_t

  interface cloud_species_t
    module procedure :: cloud_species_constructor
  end interface cloud_species_t

  type(cloud_species_t) :: so2, nh3, hno3, h2o2, o3, ho2, msa, so4, h2so4

  ! TODO: Figure out what this flag is for
  logical :: cloud_borne = .false.

  ! Constants that should be moved to a common module
  real(r8), parameter :: AVOGADRO = 6.02214129e23_r8 ! mol-1
  real(r8), parameter :: PASCAL_TO_ATM = 1.0_r8/101325.0_r8
  real(r8), parameter :: GAS_CONSTANT_L_ATM_MOL_K = 8314.46261815324_r8*PASCAL_TO_ATM
  real(r8), parameter :: BOLTZMANN = 1.380649e-23_r8 ! J K-1
  real(r8), parameter :: GAS_CONSTANT_DRY_AIR_J_KG_K = 287.0_r8 ! J kg-1 K-1

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine initialize()
    !-----------------------------------------------------------------------
    !	... prepare for cloud aqueous chemistry
    ! - Look up cloud chemistry species
    ! - Determine if enough species are present to perform cloud chemistry
    !-----------------------------------------------------------------------

    use spmd_utils,      only : masterproc
    use phys_control,    only : phys_getopts
    use carma_flags_mod, only : carma_do_cloudborne
#ifdef USE_MAM
    use mam_clouds,      only : sox_cldaero_init
#endif
#ifdef USE_CARMA
    use carma_clouds,    only : sox_cldaero_init
#endif

    logical :: is_modal_aerosols

    call phys_getopts( prog_modal_aero_out=is_modal_aerosols )
    cloud_borne = is_modal_aerosols .or. carma_do_cloudborne

    so2 = cloud_species_t( 'SO2' )
    nh3 = cloud_species_t( 'NH3' )
    hno3 = cloud_species_t( 'HNO3' )
    h2o2 = cloud_species_t( 'H2O2' )
    o3 = cloud_species_t( 'O3' )
    ho2 = cloud_species_t( 'HO2' )
    msa = cloud_species_t( 'MSA' )
    so4 = cloud_species_t( 'SO4' )
    h2so4 = cloud_species_t( 'H2SO4' )

    do_cloud_aqueous_chemistry = so2%exists() .and. h2o2%exists() .and. &
                                 o3%exists() .and. ho2%exists()
    if (do_cloud_aqueous_chemistry) then
      if (cloud_borne) then
        if (.not. h2so4%exists()) then
          do_cloud_aqueous_chemistry = .false.
        endif
      else
        if (.not. (so4%exists() .and. nh3%exists())) then
          do_cloud_aqueous_chemistry = .false.
        endif
      endif
    endif

    if (masterproc) then
      if( do_cloud_aqueous_chemistry ) then
         write(iulog,*) '-----------------------------------------'
         write(iulog,*) ' cloud aqueous chemistry is active'
         write(iulog,*) '-----------------------------------------'
      else
         write(iulog,*) '-----------------------------------------'
         write(iulog,*) ' cloud aqueous chemistry is inactive'
         write(iulog,*) '-----------------------------------------'
      end if
    end if
    if (.not. do_cloud_aqueous_chemistry) return
    
    call sox_cldaero_init()

  end subroutine initialize

!-----------------------------------------------------------------------
! Calculates the formation of sulfate and updates sulfate concentrations
! due to cloud aqueous chemistry. Also outputs the production rates
! (kg m-2 s-1) of various reactions of sulfur species with various
! oxidants (e.g., H2O2(aq), H2SO4(aq)). 
!-----------------------------------------------------------------------
  subroutine calculate( state,         &
       pbuf,                           &
       ncol,                           &
       lchnk,                          &
       loffset,                        &
       time_step,                      &
       midpoint_pressure,              &
       pressure_thickness,             &
       temperature,                    &
       mean_mass,                      &
       cloud_water,                    &
       cloud_fraction,                 &
       cloud_droplet_number,           &
       air_number_density,             &
       fixed_concentrations,           &
       cloud_borne_aerosol_vmr,        &
       species_vmr,                    &
       ph_times_cloud_water,           &
       aq_so4_production,              &
       aq_h2so4_production,            &
       aq_so4_production_from_h2o2,    &
       aq_so4_production_from_o3,      &
       specified_ph,                   &
       aq_so4_production_from_h2o2_3d, &
       aq_so4_production_from_o3_3d    &
       )

    !-----------------------------------------------------------------------
    !          ... Compute heterogeneous reactions of SOX
    !
    !       (0) using initial PH to calculate PH
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) PH values
    !
    !       (1) using new PH to repeat
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) REACTION rates
    !           (d) PREDICTION
    !-----------------------------------------------------------------------
    !
    use ppgrid,          only : pver
#ifdef USE_MAM
    use mam_clouds,      only : sox_cldaero_update, sox_cldaero_create_obj, &
                                sox_cldaero_destroy_obj
#endif
#ifdef USE_CARMA
    use carma_clouds,    only : sox_cldaero_update, sox_cldaero_create_obj, &
                                sox_cldaero_destroy_obj
#endif
    use cloud_utilities, only : cldaero_conc_t

    !
    !-----------------------------------------------------------------------
    !      ... Dummy arguments
    !-----------------------------------------------------------------------
    type(physics_state),                intent(in)    :: state         ! Physics state variables
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)       ! Physics buffer
    integer,          intent(in)    :: ncol                            ! num of columns in chunk
    integer,          intent(in)    :: lchnk                           ! chunk id
    integer,          intent(in)    :: loffset                         ! offset of chem tracers in the advected tracers array
    real(r8),         intent(in)    :: time_step                       ! time step (sec)
    real(r8),         intent(in)    :: midpoint_pressure(:,:)          ! midpoint pressure (Pa)
    real(r8),         intent(in)    :: pressure_thickness(:,:)         ! pressure thickness of levels (Pa) [pdel elsewhere in CAM]
    real(r8),         intent(in)    :: temperature(:,:)                ! temperature (K) [tfld elsewhere in CAM]
    real(r8),         intent(in)    :: mean_mass(:,:)                  ! mean wet atmospheric mass (amu)
    real(r8), target, intent(in)    :: cloud_water(:,:)                ! cloud liquid water content (kg/kg)
    real(r8), target, intent(in)    :: cloud_fraction(:,:)             ! cloud fraction (unitless)
    real(r8),         intent(in)    :: cloud_droplet_number(:,:)       ! droplet number concentration (#/kg)
    real(r8),         intent(in)    :: air_number_density(:,:)         ! total atmospheric number density (/cm**3)
    real(r8),         intent(in)    :: fixed_concentrations(:,:,:)     ! fixed concentrations (/cm**3) [invariants elsewhere in CAM]
    real(r8), target, intent(inout) :: cloud_borne_aerosol_vmr(:,:,:)  ! cloud-borne aerosol (vmr) mol/mol = m3/m3 [qcw elsewhere in CAM]
    real(r8),         intent(inout) :: species_vmr(:,:,:)              ! transported species (vmr) mol/mol = m3/m3 [qin elsewhere in CAM]
    real(r8),         intent(out)   :: ph_times_cloud_water(:,:)       ! pH value multiplied by cloud liquid water content

    real(r8),         intent(out)   :: aq_so4_production(:,:)         ! aqueous phase production of SO4 (kg/m2/s)
    real(r8),         intent(out)   :: aq_h2so4_production(:,:)       ! aqueous phase production of H2SO4 (kg/m2/s)
    real(r8),         intent(out)   :: aq_so4_production_from_h2o2(:) ! SO4 aqueous phase production due to H2O2 (kg/m2/s)
    real(r8),         intent(out)   :: aq_so4_production_from_o3(:)   ! SO4 aqueous phase production due to O3 (kg/m2/s)
    real(r8),         intent(in),  optional :: specified_ph                         ! specified pH value. If present, this value will be used instead of calculated pH
    real(r8),         intent(out), optional :: aq_so4_production_from_h2o2_3d(:, :) ! 3D SO4 aqueous phase production due to H2O2 (kg/m2/s)
    real(r8),         intent(out), optional :: aq_so4_production_from_o3_3d(:, :)   ! 3D SO4 aqueous phase production due to O3 (kg/m2/s)


    !-----------------------------------------------------------------------
    !      ... Local variables
    !-----------------------------------------------------------------------
    integer,  parameter :: MAX_ITERATIONS = 20
    real(r8), parameter :: INITIAL_PH = 5.0_r8  ! Initial pH value
    real(r8), parameter :: MINIMUM_CLOUD_LIQUID_WATER = 1.e-8_r8 ! Minimum cloud liquid water content (kg/kg)

    ! Effective Henry's Law constants for HO2 partitioning
    ! TODO: skipping remnaming of these in anticipation of a partitioning struct
    real(r8), parameter :: kh0 = 9.e3_r8            ! HO2(g)          -> Ho2(a)
    real(r8), parameter :: kh1 = 2.05e-5_r8         ! HO2(a)          -> H+ + O2-
    real(r8), parameter :: kh2 = 8.6e5_r8           ! HO2(a) + ho2(a) -> h2o2(a) + o2
    real(r8), parameter :: kh3 = 1.e8_r8            ! HO2(a) + o2-    -> h2o2(a) + o2

    ! Water dissociation constant (mol2/L2) [H+][OH-]
    real(r8), parameter :: water_dissociation_constant = 1.e-14_r8

    ! Change in aqueous sulfate volume mixing ratio over current time step (mol mol-1)
    real(r8) :: change_in_aq_so4_mixing_ratio(ncol,pver)

    ! TODO: Skipping renaming of these in anticipation of partitioning/pH estimation structs
    integer  :: k, i, iter
    real(r8) :: xk, xe, x2
    real(r8) :: xl, px, patm
    real(r8) :: Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
    real(r8) :: so2g, h2o2g, co2g, o3g
    real(r8) :: k_siv_h2o2  ! rate constant for reaction of S(IV) with H2O2
    real(r8) :: k_siv_o3    ! rate constant for reaction of S(IV) with O3
    real(r8) :: dso4_dt     ! rate of change of SO4
    real(r8) :: delta_concentration

    real(r8) :: hno3g(ncol,pver), nh3g(ncol,pver)
    !
    !-----------------------------------------------------------------------
    !            for Ho2(g) -> H2o2(a) formation
    !            schwartz JGR, 1984, 11589
    !-----------------------------------------------------------------------
    real(r8) :: ho2s   ! ho2s = ho2(a)+o2-
    real(r8) :: r1h2o2 ! prod(h2o2) by ho2 in mole/L(w)/s
    real(r8) :: r2h2o2 ! prod(h2o2) by ho2 in mix/s

    ! volume mixing ratios for cloud chemistry species
    real(r8), dimension(ncol,pver) ::  xhno3, xh2o2, xso2, xso4, xno3, &
         xnh3, xnh4, xo3, xph, xho2, xh2so4, xmsa, xso4_init

    real(r8), dimension(ncol,pver) :: air_mass_density_kg_l ! kg L-1
    
    ! Effective Henry's Law constants
    real(r8), dimension(ncol,pver) :: hehno3, heh2o2, heso2, henh3, heo3

    ! TODO: Figure out what this actually represents (it differs based on the value of cloud_borne)
    real(r8) :: patm_x

    real(r8), dimension(ncol)  :: work1
    logical :: converged

    ! Cloud-borne species volume mixing ratios
    real(r8), pointer :: xso4c(:,:)
    real(r8), pointer :: xnh4c(:,:)
    real(r8), pointer :: xno3c(:,:)
    type(cldaero_conc_t), pointer :: cldconc

    ! TODO: Skipping renaming of these in anticipation of partitioning/pH estimation structs
    real(r8) :: fact1_hno3, fact2_hno3, fact3_hno3
    real(r8) :: fact1_so2, fact2_so2, fact3_so2, fact4_so2
    real(r8) :: fact1_nh3, fact2_nh3, fact3_nh3
    real(r8) :: tmp_hp, tmp_hso3, tmp_hco3, tmp_nh4, tmp_no3
    real(r8) :: tmp_oh, tmp_so3, tmp_so4
    real(r8) :: tmp_neg, tmp_pos
    real(r8) :: yph, yph_lo, yph_hi
    real(r8) :: ynetpos, ynetpos_lo, ynetpos_hi

    !==================================================================
    !       ... First set the PH
    !==================================================================
    !      ... Initial values
    !           The values of so2, so4 are after (1) SLT, and CHEM
    !-----------------------------------------------------------------
    do k = 1,pver
       air_mass_density_kg_l(:,k) = &
              air_number_density(:,k)               & ! molecules(air) cm-3
            * 1.e3_r8                               & ! molecules(air) L-1
            * BOLTZMANN/GAS_CONSTANT_DRY_AIR_J_KG_K   ! kg(air) L-1
    end do

    cldconc => sox_cldaero_create_obj( cloud_fraction, cloud_borne_aerosol_vmr, &
                                       cloud_water, air_mass_density_kg_l, ncol, loffset )
    xso4c => cldconc%so4c
    xnh4c => cldconc%nh4c
    xno3c => cldconc%no3c
    xso4(:,:) = 0._r8
    xno3(:,:) = 0._r8
    xnh4(:,:) = 0._r8
    call so2%mixing_ratio(   species_vmr, fixed_concentrations, air_number_density, xso2   )
    call nh3%mixing_ratio(   species_vmr, fixed_concentrations, air_number_density, xnh3   )
    call hno3%mixing_ratio(  species_vmr, fixed_concentrations, air_number_density, xhno3  )
    call h2o2%mixing_ratio(  species_vmr, fixed_concentrations, air_number_density, xh2o2  )
    call o3%mixing_ratio(    species_vmr, fixed_concentrations, air_number_density, xo3    )
    call ho2%mixing_ratio(   species_vmr, fixed_concentrations, air_number_density, xho2   )
    call msa%mixing_ratio(   species_vmr, fixed_concentrations, air_number_density, xmsa   )
    call so4%mixing_ratio(   species_vmr, fixed_concentrations, air_number_density, xso4   )
    call h2so4%mixing_ratio( species_vmr, fixed_concentrations, air_number_density, xh2so4 )

    xph(:,:) = 10._r8**(-INITIAL_PH)                                ! initial PH value

    !-----------------------------------------------------------------
    !       ... Temperature dependent Henry constants
    !-----------------------------------------------------------------
    ver_loop0: do k = 1,pver                               !! pver loop for STEP 0
       col_loop0: do i = 1,ncol

          if (cloud_borne .and. cloud_fraction(i,k)>0._r8) then
             xso4(i,k) = xso4c(i,k) / cloud_fraction(i,k)
             xnh4(i,k) = xnh4c(i,k) / cloud_fraction(i,k)
             xno3(i,k) = xno3c(i,k) / cloud_fraction(i,k)
          endif
          xl = cldconc%xlwc(i,k)

          if( xl >= MINIMUM_CLOUD_LIQUID_WATER ) then
             work1(i) = 1._r8 / temperature(i,k) - 1._r8 / 298._r8

             !-----------------------------------------------------------------
             ! 21-mar-2011 changes by rce
             ! ph calculation now uses bisection method to solve the electro-neutrality equation
             ! 3-mode aerosols (where so4 is assumed to be nh4hso4)
             !    old code set xnh4c = so4c
             !    new code sets xnh4c = 0, then uses a -1 charge (instead of -2)
             !       for so4 when solving the electro-neutrality equation
             !-----------------------------------------------------------------

             !-----------------------------------------------------------------
             !  calculations done before iterating
             !-----------------------------------------------------------------

             !-----------------------------------------------------------------
             ! This should be divided by 101325, not 101300, but fixing this breaks the tests
             patm = midpoint_pressure(i,k)/101300._r8

             !-----------------------------------------------------------------
             !        ... hno3
             !-----------------------------------------------------------------
             ! previous code
             !    hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))
             !    px = hehno3(i,k) * GAS_CONSTANT_L_ATM_MOL_K * tz * xl
             !    hno3g = xhno3(i,k)/(1._r8 + px)
             !    Ehno3 = xk*xe*hno3g *patm
             ! equivalent new code
             !    hehno3 = xk + xk*xe/hplus
             !    hno3g = xhno3/(1 + px)
             !          = xhno3/(1 + hehno3*ra*tz*xl)
             !          = xhno3/(1 + xk*ra*tz*xl*(1 + xe/hplus)
             !    ehno3 = hno3g*xk*xe*patm
             !          = xk*xe*patm*xhno3/(1 + xk*ra*tz*xl*(1 + xe/hplus)
             !          = ( fact1_hno3    )/(1 + fact2_hno3 *(1 + fact3_hno3/hplus)
             !    [hno3-] = ehno3/hplus
             xk = 2.1e5_r8 *EXP( 8700._r8*work1(i) )
             xe = 15.4_r8
             fact1_hno3 = xk*xe*patm*xhno3(i,k)
             fact2_hno3 = xk*GAS_CONSTANT_L_ATM_MOL_K*temperature(i,k)*xl
             fact3_hno3 = xe

             !-----------------------------------------------------------------
             !          ... so2
             !-----------------------------------------------------------------
             ! previous code
             !    heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))
             !    px = heso2(i,k) * GAS_CONSTANT_L_ATM_MOL_K * tz * xl
             !    so2g =  xso2(i,k)/(1._r8+ px)
             !    Eso2 = xk*xe*so2g *patm
             ! equivalent new code
             !    heso2 = xk + xk*xe/hplus * xk*xe*x2/hplus**2
             !    so2g = xso2/(1 + px)
             !         = xso2/(1 + heso2*ra*tz*xl)
             !         = xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
             !    eso2 = so2g*xk*xe*patm
             !          = xk*xe*patm*xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
             !          = ( fact1_so2    )/(1 + fact2_so2 *(1 + (fact3_so2/hplus)*(1 + fact4_so2/hplus)
             !    [hso3-] + 2*[so3--] = (eso2/hplus)*(1 + 2*x2/hplus)
             xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
             xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
             x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )
             fact1_so2 = xk*xe*patm*xso2(i,k)
             fact2_so2 = xk*GAS_CONSTANT_L_ATM_MOL_K*temperature(i,k)*xl
             fact3_so2 = xe
             fact4_so2 = x2

             !-----------------------------------------------------------------
             !          ... nh3
             !-----------------------------------------------------------------
             ! previous code
             !    henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/water_dissociation_constant)
             !    px = henh3(i,k) * GAS_CONSTANT_L_ATM_MOL_K * tz * xl
             !    nh3g = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
             !    Enh3 = xk*xe*nh3g/water_dissociation_constant *patm
             ! equivalent new code
             !    henh3 = xk + xk*xe*hplus/water_dissociation_constant
             !    nh3g = xnh34/(1 + px)
             !         = xnh34/(1 + henh3*ra*tz*xl)
             !         = xnh34/(1 + xk*ra*tz*xl*(1 + xe*hplus/water_dissociation_constant)
             !    enh3 = nh3g*xk*xe*patm/water_dissociation_constant
             !          = ((xk*xe*patm/water_dissociation_constant)*xnh34)/
             !              (1 + xk*ra*tz*xl*(1 + xe*hplus/water_dissociation_constant)
             !          = ( fact1_nh3            )/(1 + fact2_nh3  *(1 + fact3_nh3*hplus)
             !    [nh4+] = enh3*hplus
             xk = 58._r8   *EXP( 4085._r8*work1(i) )
             xe = 1.7e-5_r8*EXP( -4325._r8*work1(i) )

             fact1_nh3 = (xk*xe*patm/water_dissociation_constant)*(xnh3(i,k)+xnh4(i,k))
             fact2_nh3 = xk*GAS_CONSTANT_L_ATM_MOL_K*temperature(i,k)*xl
             fact3_nh3 = xe/water_dissociation_constant

             !-----------------------------------------------------------------
             !        ... h2o effects
             !-----------------------------------------------------------------
             Eh2o = water_dissociation_constant

             !-----------------------------------------------------------------
             !        ... co2 effects
             !-----------------------------------------------------------------
             co2g = 330.e-6_r8                            !330 ppm = 330.e-6 atm
             xk = 3.1e-2_r8*EXP( 2423._r8*work1(i) )
             xe = 4.3e-7_r8*EXP(-913._r8 *work1(i) )
             Eco2 = xk*xe*co2g  *patm

             !-----------------------------------------------------------------
             !         ... so4 effect
             !-----------------------------------------------------------------
             Eso4 = xso4(i,k)*air_number_density(i,k)   &         ! /cm3(a)
                  * (1.e3_r8/AVOGADRO) / xl


             !-----------------------------------------------------------------
             ! now use bisection method to solve electro-neutrality equation
             !
             ! during the iteration loop,
             !    yph_lo = lower ph value that brackets the root (i.e., correct ph)
             !    yph_hi = upper ph value that brackets the root (i.e., correct ph)
             !    yph    = current ph value
             !    yposnet_lo and yposnet_hi = net positive ions for
             !       yph_lo and yph_hi
             !-----------------------------------------------------------------
             do iter = 1, MAX_ITERATIONS

                if (.not. present(specified_ph)) then
                   if (iter == 1) then
                      ! 1st iteration ph = lower bound value
                      yph_lo = 2.0_r8
                      yph_hi = yph_lo
                      yph = yph_lo
                   else if (iter == 2) then
                      ! 2nd iteration ph = upper bound value
                      yph_hi = 7.0_r8
                      yph = yph_hi
                   else
                      ! later iteration ph = mean of the two bracketing values
                      yph = 0.5_r8*(yph_lo + yph_hi)
                   end if
                else
                   yph = specified_ph
                end if

                ! calc current [H+] from ph
                xph(i,k) = 10.0_r8**(-yph)


                !-----------------------------------------------------------------
                !        ... hno3
                !-----------------------------------------------------------------
                Ehno3 = fact1_hno3/(1.0_r8 + fact2_hno3*(1.0_r8 + fact3_hno3/xph(i,k)))

                !-----------------------------------------------------------------
                !          ... so2
                !-----------------------------------------------------------------
                Eso2 = fact1_so2/(1.0_r8 + fact2_so2*(1.0_r8 + (fact3_so2/xph(i,k)) &
                     *(1.0_r8 +  fact4_so2/xph(i,k))))

                !-----------------------------------------------------------------
                !          ... nh3
                !-----------------------------------------------------------------
                Enh3 = fact1_nh3/(1.0_r8 + fact2_nh3*(1.0_r8 + fact3_nh3*xph(i,k)))

                tmp_nh4  = Enh3 * xph(i,k)
                tmp_hso3 = Eso2 / xph(i,k)
                tmp_so3  = tmp_hso3 * 2.0_r8*fact4_so2/xph(i,k)
                tmp_hco3 = Eco2 / xph(i,k)
                tmp_oh   = Eh2o / xph(i,k)
                tmp_no3  = Ehno3 / xph(i,k)
                tmp_so4 = cldconc%so4_fact*Eso4
                tmp_pos = xph(i,k) + tmp_nh4
                tmp_neg = tmp_oh + tmp_hco3 + tmp_no3 + tmp_hso3 + tmp_so3 + tmp_so4

                ynetpos = tmp_pos - tmp_neg


                ! yposnet = net positive ions/charge
                ! if the correct ph is bracketed by yph_lo and yph_hi (with yph_lo < yph_hi),
                !    then you will have yposnet_lo > 0 and yposnet_hi < 0
                converged = .false.
                if (iter > 2) then
                   if (ynetpos == 0.0_r8) then
                      ! the exact solution was found (very unlikely)
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   else if (ynetpos >= 0.0_r8) then
                      ! net positive ions are >= 0 for both yph and yph_lo
                      !    so replace yph_lo with yph
                      yph_lo = yph
                      ynetpos_lo = ynetpos
                   else
                      ! net positive ions are <= 0 for both yph and yph_hi
                      !    so replace yph_hi with yph
                      yph_hi = yph
                      ynetpos_hi = ynetpos
                   end if

                   if (abs(yph_hi - yph_lo) .le. 0.005_r8) then
                      ! |yph_hi - yph_lo| <= convergence criterion, so set
                      !    final ph to their midpoint and exit
                      ! (.005 absolute error in pH gives .01 relative error in H+)
                      tmp_hp = xph(i,k)
                      yph = 0.5_r8*(yph_hi + yph_lo)
                      xph(i,k) = 10.0_r8**(-yph)
                      converged = .true.
                      exit
                   else
                      ! do another iteration
                      converged = .false.
                   end if

                else if (iter == 1) then
                   if (ynetpos <= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use the lower bound
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_lo = ynetpos

                else ! (iter == 2)
                   if (ynetpos >= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use they upper bound
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_hi = ynetpos
                end if

             end do ! iter

             if( .not. converged ) then
                write(iulog,*) 'setsox: pH failed to converge @ (',i,',',k,')'
             end if
          else
             xph(i,k) =  1.e-7_r8
          end if
       end do col_loop0
    end do ver_loop0 ! end pver loop for STEP 0

    !==============================================================
    !          ... Now use the actual PH
    !==============================================================
    ver_loop1: do k = 1,pver
       col_loop1: do i = 1,ncol
          work1(i) = 1._r8 / temperature(i,k) - 1._r8 / 298._r8
          xl = cldconc%xlwc(i,k)

          ! This should be dividing by 101325, not 101300, but changing it breaks the tests
          patm = midpoint_pressure(i,k) / 101300._r8        ! press is in pascal

          !-----------------------------------------------------------------------
          !        ... hno3
          !-----------------------------------------------------------------------
          xk = 2.1e5_r8 *EXP( 8700._r8*work1(i) )
          xe = 15.4_r8
          hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !        ... h2o2
          !-----------------------------------------------------------------
          xk = 7.4e4_r8   *EXP( 6621._r8*work1(i) )
          xe = 2.2e-12_r8 *EXP(-3730._r8*work1(i) )
          heh2o2(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !         ... so2
          !-----------------------------------------------------------------
          xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
          xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
          x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )

          heso2(i,k)  = xk*(1._r8 + xe/xph(i,k)*(1._r8 + x2/xph(i,k)))

          !-----------------------------------------------------------------
          !          ... nh3
          !-----------------------------------------------------------------
          xk = 58._r8   *EXP( 4085._r8*work1(i) )
          xe = 1.7e-5_r8*EXP(-4325._r8*work1(i) )
          henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/water_dissociation_constant)

          !-----------------------------------------------------------------
          !        ... o3
          !-----------------------------------------------------------------
          xk = 1.15e-2_r8 *EXP( 2560._r8*work1(i) )
          heo3(i,k) = xk

          !------------------------------------------------------------------------
          !       ... for Ho2(g) -> H2o2(a) formation
          !           schwartz JGR, 1984, 11589
          !------------------------------------------------------------------------
          ho2s = kh0*xho2(i,k)*patm*(1._r8 + kh1/xph(i,k))  ! ho2s = ho2(a)+o2-
          r1h2o2 = (kh2 + kh3*kh1/xph(i,k)) / ((1._r8 + kh1/xph(i,k))**2)*ho2s*ho2s ! prod(h2o2) in mole/L(w)/s

          if ( cloud_borne ) then
             r2h2o2 = r1h2o2*xl                  & ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  / (1.e3_r8/AVOGADRO)*1.e+6_r8  & ! correct a bug here ????
                  / (midpoint_pressure(i,k)/(BOLTZMANN*temperature(i,k)))
          else
             r2h2o2 = r1h2o2*xl         & ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  * (1.e3_r8/AVOGADRO)  & ! mole/fm3(a)/s * 1.e-3       = mole/cm3(a)/s
                  / (midpoint_pressure(i,k)/(BOLTZMANN*temperature(i,k))) ! /cm3(a)/s    / air-den     = mix-ratio/s
          endif

          if ( .not. cloud_borne) then    ! this seems to be specific to aerosols that are not cloud borne
             xh2o2(i,k) = xh2o2(i,k) + r2h2o2*time_step         ! updated h2o2 by het production
          endif

          !-----------------------------------------------
          !       ... Partioning
          !-----------------------------------------------

          !-----------------------------------------------------------------
          !        ... hno3
          !-----------------------------------------------------------------
          px = hehno3(i,k) * GAS_CONSTANT_L_ATM_MOL_K * temperature(i,k) * xl
          hno3g(i,k) = (xhno3(i,k)+xno3(i,k))/(1._r8 + px)

          !------------------------------------------------------------------------
          !        ... h2o2
          !------------------------------------------------------------------------
          px = heh2o2(i,k) * GAS_CONSTANT_L_ATM_MOL_K * temperature(i,k) * xl
          h2o2g =  xh2o2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... so2
          !------------------------------------------------------------------------
          px = heso2(i,k) * GAS_CONSTANT_L_ATM_MOL_K * temperature(i,k) * xl
          so2g =  xso2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... o3
          !------------------------------------------------------------------------
          px = heo3(i,k) * GAS_CONSTANT_L_ATM_MOL_K * temperature(i,k) * xl
          o3g =  xo3(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... nh3
          !------------------------------------------------------------------------
          px = henh3(i,k) * GAS_CONSTANT_L_ATM_MOL_K * temperature(i,k) * xl
          if (nh3%exists()) then
             nh3g(i,k) = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
          else
             nh3g(i,k) = 0._r8
          endif

          !-----------------------------------------------
          !       ... Aqueous phase reaction rates
          !           SO2 + H2O2 -> SO4
          !           SO2 + O3   -> SO4
          !-----------------------------------------------

          !------------------------------------------------------------------------
          !       ... S(IV) (HSO3) + H2O2
          !------------------------------------------------------------------------
          k_siv_h2o2 = 8.e4_r8 * EXP( -3650._r8*work1(i) )  &
               / (.1_r8 + xph(i,k))

          !------------------------------------------------------------------------
          !        ... S(IV)+ O3
          !------------------------------------------------------------------------
          k_siv_o3   = 4.39e11_r8 * EXP(-4131._r8/temperature(i,k))  &
               + 2.56e3_r8  * EXP(-996._r8 /temperature(i,k)) /xph(i,k)

          !-----------------------------------------------------------------
          !       ... Prediction after aqueous phase
          !       so4
          !       When Cloud is present
          !
          !       S(IV) + H2O2 = S(VI)
          !       S(IV) + O3   = S(VI)
          !
          !       reference:
          !           (1) Seinfeld
          !           (2) Benkovitz
          !-----------------------------------------------------------------

          !............................
          !       S(IV) + H2O2 = S(VI)
          !............................

          IF (XL .ge. MINIMUM_CLOUD_LIQUID_WATER) THEN    !! WHEN CLOUD IS PRESENTED

             if (cloud_borne) then
                patm_x = patm
             else
                patm_x = 1._r8
             endif

             if (cloud_borne) then
                dso4_dt = k_siv_h2o2 * 7.4e4_r8*EXP(6621._r8*work1(i)) * h2o2g * patm_x &
                     * 1.23_r8 *EXP(3120._r8*work1(i)) * so2g * patm_x
             else
                dso4_dt = k_siv_h2o2 * heh2o2(i,k) * h2o2g * patm_x  &
                     * heso2(i,k)  * so2g  * patm_x    ! [M/s]
             endif

             dso4_dt = dso4_dt & ! [M/s] = [mole/L(w)/s]
                  * xl & ! [mole/L(a)/s]
                  / (1.e3_r8/AVOGADRO) & ! [/L(a)/s]
                  / air_number_density(i,k)


             delta_concentration = dso4_dt*time_step
             delta_concentration = max(delta_concentration, 1.e-30_r8)

             xso4_init(i,k)=xso4(i,k)

             IF (xh2o2(i,k) .gt. xso2(i,k)) THEN
                if (delta_concentration .gt. xso2(i,k)) then
                   xso4(i,k)=xso4(i,k)+xso2(i,k)
                   if (cloud_borne) then
                      xh2o2(i,k)=xh2o2(i,k)-xso2(i,k)
                      xso2(i,k)=1.e-20_r8
                   else       ! ???? bug ????
                      xso2(i,k)=1.e-20_r8
                      xh2o2(i,k)=xh2o2(i,k)-xso2(i,k)
                   endif
                else
                   xso4(i,k)  = xso4(i,k)  + delta_concentration
                   xh2o2(i,k) = xh2o2(i,k) - delta_concentration
                   xso2(i,k)  = xso2(i,k)  - delta_concentration
                end if

             ELSE
                if (delta_concentration  .gt. xh2o2(i,k)) then
                   xso4(i,k)=xso4(i,k)+xh2o2(i,k)
                   xso2(i,k)=xso2(i,k)-xh2o2(i,k)
                   xh2o2(i,k)=1.e-20_r8
                else
                   xso4(i,k)  = xso4(i,k)  + delta_concentration
                   xh2o2(i,k) = xh2o2(i,k) - delta_concentration
                   xso2(i,k)  = xso2(i,k)  - delta_concentration
                end if
             END IF

             if (cloud_borne) then
                change_in_aq_so4_mixing_ratio(i,k)  =  xso4(i,k) - xso4_init(i,k)
             endif
             !...........................
             !       S(IV) + O3 = S(VI)
             !...........................

             dso4_dt = k_siv_o3 * heo3(i,k)*o3g*patm_x * heso2(i,k)*so2g*patm_x  ! [M/s]

             dso4_dt = dso4_dt               &    ! [M/s] =  [mole/L(w)/s]
                  * xl                 &    ! [mole/L(a)/s]
                  / (1.e3_r8/AVOGADRO) &    ! [/L(a)/s]
                  / air_number_density(i,k) ! [mixing ratio/s]

             delta_concentration = dso4_dt*time_step
             delta_concentration = max(delta_concentration, 1.e-30_r8)

             xso4_init(i,k)=xso4(i,k)

             if (delta_concentration .gt. xso2(i,k)) then
                xso4(i,k) = xso4(i,k) + xso2(i,k)
                xso2(i,k) = 1.e-20_r8
             else
                xso4(i,k) = xso4(i,k) + delta_concentration
                xso2(i,k) = xso2(i,k) - delta_concentration
             end if

          END IF !! WHEN CLOUD IS PRESENTED

       end do col_loop1
    end do ver_loop1

    call sox_cldaero_update( state, &
          pbuf, ncol, lchnk, loffset, time_step, mean_mass, pressure_thickness, midpoint_pressure, temperature, cloud_droplet_number, cloud_fraction, air_mass_density_kg_l, cldconc%xlwc, &
          change_in_aq_so4_mixing_ratio, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, cloud_borne_aerosol_vmr, species_vmr, &
          aq_so4_production, aq_h2so4_production, aq_so4_production_from_h2o2, aq_so4_production_from_o3, aqso4_h2o2_3d=aq_so4_production_from_h2o2_3d, aqso4_o3_3d=aq_so4_production_from_o3_3d )

    ph_times_cloud_water(:,:) = 0._r8
    do k = 1, pver
       do i = 1, ncol
          if (cloud_fraction(i,k)>=1.e-5_r8 .and. cloud_water(i,k)>=1.e-8_r8) then
             ph_times_cloud_water(i,k) = -1._r8*log10(xph(i,k)) * cloud_water(i,k)
          endif
       end do
    end do

    call sox_cldaero_destroy_obj(cldconc)

  end subroutine calculate

!-------------------------------------------------------------------------------
!
! Support routines to be moved to separate modules
!
!-------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------
   ! Creates a cloud species object with the given name
   !
   ! The name is used to determine the species index in the state arrays.
   ! If the species is not found, the index is set to CLOUD_INDEX_UNDEFINED.
   function cloud_species_constructor( species_name ) result( this )
   
      use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
   
      type(cloud_species_t) :: this
      character(len=*), intent(in) :: species_name
   
      this%name_ = species_name
      this%state_index_ = get_inv_ndx( species_name )
      this%is_constant_ = this%state_index_ > 0
      if ( .not. this%is_constant_ ) &
          this%state_index_ = get_spc_ndx( species_name )
      if ( this%state_index_ <= 0 ) this%state_index_ = CLOUD_INDEX_UNDEFINED
   
   end function cloud_species_constructor

   !-------------------------------------------------------------------------------
   ! Returns whether a cloud species is defined
   logical function cloud_species_exists( this )
   
      class(cloud_species_t), intent(in) :: this
   
      cloud_species_exists = this%state_index_ .ne. CLOUD_INDEX_UNDEFINED
   
   end function cloud_species_exists

   !-------------------------------------------------------------------------------
   ! Returns the mixing ratio for a cloud species
   !
   ! Constant species use the fixed number concentration and the air density
   ! to calculate the mixing ratio. Time-varying species simply return the
   ! mixing ratio from the state array.
   !
   ! For species that are not defined, the mixing ratio is set to zero.
   subroutine cloud_species_get_mixing_ratio( this, mixing_ratios, &
       fixed_concentrations, air_number_density, mixing_ratio )

      class(cloud_species_t), intent(in) :: this
      real(r8), intent(in)  :: mixing_ratios(:,:,:)        ! all mixing ratios (mol mol-1) [column, layer, species]
      real(r8), intent(in)  :: fixed_concentrations(:,:,:) ! all fixed concentrations (# cm-3) [column, layer, species]
      real(r8), intent(in)  :: air_number_density(:,:)     ! air density (# cm-3) [column, layer]
      real(r8), intent(out) :: mixing_ratio(:,:)           ! species mixing ratio (mol mol-1) [column, layer]
      
      if ( this%state_index_ == CLOUD_INDEX_UNDEFINED ) then
         mixing_ratio(:,:) = 0._r8
         return
      end if
      if ( this%is_constant_ ) then
         mixing_ratio(:,:) = fixed_concentrations(:,:,this%state_index_) &
                             / air_number_density(:,:)
      else
         mixing_ratio(:,:) = mixing_ratios(:,:,this%state_index_)
      end if

   end subroutine cloud_species_get_mixing_ratio

end module cloud_aqueous_chemistry

