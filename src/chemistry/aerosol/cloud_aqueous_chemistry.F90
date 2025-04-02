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
  use chemistry_test_data

  implicit none

  private
  public :: initialize, calculate
  public :: do_cloud_aqueous_chemistry

  logical :: do_cloud_aqueous_chemistry = .false.

  integer, parameter :: CLOUD_INDEX_UNDEFINED = -1

  !> @brief Cloud chemistry species information
  type :: cloud_species_t
    character(len=:), allocatable :: name_
    integer :: state_index_ = CLOUD_INDEX_UNDEFINED ! index in the state vector or the fixed concentrations array
    logical :: is_constant_ = .false.
    real(r8) :: default_mixing_ratio_ = 0.0_r8 ! mol mol-1
  contains
    procedure :: exists => cloud_species_exists
    procedure :: mixing_ratio => cloud_species_get_mixing_ratio
  end type cloud_species_t

  interface cloud_species_t
    module procedure :: cloud_species_constructor
  end interface cloud_species_t

  !> @brief Van 't Hoff equation parameters
  !! Used to calculate the temperature dependence of Henry's Law constants
  !! The equation is given by:
  !! A * exp(-B * (1/T - 1/T_ref))
  !! where A and B are the parameters, T_ref is the reference temperature,
  !! and T is the current temperature. 
  type :: van_t_hoff_t
    real(r8) :: A_ = 0.0_r8
    real(r8) :: B_ = 0.0_r8
  end type van_t_hoff_t

  !> @brief Henry's Law types
  integer, parameter :: HENRYS_LAW_UNDEFINED = 0
  integer, parameter :: HENRYS_LAW_MONOPROTIC_ACID = 1
  integer, parameter :: HENRYS_LAW_DIPROTIC_ACID = 2
  integer, parameter :: HENRYS_LAW_BASE = 3
  integer, parameter :: HENRYS_LAW_NEUTRAL = 4
  integer, parameter :: HENRYS_LAW_HO2 = 5
  integer, parameter :: HENRYS_LAW_CO2 = 6

  !> @brief Henry's Law parameters for acids and bases
  !! Calculates an equilibrium concentration of the dissociated acid based on the
  !! pH of the solution. The equilibrium constant is given by:
  !! K_eq = [H+][A-] (mono-protic acid)
  !! K_eq = [H+]([HA-] + 2[A--]) (di-protic acid)
  !! K_eq = [BH+]/[H+] (base)
  !! and has units of mol^2 L^-2.
  !! Also, calculates the gas-phase mixing ratio (mol mol-1) of the acid based
  !! on the Henry's Law constant and the total mixing ratio.
  type :: henrys_law_t
    integer :: type_ = HENRYS_LAW_UNDEFINED
    type(van_t_hoff_t) :: partitioning_factor_       ! A = mol L-1 atm-1, B = K
    type(van_t_hoff_t) :: first_dissociation_factor_ ! A = mol L-1, B = K
    type(van_t_hoff_t) :: second_dissociation_factor_! A = mol L-1, B = K
    type(van_t_hoff_t) :: protonation_factor_        ! A = mol L-1, B = K
    type(van_t_hoff_t) :: ho2__hplus_o2minus_        ! A = mol2 L-2, B = K
    type(van_t_hoff_t) :: ho2_ho2__h2o2_o2_          ! A = mol2 L-2, B = K
    type(van_t_hoff_t) :: ho2_o2minus__h2o2_o2_      ! A = mol2 L-2, B = K
    real(r8) :: reference_temperature_ = 298.0_r8    ! K
    real(r8), allocatable :: terms_(:,:,:) ! 1: mol2 L-2, 2: unitless, 3: mol L-1, 4: mol L-1 (term, column, level)
  contains
    procedure :: set_conditions => henrys_law_set_conditions
    procedure :: equilibrium_constant => henrys_law_equilibrium_constant
    procedure :: gas_phase_mixing_ratio => henrys_law_gas_phase_mixing_ratio
    procedure :: partitioned_species_concentration => henrys_law_partitioned_species_concentration
    procedure :: effective_henrys_law_constant => henrys_law_effective_henrys_law_constant
    procedure :: reaction_rate => henrys_law_reaction_rate
  end type henrys_law_t

  interface henrys_law_t
    module procedure :: henrys_law_constructor
  end interface henrys_law_t
    
  type(cloud_species_t) :: so2, nh3, hno3, h2o2, o3, ho2, msa, so4, h2so4, co2

  ! TODO: Figure out what this flag is for
  logical :: cloud_borne = .false.

  ! Constants that should be moved to a common module
  real(r8), parameter :: AVOGADRO = 6.02214076e23_r8 ! mol-1
  real(r8), parameter :: BOLTZMANN = 1.380649e-23_r8 ! J K-1
  real(r8), parameter :: PASCAL_TO_ATM = 1.0_r8/101325.0_r8 ! atm Pa-1
  real(r8), parameter :: M3_TO_L = 1.0e3_r8 ! L m-3
  real(r8), parameter :: GAS_CONSTANT = BOLTZMANN*AVOGADRO ! J K-1 mol-1
  real(r8), parameter :: GAS_CONSTANT_L_ATM_MOL_K = BOLTZMANN*AVOGADRO*M3_TO_L*PASCAL_TO_ATM ! L atm K-1 mol-1
  real(r8), parameter :: GAS_CONSTANT_DRY_AIR_J_KG_K = 287.0_r8 ! J kg-1 K-1
  real(r8), parameter :: SMALL_NUMBER = 1.0e-30_r8 ! unitless
  real(r8), parameter :: WATER_DISSOCIATION_CONSTANT = 1.0e-14_r8 ! mol2/L2  [H+][OH-]

contains

  !-----------------------------------------------------------------------
  !	Prepares for cloud aqueous chemistry
  ! - Looks up cloud chemistry species
  ! - Determines if enough species are present to perform cloud chemistry
  !-----------------------------------------------------------------------
  subroutine initialize()

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

    so2   = cloud_species_t( 'SO2'   )
    nh3   = cloud_species_t( 'NH3'   )
    hno3  = cloud_species_t( 'HNO3'  )
    h2o2  = cloud_species_t( 'H2O2'  )
    o3    = cloud_species_t( 'O3'    )
    ho2   = cloud_species_t( 'HO2'   )
    msa   = cloud_species_t( 'MSA'   )
    so4   = cloud_species_t( 'SO4'   )
    h2so4 = cloud_species_t( 'H2SO4' )
    co2   = cloud_species_t( 'CO2'   )

    do_cloud_aqueous_chemistry = so2%exists() .and. h2o2%exists() .and. &
                                 o3%exists() .and. ho2%exists() .and. co2%exists()
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
    real(r8), target, intent(in)    :: cloud_water(:,:)                ! cloud liquid water content (kg_water/kg_air)
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

    ! Change in aqueous sulfate volume mixing ratio over current time step (mol mol-1)
    real(r8) :: change_in_aq_so4_mixing_ratio(ncol,pver)

    ! Henry's Law Partitioning calculators
    type(henrys_law_t) :: hl_hno3, hl_so2, hl_nh3, hl_co2, hl_h2o2, hl_ho2, hl_o3

    ! Equilibrium constants used to determine condensed phase ion concentrations [various units]
    real(r8) :: Eso2, Ehno3, Eco2, Enh3

    ! Calculated gas-phase mixing ratios (mol mol-1)
    real(r8) :: hno3g(ncol,pver), nh3g(ncol,pver), so2g, h2o2g, o3g

    ! Calculated conversion from mol(X) L(H2O)-1 to mol(X) mol(air)-1 based on
    ! cloud liquid water content
    real(r8), dimension(ncol,pver) :: molar_to_mixing_ratio ! L mol-1

    ! Mid-point pressure in atm
    real(r8) :: patm
    
    ! Oxidation reaction parameters
    real(r8) :: dh2o2_dt    ! reaction rate for formation of H2O2(aq) (mol mol-1 s-1)
    real(r8) :: k_siv_h2o2  ! rate constant for reaction of S(IV) with H2O2 (L mol-1 s-1)
    real(r8) :: k_siv_o3    ! rate constant for reaction of S(IV) with O3 (L mol-1 s-1)
    real(r8) :: dso4_dt     ! rate of change of SO4 (mol(SIV) mol(air)-1 s-1)
    real(r8) :: delta_concentration ! change in concentration of S(IV) over the timestep (mol(SIV) mol(air)-1)

    ! index variables and couters
    integer  :: k, i, iter

    ! volume mixing ratios for cloud chemistry species (mol mol-1)
    real(r8), dimension(ncol,pver) ::  xhno3, xh2o2, xso2, xso4, xno3, &
         xnh3, xnh4, xo3, xho2, xh2so4, xmsa, xco2, xso4_init

    ! [H+] concentration in the aqueous phase (mol L-1)
    real(r8), dimension(ncol,pver) :: h_plus_conc

    ! air mass density (kg L-1)
    real(r8), dimension(ncol,pver) :: air_mass_density_kg_l
    
    ! Effective Henry's Law constants (mol L-1 atm-1)
    real(r8), dimension(ncol,pver) :: Heff_h2o2, Heff_so2, Heff_o3

    ! Aqueous-phase concentrations [mol(X) L(air)-1]
    real(r8) :: h2o2_aq, so2_aq, o3_aq, so4_aq

    logical :: converged

    ! Information about cloud composition
    type(cldaero_conc_t), pointer :: cloud_composition

    ! TODO: Skipping renaming of these in anticipation of partitioning/pH estimation structs
    real(r8) :: tmp_hp, tmp_hso3, tmp_hco3, tmp_nh4, tmp_no3
    real(r8) :: tmp_oh, tmp_so3, tmp_so4
    real(r8) :: tmp_neg, tmp_pos
    real(r8) :: yph, yph_lo, yph_hi
    real(r8) :: ynetpos, ynetpos_lo, ynetpos_hi

    !==================================================================
    ! Set partitioning parameters
    !==================================================================

    ! HNO3 partitioning parameters
    ! NOTE: The partitioning factor's A value is in mol m-3 Pa-1 in the reference below.
    !       This appears to have been roughly converted to mol L-1 atm-1 by multiplying by 1.e2,
    !       but this really should have been multiplied by 101325 / 1000 (Pa m3 atm-1 L-1).
    !       FUTURE_ANSWER_CHANGING_MODIFICATION
    !
    ! Sander, R., 2015. Compilation of Henry's law constants (version 4) for water as solvent.
    ! Atmos. Chem. Phys., 15, 4399–4981, 2015. DOI: 10.5194/acp-15-4399-2015
    !
    ! However, this reference cites Schwartz and White (1981) for these values, but gives their
    ! units in M atm-1 (Schwartz and White, 1981 doesn't appear to be available online):
    ! Xiaohong Yao, Tsz Yan Ling, Ming Fang, Chak K. Chan,
    ! Comparison of thermodynamic predictions for in situ pH in PM2.5,
    ! Atmospheric Environment, Volume 40, Issue 16, 2006, Pages 2835-2844, ISSN 1352-2310,
    ! DOI: https://doi.org/10.1016/j.atmosenv.2006.01.006.
    !
    ! JPL 19-5 () Lists several models for determining H_eff for HNO3 in acidic solutions,
    ! but doesn't seem to include equations in the text:
    !   "These models and the Carslaw et al. review should be consulted for plots/predictions
    !    of H* for HNO3 and HCl in strong acid solutions over the atmospheric temperature range."
    hl_hno3 = henrys_law_t( ncol, pver )
    hl_hno3%type_ = HENRYS_LAW_MONOPROTIC_ACID
    hl_hno3%partitioning_factor_%A_ = 2.1e5_r8      ! mol m-3 Pa-1 * 100 OR mol L-1 atm-1 ?
    hl_hno3%partitioning_factor_%B_ = 8700.0_r8     ! K
    hl_hno3%first_dissociation_factor_%A_ = 15.4_r8 ! M NOTE: This is 15.0 in Yao et al. (2006)
    hl_hno3%first_dissociation_factor_%B_ = 0.0_r8  ! K

    ! SO2 partitioning parameters
    !
    ! The values used here as the same as those used in the following paper (probably from some other source):
    ! Lu et al. Assessmment of the impacts of cloud chemistry on surface SO2 and sulfate levels in typical regions of China
    ! Atmos. Chem. Phys., 23, 8021–8037, 2023
    ! https://doi.org/10.5194/acp-23-8021-2023
    ! 
    ! Chameides, W. L. (1984), The photochemistry of a remote marine stratiform cloud,
    ! J. Geophys. Res., 89(D3), 4739–4755, doi:10.1029/JD089iD03p04739.
    ! Original source: National Bureau of Standards [1965]
    !
    ! All the values for the Henry's Law constants in:
    ! Sander, R., 2015. Compilation of Henry's law constants (version 4) for water as solvent.
    ! Atmos. Chem. Phys., 15, 4399–4981, 2015. DOI: 10.5194/acp-15-4399-2015
    hl_so2 = henrys_law_t( ncol, pver )
    hl_so2%type_ = HENRYS_LAW_DIPROTIC_ACID
    hl_so2%partitioning_factor_%A_ = 1.23_r8           ! M atm-1
    hl_so2%partitioning_factor_%B_ = 3120.0_r8         ! K
    hl_so2%first_dissociation_factor_%A_ = 1.7e-2_r8   ! M
    hl_so2%first_dissociation_factor_%B_ = 2090.0_r8   ! K
    hl_so2%second_dissociation_factor_%A_ = 6.0e-8_r8  ! M
    hl_so2%second_dissociation_factor_%B_ = 1120.0_r8  ! K

    ! NH3 partitioning parameters
    !
    ! Chameides, W. L. (1984), The photochemistry of a remote marine stratiform cloud,
    ! J. Geophys. Res., 89(D3), 4739–4755, doi:10.1029/JD089iD03p04739.
    ! Original source: National Bureau of Standards [1965]
    hl_nh3 = henrys_law_t( ncol, pver )
    hl_nh3%type_ = HENRYS_LAW_BASE
    hl_nh3%partitioning_factor_%A_ = 58.0_r8   ! M atm-1
    hl_nh3%partitioning_factor_%B_ = 4085.0_r8 ! K
    hl_nh3%protonation_factor_%A_ = 1.7e-5_r8  ! M
    hl_nh3%protonation_factor_%B_ = -4325.0_r8 ! K

    ! CO2 partitioning parameters
    !
    ! Chameides, W. L. (1984), The photochemistry of a remote marine stratiform cloud,
    ! J. Geophys. Res., 89(D3), 4739–4755, doi:10.1029/JD089iD03p04739.
    ! Original source: National Bureau of Standards [1965]
    hl_co2 = henrys_law_t( ncol, pver )
    hl_co2%type_ = HENRYS_LAW_CO2
    hl_co2%partitioning_factor_%A_ = 3.1e-2_r8        ! M atm-1 NOTE: 3.11e-2 in Chameides (1984)
    hl_co2%partitioning_factor_%B_ = 2423.0_r8        ! K
    hl_co2%first_dissociation_factor_%A_ = 4.3e-7_r8  ! M
    hl_co2%first_dissociation_factor_%B_ = -913.0_r8  ! K

    ! H2O2 paritioning parameters
    !
    ! NOTE: I can't find a source for these values.
    hl_h2o2 = henrys_law_t( ncol, pver )
    hl_h2o2%type_ = HENRYS_LAW_MONOPROTIC_ACID
    hl_h2o2%partitioning_factor_%A_ = 7.4e4_r8
    hl_h2o2%partitioning_factor_%B_ = 6621.0_r8
    hl_h2o2%first_dissociation_factor_%A_ = 2.2e-12_r8
    hl_h2o2%first_dissociation_factor_%B_ = -3730.0_r8

    ! HO2 partitioning parameters
    !       ... for Ho2(g) -> H2o2(a) formation
    ! Chameides, W. L. (1984), The photochemistry of a remote marine stratiform cloud,
    ! J. Geophys. Res., 89(D3), 4739–4755, doi:10.1029/JD089iD03p04739.
    ! Original sources: National Bureau of Standards [1965], Bielski [1978],
    !                  and Chalmeides and Davis [1982]
    ! Schwartz, S. E. (1984), Gas- and aqueous-phase chemistry of HO2 in Liquid Water Clouds,
    ! J. Geophys. Res., 89(D7), 11589–11598
    ! Original source: JPL-82, Bieliski [1978]
    hl_ho2 = henrys_law_t( ncol, pver )
    hl_ho2%type_ = HENRYS_LAW_HO2
    hl_ho2%partitioning_factor_%A_ = 9.0e3_r8    ! M atm-1 NOTE: This is from Chameides (1984)
    hl_ho2%partitioning_factor_%B_ = 0.0_r8      ! K
    hl_ho2%ho2__hplus_o2minus_%A_ = 2.05e-5_r8   ! M NOTE: This is in both Chameides (1984) and Schwartz (1984)
    hl_ho2%ho2__hplus_o2minus_%B_ = 0.0_r8       ! K
    hl_ho2%ho2_ho2__h2o2_o2_%A_ = 8.6e5_r8       ! L2 mol-2 s-1 NOTE: This is from Schwartz (1984)
    hl_ho2%ho2_ho2__h2o2_o2_%B_ = 0.0_r8         ! K
    hl_ho2%ho2_o2minus__h2o2_o2_%A_ = 1.0e8_r8   ! L mol-1 s-1 NOTE: This is from Schwartz (1984)
    hl_ho2%ho2_o2minus__h2o2_o2_%B_ = 0.0_r8     ! K

    ! O3 paritioning parameters
    !
    ! The values used here as the same as those used in the following paper (probably from some other source):
    ! Lu et al. Assessmment of the impacts of cloud chemistry on surface SO2 and sulfate levels in typical regions of China
    ! Atmos. Chem. Phys., 23, 8021–8037, 2023
    ! https://doi.org/10.5194/acp-23-8021-2023
    !
    ! Chameides, W. L. (1984), The photochemistry of a remote marine stratiform cloud,
    ! J. Geophys. Res., 89(D3), 4739–4755, doi:10.1029/JD089iD03p04739.
    ! Original source: National Bureau of Standards [1965]
    hl_o3 = henrys_law_t( ncol, pver )
    hl_o3%type_ = HENRYS_LAW_NEUTRAL
    hl_o3%partitioning_factor_%A_ = 1.15e-2_r8 ! M atm-1
    hl_o3%partitioning_factor_%B_ = 2560.0_r8  ! K

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

    cloud_composition => sox_cldaero_create_obj( cloud_fraction, cloud_borne_aerosol_vmr, &
                                         cloud_water, air_mass_density_kg_l, ncol, loffset )
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
    call co2%mixing_ratio(   species_vmr, fixed_concentrations, air_number_density, xco2   )

    h_plus_conc(:,:) = 10._r8**(-INITIAL_PH) ! initial [H+] (mol L-1)

    !-----------------------------------------------------------------
    !       ... Temperature dependent Henry constants
    !-----------------------------------------------------------------
    ver_loop0: do k = 1,pver                               !! pver loop for STEP 0
       col_loop0: do i = 1,ncol

          if (cloud_borne .and. cloud_fraction(i,k)>0._r8) then
             xso4(i,k) = cloud_composition%so4c(i,k) / cloud_fraction(i,k)
             xnh4(i,k) = cloud_composition%nh4c(i,k) / cloud_fraction(i,k)
             xno3(i,k) = cloud_composition%no3c(i,k) / cloud_fraction(i,k)
          endif

          molar_to_mixing_ratio(i,k) = cloud_water(i,k) &
                                       * M3_TO_L * GAS_CONSTANT &
                                       * temperature(i,k) / midpoint_pressure(i,k)

          if( cloud_composition%xlwc(i,k) >= MINIMUM_CLOUD_LIQUID_WATER ) then

             patm = midpoint_pressure(i,k) * PASCAL_TO_ATM ! atm

             !-----------------------------------------------------------------
             ! Update Henry's Law calculators with current conditions
             !-----------------------------------------------------------------
             associate( lwc => cloud_composition%xlwc(i,k) )
               call hl_hno3%set_conditions( i, k, temperature(i,k), patm, lwc, xhno3(i,k)            )
               call hl_so2%set_conditions(  i, k, temperature(i,k), patm, lwc, xso2(i,k)             )
               call hl_nh3%set_conditions(  i, k, temperature(i,k), patm, lwc, xnh3(i,k) + xnh4(i,k) )
               call hl_co2%set_conditions(  i, k, temperature(i,k), patm, lwc, xco2(i,k)             )
               call hl_h2o2%set_conditions( i, k, temperature(i,k), patm, lwc, xh2o2(i,k)            )
               call hl_ho2%set_conditions(  i, k, temperature(i,k), patm, lwc, xho2(i,k)             )
               call hl_o3%set_conditions(   i, k, temperature(i,k), patm, lwc, xo3(i,k)              )
             end associate ! lwc
             !-----------------------------------------------------------------
             !         ... so4 effect
             !-----------------------------------------------------------------
             so4_aq = xso4(i,k) / molar_to_mixing_ratio(i,k) ! mol(so4)/L(w)

             !-----------------------------------------------------------------
             ! 21-mar-2011 chamges by rce
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
                h_plus_conc(i,k) = 10.0_r8**(-yph)

                !-----------------------------------------------------------------
                ! Get equilibrium constants for the current pH
                !-----------------------------------------------------------------
                call hl_hno3%equilibrium_constant( i, k, h_plus_conc(i,k), Ehno3 )
                call hl_so2%equilibrium_constant(  i, k, h_plus_conc(i,k), Eso2  )
                call hl_nh3%equilibrium_constant(  i, k, h_plus_conc(i,k), Enh3  )
                call hl_co2%equilibrium_constant(  i, k, h_plus_conc(i,k), Eco2  )

                tmp_nh4  = Enh3 * h_plus_conc(i,k)
                tmp_hso3 = Eso2 / h_plus_conc(i,k)
                tmp_so3  = tmp_hso3 * 2.0_r8*hl_so2%terms_(5,i,k) / h_plus_conc(i,k)
                tmp_hco3 = Eco2 / h_plus_conc(i,k)
                tmp_oh   = WATER_DISSOCIATION_CONSTANT / h_plus_conc(i,k)
                tmp_no3  = Ehno3 / h_plus_conc(i,k)
                tmp_so4 = cloud_composition%so4_fact*so4_aq
                tmp_pos = h_plus_conc(i,k) + tmp_nh4
                tmp_neg = tmp_oh + tmp_hco3 + tmp_no3 + tmp_hso3 + tmp_so3 + tmp_so4

                ynetpos = tmp_pos - tmp_neg

                ! yposnet = net positive ions/charge
                ! if the correct ph is bracketed by yph_lo and yph_hi (with yph_lo < yph_hi),
                !    then you will have yposnet_lo > 0 and yposnet_hi < 0
                converged = .false.
                if (iter > 2) then
                   if (ynetpos == 0.0_r8) then
                      ! the exact solution was found (very unlikely)
                      tmp_hp = h_plus_conc(i,k)
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
                      tmp_hp = h_plus_conc(i,k)
                      yph = 0.5_r8*(yph_hi + yph_lo)
                      h_plus_conc(i,k) = 10.0_r8**(-yph)
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
                      tmp_hp = h_plus_conc(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_lo = ynetpos

                else ! (iter == 2)
                   if (ynetpos >= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use they upper bound
                      tmp_hp = h_plus_conc(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_hi = ynetpos
                end if

             end do ! iter

             if( .not. converged ) then
                write(iulog,*) 'Cloud aqueous chemistry: pH failed to converge @ (',i,',',k,')'
             end if
          else
             h_plus_conc(i,k) =  1.e-7_r8 ! neutral pH
          end if
       end do col_loop0
    end do ver_loop0 ! end pver loop for STEP 0

    !==============================================================
    !          ... Now use the actual PH
    !==============================================================
    ver_loop1: do k = 1,pver
       col_loop1: do i = 1,ncol

          patm = midpoint_pressure(i,k) * PASCAL_TO_ATM ! atm

          ! Remove NH3/NH4 associated with NH4HSO4 in clouds
          if (cloud_borne .and. cloud_fraction(i,k)>0._r8) then
             xnh4(i,k) = xnh4(i,k) - cloud_composition%nh4c(i,k) / cloud_fraction(i,k)
          endif
          ! Account for H2O2(aq) in cloud water from HO2 uptake
          ! TODO: Investigate whether this should be done for cloud_borne aerosols
          if ( .not. cloud_borne ) then
             call hl_ho2%reaction_rate( i, k, h_plus_conc(i,k), xho2(i,k), dh2o2_dt )
             xh2o2(i,k) = xh2o2(i,k) + dh2o2_dt * time_step
          endif

          !-----------------------------------------------
          !       ... Partioning
          !-----------------------------------------------
          call hl_nh3%gas_phase_mixing_ratio(  i, k, h_plus_conc(i,k), xnh3(i,k)+xnh4(i,k),  nh3g(i,k)  )
          call hl_hno3%gas_phase_mixing_ratio( i, k, h_plus_conc(i,k), xhno3(i,k)+xno3(i,k), hno3g(i,k) )
          call hl_so2%gas_phase_mixing_ratio(  i, k, h_plus_conc(i,k), xso2(i,k),            so2g       )
          call hl_h2o2%gas_phase_mixing_ratio( i, k, h_plus_conc(i,k), xh2o2(i,k),           h2o2g      )
          call hl_o3%gas_phase_mixing_ratio(   i, k, h_plus_conc(i,k), xo3(i,k),             o3g        )

          !-----------------------------------------------
          !       ... Aqueous phase reaction rates
          !           SO2 + H2O2 -> SO4
          !           SO2 + O3   -> SO4
          !       reference:
          !           (1) Seinfeld
          !           (2) Benkovitz
          !-----------------------------------------------
          IF (cloud_composition%xlwc(i,k) .ge. MINIMUM_CLOUD_LIQUID_WATER) THEN    !! WHEN CLOUD IS PRESENTED

             call hl_so2%effective_henrys_law_constant( i, k, h_plus_conc(i,k), Heff_so2(i,k) )
             call hl_o3%effective_henrys_law_constant(  i, k, h_plus_conc(i,k), Heff_o3(i,k)  )
             
             call hl_h2o2%partitioned_species_concentration( i, k, h_plus_conc(i,k), xh2o2(i,k), patm, h2o2_aq )
             call hl_so2%partitioned_species_concentration(  i, k, h_plus_conc(i,k), xso2(i,k),  patm, so2_aq  )
             call hl_o3%partitioned_species_concentration(   i, k, h_plus_conc(i,k), xo3(i,k),   patm, o3_aq   )

             !------------------------------------------------------------------------
             !       ... S(IV) (HSO3) + H2O2
             ! (https://acp.copernicus.org/articles/24/1467/2024/acp-24-1467-2024.pdf)
             !------------------------------------------------------------------------
             ! FUTURE_ANSWER_CHANGING_MODIFICATION
             ! It is unclear why different logic is used here for cloud_borne/not cloud_borne
             ! Specifically, there is no explanation for why pressure should be 1.0 atm in the
             ! non-cloud-borne case, but the actual pressure should be used in the cloud-borne case.
             ! Also, in the rate calculations, the cloud-borne logic seems more consistent
             ! with the description in the reference, where the [SO2](aq) and [H2O2](aq)
             ! are calculated using the actual Henry's Law constants (not effective HLC) and
             ! the gas-phase mixing ratios, because the gas-phase concentrations were already
             ! calculated using H_eff:
             ! xH2O2(gas) = xH2O2(total) / (1 + H_eff*MIXING_RATIO_TO_ATM)
             ! [H2O2](aq) = xH2O2(gas) * H
             !
             ! In the S(IV) + O3 reaction, it also appears to be inconsistently using the
             ! total of [SO2](aq), [HSO3-](aq), and [SO3--](aq) in the reaction rate calculation
             k_siv_h2o2 = 8.e4_r8 * EXP( -3650._r8*(1._r8 / temperature(i,k) - 1._r8 / 298._r8) )  &
                          / (.1_r8 + h_plus_conc(i,k))
             if (cloud_borne) then
                ! NOTE: uses,
                !    [H2O2](aq) = [H2O2](aq)
                !    [SO2](aq)  = [SO2](aq)
                dso4_dt = k_siv_h2o2 * h2o2_aq * so2_aq * molar_to_mixing_ratio(i,k) ! mol mol-1 s-1
             else
                ! NOTE: uses,
                !    [H2O2](aq) = [H2O2](aq) + [HO2-](aq)
                !    [SO2](aq)  = [SO2](aq) + [HSO3-](aq) + [SO3--](aq)
                ! For some reason, assumes a pressure of 1 atm
                call hl_h2o2%effective_henrys_law_constant( i, k, h_plus_conc(i,k), Heff_h2o2(i,k) )
                dso4_dt = k_siv_h2o2 * Heff_h2o2(i,k) * h2o2g * Heff_so2(i,k) * so2g &
                          * molar_to_mixing_ratio(i,k)    ! [mol mol-1 s-1]
             endif
             delta_concentration = max(dso4_dt*time_step, SMALL_NUMBER)

! FUTURE_ANSWER_CHANGING_MODIFICATION - simplify and use consistent logic
#if 0
             delta_concentration = min(delta_concentration, min(xh2o2(i,k), xso2(i,k)))
             xso4(i,k)  = max(xso4(i,k)  + delta_concentration, SMALL_NUMBER)
             xh2o2(i,k) = max(xh2o2(i,k) - delta_concentration, SMALL_NUMBER)
             xso2(i,k)  = max(xso2(i,k)  - delta_concentration, 1.0e-20_r8) ! FUTURE_ANSWER_CHANGING_MODIFICATION Use SMALL_NUMBER?
             change_in_aq_so4_mixing_ratio(i,k) = delta_concentration
#else
             xso4_init(i,k)=xso4(i,k)

             IF (xh2o2(i,k) .gt. xso2(i,k)) THEN
                if (delta_concentration .gt. xso2(i,k)) then
                   xso4(i,k)=xso4(i,k)+xso2(i,k)
                   if (cloud_borne) then
                      xh2o2(i,k)=xh2o2(i,k)-xso2(i,k)
                      xso2(i,k)=1.e-20_r8 ! TODO: See if SMALL_NUMBER is more appropriate
                   else       ! ???? bug ????
                      ! FUTURE_ANSWER_CHANGING_MODIFICATION - this appears to be incorrect
                      xso2(i,k)=1.e-20_r8 ! TODO: See if SMALL_NUMBER is more appropriate
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

             change_in_aq_so4_mixing_ratio(i,k)  =  xso4(i,k) - xso4_init(i,k)
#endif
             !------------------------------------------------------------------------
             !       S(IV) + O3 = S(VI)
             ! (https://acp.copernicus.org/articles/24/1467/2024/acp-24-1467-2024.pdf)
             !------------------------------------------------------------------------
             ! NOTE: uses,
             !   [SO2](aq) = [SO2](aq) + [HSO3-](aq) + [SO3--](aq)
             !   [O3](aq)  = [O3](aq)
             ! FUTURE_ANSWER_CHANGING_MODIFICATION
             !    see if we should use [SO2](aq) here instead
             !    also, unclear why a pressure of 1 atm is used for non-cloud-borne case
             k_siv_o3   = 4.39e11_r8 * EXP(-4131._r8/temperature(i,k))  &
                          + 2.56e3_r8  * EXP(-996._r8 /temperature(i,k)) / h_plus_conc(i,k)
             if (cloud_borne) then
               dso4_dt = k_siv_o3 * o3_aq * Heff_so2(i,k)*so2g*patm * molar_to_mixing_ratio(i,k) ! [mol mol-1 s-1]
             else
               dso4_dt = k_siv_o3 * o3_aq * Heff_so2(i,k)*so2g * molar_to_mixing_ratio(i,k)      ! [mol mol-1 s-1]
             endif
             delta_concentration = max(dso4_dt*time_step, SMALL_NUMBER)
             delta_concentration = min(delta_concentration, xso2(i,k))
             xso4_init(i,k) = xso4(i,k)
             xso4(i,k) = max(xso4(i,k) + delta_concentration, SMALL_NUMBER)
             xso2(i,k) = max(xso2(i,k) - delta_concentration, 1.0e-20_r8) ! FUTURE_ANSWER_CHANGING_MODIFICATION Use SMALL_NUMBER?

          END IF !! WHEN CLOUD IS PRESENTED

       end do col_loop1
    end do ver_loop1

    call sox_cldaero_update( state, &
          pbuf, ncol, lchnk, loffset, time_step, mean_mass, pressure_thickness, midpoint_pressure, temperature, cloud_droplet_number, cloud_fraction, air_mass_density_kg_l, cloud_composition%xlwc, &
          change_in_aq_so4_mixing_ratio, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, cloud_composition%nh4c, cloud_composition%no3c, xmsa, xso2, xh2o2, cloud_borne_aerosol_vmr, species_vmr, &
          aq_so4_production, aq_h2so4_production, aq_so4_production_from_h2o2, aq_so4_production_from_o3, aqso4_h2o2_3d=aq_so4_production_from_h2o2_3d, aqso4_o3_3d=aq_so4_production_from_o3_3d )

    ph_times_cloud_water(:,:) = 0._r8
    do k = 1, pver
       do i = 1, ncol
          if (cloud_fraction(i,k)>=1.e-5_r8 .and. cloud_water(i,k)>=1.e-8_r8) then
             ph_times_cloud_water(i,k) = -1._r8*log10(h_plus_conc(i,k)) * cloud_water(i,k)
          endif
       end do
    end do

    call sox_cldaero_destroy_obj(cloud_composition)

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
   ! If a default mixing ratio is provided, it is used when the species is not
   ! found in the state arrays. Otherwise, the mixing ratio is set to zero.
   function cloud_species_constructor( species_name, default_mixing_ratio ) &
         result( this )
   
      use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
   
      type(cloud_species_t)          :: this
      character(len=*),   intent(in) :: species_name
      real(r8), optional, intent(in) :: default_mixing_ratio ! mol mol-1
   
      this%name_ = species_name
      this%state_index_ = get_inv_ndx( species_name )
      this%is_constant_ = this%state_index_ > 0
      if ( .not. this%is_constant_ ) &
          this%state_index_ = get_spc_ndx( species_name )
      if ( this%state_index_ <= 0 ) this%state_index_ = CLOUD_INDEX_UNDEFINED
      if ( present(default_mixing_ratio) ) then
         this%default_mixing_ratio_ = default_mixing_ratio
      end if
   
   end function cloud_species_constructor

   !-------------------------------------------------------------------------------
   ! Returns whether a cloud species is defined
   logical function cloud_species_exists( this )
   
      class(cloud_species_t), intent(in) :: this
   
      cloud_species_exists = (this%state_index_ .ne. CLOUD_INDEX_UNDEFINED) &
                             .or. (this%default_mixing_ratio_ > 0.0_r8)
   
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
         mixing_ratio(:,:) = this%default_mixing_ratio_
         return
      end if
      if ( this%is_constant_ ) then
         mixing_ratio(:,:) = fixed_concentrations(:,:,this%state_index_) &
                             / air_number_density(:,:)
      else
         mixing_ratio(:,:) = mixing_ratios(:,:,this%state_index_)
      end if

   end subroutine cloud_species_get_mixing_ratio

   !-------------------------------------------------------------------------------
   ! Constructor for the Henry's Law acid object
   !-------------------------------------------------------------------------------
   function henrys_law_constructor( number_of_columns, number_of_layers) &
         result( this )

      type(henrys_law_t) :: this
      integer, intent(in) :: number_of_columns
      integer, intent(in) :: number_of_layers

      ! TODO: Consider allocating the specific number of terms based on the
      !       type of partitioning (e.g., monoprotic, diprotic, etc.)
      allocate( this%terms_(6,number_of_columns,number_of_layers) )

   end function henrys_law_constructor

   !-------------------------------------------------------------------------------
   ! Updates the partitioning terms for the current conditions
   !
   ! Updates the partitioning terms for the current conditions. The partitioning
   ! terms are calculated based on the current temperature, cloud liquid water
   ! content, and the total mixing ratio of the species (gas and aqueous).
   !
   ! The partitioning terms are stored in the object for later use for
   ! calculating the gas-phase mixing ratio and the dissociated acid concentration.
   !
   ! For a monoprotic acid, the partitioning and first dissociation factors should
   ! be non-zero, and the second dissociation factor should be zero. For a diprotic
   ! acid, the partitioning, first, and second dissociation factors should be non-zero.
   ! For a base, the partitioning factor and protonation factor should be non-zero.
   !
   ! The terms as used as follows:
   !     H_eff = vth1 * (1 + vth2 / [H+] * (1 + vth3 / [H+])) * (1 + vth4 * [H+] / K_w)
   !     X_gas = X_total / (1 + H_eff * R * T * LWC)
   !      K_eq = X_gas * vth1 * (vth2 + vth4 / K_w) * P
   !           = X_total * vth1 * (vth2 + vth4 / K_w) * P / (1 + H_eff * R * T * LWC)
   !      [A-] = K_eq / [H+] (mono-protic acid)
   !     [HA-] = K_eq / [H+] (di-protic acid)
   !     [A--] = K_eq / [H+] * (1 + 2 * vth3 / [H+]) (di-protic acid)
   !     [BH+] = K_eq * [H+] (base)
   !
   ! where:
   !    H_eff is the effective Henry's Law constant (mol L-1 atm-1)
   !    X_gas is the gas-phase mixing ratio of the acid (mol mol-1)
   !    X_total is the total mixing ratio of the acid (mol mol-1)
   !    K_eq = [A-][H+] is the equilibrium constant (mol^2 L-2)
   !    P is the pressure (atm)
   !    R is the gas constant (L atm mol-1 K-1)
   !    T is the temperature (K)
   !    LWC is the cloud liquid water content (L_water L_air-1)
   !    [H+] is the hydrogen ion concentration (mol L-1)
   !    [A-] is the dissociated acid concentration (mol L-1)
   !    vth1 is the partitioning van 't Hoff factor
   !    vth2 is the first dissociation van 't Hoff factor
   !    vth3 is the second dissociation van 't Hoff factor
   !    vth4 is the protonation van 't Hoff factor
   !    K_w is the water dissociation constant (mol^2 L-2)
   !
   ! NOTE: If vth4 in non-zero, vth2 and vth3 should be zero.
   elemental subroutine henrys_law_set_conditions( this, i_column, &
         i_layer, temperature, pressure, cloud_water, total_mixing_ratio )

      class(henrys_law_t), intent(inout) :: this
      integer,             intent(in)    :: i_column
      integer,             intent(in)    :: i_layer
      real(r8),            intent(in)    :: temperature        ! K
      real(r8),            intent(in)    :: pressure           ! atm
      real(r8),            intent(in)    :: cloud_water        ! L_water L_air-1
      real(r8),            intent(in)    :: total_mixing_ratio ! mol mol-1

      real(r8) :: v1, v2, v3, v4, v4_kw, temp_delta

      temp_delta = 1._r8 / temperature - 1._r8 / this%reference_temperature_     ! K-1

      v1 = this%partitioning_factor_%A_ &
           * exp( this%partitioning_factor_%B_ * temp_delta )                    ! mol L-1 atm-1
      if (this%type_ == HENRYS_LAW_MONOPROTIC_ACID .or. &
          this%type_ == HENRYS_LAW_DIPROTIC_ACID .or. &
          this%type_ == HENRYS_LAW_NEUTRAL .or. &
          this%type_ == HENRYS_LAW_CO2) then
         v2 = this%first_dissociation_factor_%A_ &
              * exp( this%first_dissociation_factor_%B_ * temp_delta )           ! mol L-1
         if (this%type_ == HENRYS_LAW_DIPROTIC_ACID) then
            v3 = this%second_dissociation_factor_%A_ &
                 * exp( this%second_dissociation_factor_%B_ * temp_delta )       ! mol L-1
         else
            v3 = 0._r8
         end if
      else if (this%type_ == HENRYS_LAW_HO2) then
         v2 = this%ho2__hplus_o2minus_%A_ &
              * exp( this%ho2__hplus_o2minus_%B_ * temp_delta )                  ! mol L-1
         v3 = this%ho2_ho2__h2o2_o2_%A_ &
              * exp( this%ho2_ho2__h2o2_o2_%B_ * temp_delta )                    ! mol L-1
         v4 = this%ho2_o2minus__h2o2_o2_%A_ &
              * exp( this%ho2_o2minus__h2o2_o2_%B_ * temp_delta )                ! mol L-1
      else
         v2 = 0._r8
         v3 = 0._r8
      end if
      if (this%type_ == HENRYS_LAW_BASE) then
         v4_kw = this%protonation_factor_%A_ &
                 * exp( this%protonation_factor_%B_ * temp_delta ) / &
                 WATER_DISSOCIATION_CONSTANT                                     ! mol L-1
      else
         v4_kw = 0._r8
      end if
      if (this%type_ == HENRYS_LAW_HO2) then
         this%terms_(1,i_column,i_layer) = v1 * pressure                         ! mol_aq L-1 atm_air atm_gas-1
         this%terms_(2,i_column,i_layer) = v2                                    ! mol L-1
         this%terms_(3,i_column,i_layer) = v3                                    ! L mol-1 s-1
         this%terms_(4,i_column,i_layer) = v4                                    ! L mol-1 s-1
         this%terms_(5,i_column,i_layer) = cloud_water * &
                                           (pressure / PASCAL_TO_ATM) / &
                                           (GAS_CONSTANT * temperature)          ! m_h2o^3 molecule-1
      else
        this%terms_(1,i_column,i_layer) = v1 * (v2 + v4_kw) &
                                          * pressure * total_mixing_ratio        ! mol^2 L-2 (acid) or unitless (base)
        this%terms_(2,i_column,i_layer) = GAS_CONSTANT_L_ATM_MOL_K &
                                          * temperature * cloud_water            ! L mol-1 atm
        this%terms_(3,i_column,i_layer) = v1                                     ! mol L-1 atm-1
        this%terms_(4,i_column,i_layer) = v2                                     ! mol L-1
        this%terms_(5,i_column,i_layer) = v3                                     ! mol L-1
        this%terms_(6,i_column,i_layer) = v4_kw                                  ! mol-1 L
      end if

   end subroutine henrys_law_set_conditions

   !-------------------------------------------------------------------------------
   ! Returns the equilibrium constant for the current conditions
   !
   ! K_eq = [A-][H+] is the equilibrium constant (mol^2 L-2)
   ! See description for henrys_law_monoprotic_acid_set_conditions for details.
   elemental subroutine henrys_law_equilibrium_constant( this, &
         i_column, i_layer, h_plus_concentration, equilibrium_constant )
      
      class(henrys_law_t), intent(in)  :: this
      integer,             intent(in)  :: i_column
      integer,             intent(in)  :: i_layer
      real(r8),            intent(in)  :: h_plus_concentration ! mol L-1
      real(r8),            intent(out) :: equilibrium_constant ! mol^2 L-2

      real(r8) :: H_eff
      call this%effective_henrys_law_constant( i_column, i_layer, &
                                               h_plus_concentration, H_eff )
      equilibrium_constant = this%terms_(1,i_column,i_layer) &
                           / (1._r8 + H_eff * this%terms_(2,i_column,i_layer))

   end subroutine henrys_law_equilibrium_constant

   !-------------------------------------------------------------------------------
   ! Returns the gas-phase mixing ratio (mol mol-1) for the current conditions after
   ! partitioning.
   !
   ! See description for henrys_law_monoprotic_acid_set_conditions for details.
   elemental subroutine henrys_law_gas_phase_mixing_ratio( &
         this, i_column, i_layer, h_plus_concentration, total_mixing_ratio, &
         gas_phase_mixing_ratio )
      
      class(henrys_law_t), intent(in)  :: this
      integer,             intent(in)  :: i_column
      integer,             intent(in)  :: i_layer
      real(r8),            intent(in)  :: total_mixing_ratio     ! mol mol-1
      real(r8),            intent(in)  :: h_plus_concentration   ! mol L-1
      real(r8),            intent(out) :: gas_phase_mixing_ratio ! mol mol-1

      real(r8) :: H_eff  ! Effective Henry's Law constant (mol L-1 atm-1)

      if (this%type_ == HENRYS_LAW_HO2) then
          gas_phase_mixing_ratio = 0._r8 ! TODO: should be a failure condition
      else
        call this%effective_henrys_law_constant( i_column, i_layer, &
                                                 h_plus_concentration, H_eff )
        gas_phase_mixing_ratio = total_mixing_ratio &
                               / (1._r8 + H_eff * this%terms_(2,i_column,i_layer))
      end if

   end subroutine henrys_law_gas_phase_mixing_ratio

   !-------------------------------------------------------------------------------
   ! Returns the concentration of the neutral, unreacted aqueous-phase species
   ! concentration [mol(X) L(H2O)-1] after partitioning and dissociation/reactions.
   elemental subroutine henrys_law_partitioned_species_concentration( this, &
         i_column, i_layer, h_plus_concentration, total_mixing_ratio, &
         pressure, species_concentration )

      class(henrys_law_t), intent(in)  :: this
      integer,             intent(in)  :: i_column
      integer,             intent(in)  :: i_layer
      real(r8),            intent(in)  :: total_mixing_ratio     ! mol mol-1
      real(r8),            intent(in)  :: h_plus_concentration   ! mol L-1
      real(r8),            intent(out) :: pressure               ! atm
      real(r8),            intent(out) :: species_concentration  ! mol L-1

      real(r8) :: H_eff  ! Effective Henry's Law constant (mol L-1 atm-1)

      if (this%type_ == HENRYS_LAW_HO2) then
         species_concentration = 0._r8 ! TODO: should be a failure condition
      else
         call this%effective_henrys_law_constant( i_column, i_layer, &
                                                  h_plus_concentration, H_eff )
         species_concentration = total_mixing_ratio &                                ! mol(X) mol(air)-1 = atm(x) atm(air)-1
                               / (1._r8 + H_eff * this%terms_(2,i_column,i_layer)) & ! unitless
                               * this%terms_(3,i_column,i_layer) &                   ! mol(X) L(H2O)-1 atm(X)-1
                               * pressure                                            ! atm(air)
      end if

   end subroutine henrys_law_partitioned_species_concentration
   
   !-------------------------------------------------------------------------------
   ! Returns the rate of change in condensed-phase mixing ratio (mol mol-1 s-1)
   ! for the current conditions due to partitioning and condensed-phase reactions.
   elemental subroutine henrys_law_reaction_rate( this, i_column, i_layer, &
       total_mixing_ratio, h_plus_concentration, reaction_rate )

      class(henrys_law_t), intent(in)  :: this
      integer,             intent(in)  :: i_column
      integer,             intent(in)  :: i_layer
      real(r8),            intent(in)  :: total_mixing_ratio     ! mol mol-1
      real(r8),            intent(in)  :: h_plus_concentration   ! mol L-1
      real(r8),            intent(out) :: reaction_rate          ! mol mol-1 s-1

      real(r8) :: kh1_hplus ! KH1 / [H+] (unitless)
      real(r8) :: ho2_ss ! Steady state HO2 concentration (mol L-1)

      if (this%type_ == HENRYS_LAW_HO2) then
         kh1_hplus = this%terms_(2,i_column,i_layer) / h_plus_concentration
         ho2_ss = this%terms_(1, i_column, i_layer) * total_mixing_ratio &
                   * (1._r8 + kh1_hplus)
         reaction_rate = (this%terms_(3,i_column,i_layer) + &
                          this%terms_(4,i_column,i_layer) * kh1_hplus) &
                         / ((1._r8 + kh1_hplus)**2) * ho2_ss**2 &
                         * this%terms_(5,i_column,i_layer)
      else
         reaction_rate = 0._r8 ! This should be a failure condition
      end if

   end subroutine henrys_law_reaction_rate

   !-------------------------------------------------------------------------------
   ! Returns the Effective Henry's Law constant for the current conditions
   !
   ! See description for henrys_law_acid_set_conditions for details.
   elemental subroutine henrys_law_effective_henrys_law_constant( this, &
         i_column, i_layer, h_plus_concentration, effective_constant )
      
      class(henrys_law_t), intent(in)  :: this
      integer,             intent(in)  :: i_column
      integer,             intent(in)  :: i_layer
      real(r8),            intent(in)  :: h_plus_concentration   ! mol L-1
      real(r8),            intent(out) :: effective_constant     ! mol L-1 atm-1

      if (this%type_ == HENRYS_LAW_NEUTRAL .or. &
          this%type_ == HENRYS_LAW_MONOPROTIC_ACID .or. &
          this%type_ == HENRYS_LAW_DIPROTIC_ACID) then
         effective_constant = this%terms_(3,i_column,i_layer) &
                         * (1._r8 + this%terms_(4,i_column,i_layer) &
                            / h_plus_concentration &
                            * (1._r8 + this%terms_(5,i_column,i_layer) &
                               / h_plus_concentration))
      else if (this%type_ == HENRYS_LAW_BASE) then
         effective_constant = this%terms_(3,i_column,i_layer) &
                         * (1._r8 + this%terms_(6,i_column,i_layer) &
                            * h_plus_concentration)
      else
         effective_constant = 0._r8
      end if

   end subroutine henrys_law_effective_henrys_law_constant

end module cloud_aqueous_chemistry

