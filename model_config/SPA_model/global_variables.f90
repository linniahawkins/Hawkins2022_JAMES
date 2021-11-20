! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module gv_clim

  implicit none

  ! met-drivers used/updated each timestep..
  double precision :: atmos_press, & ! Surface atmospheric pressure (Pa)
          coa,         & ! Ambient CO2 concentration (ppm)
          par_top,     & ! PAR at top canopy layer
          snowfall,    & ! snowfall (kg.m-2.s-1)
          ppt,         & ! Precipitation (mm)
          sw_rad,      & ! Incident short wave radiation (Wm-2)
          lw_rad,      & ! Incoming longwave (W m-2)
          sw_diffuse,  & ! fraction of incoming shortwave that is diffuse
          temp_bot,    & ! Surface temperature at bottom and..
          temp_top,    & !  ..at top of canopy
          vpd_bot,     & ! Vapour Pressure Deficit at bottom..
          vpd_top,     & !  ..and at top of canopy.
          wind_spd,    & ! Wind speed (m/s)
          swc10,       & ! soil water content 10cm
          swc20,       & ! soil water content 20cm
          swc30,       & ! 
          swc40,       & !
          swc50,       & !
          swc60,       & !
          swc70,       & !
          swc80,       & !
          swc90,       & !
          swc100,      & !
          swc110,      & !
          swc120,      & !
          swc130,      & !
          swc140,      & !
          swc150         ! soil water content 150cm

  ! other meteorology-related variables.. 
  double precision    :: avtemp,      & ! average daily temperature (Celcius); used in deciduous
             lambda_bulk, & ! bulk surface latent heat of vapouration (J.kg-1)
             daypar,      & ! accumulation of PAR over a day (umol.m-2.s-1)
             dayppt,      & ! accumulation of precipitation over a day (mm.d-1)
             rnet,        & ! net radiation (W.m-2)
             wdbot,       & ! absolute water deficit at bottom of canopy (g m-3)
             wdtop          ! absolute water deficit at top of canopy (g m-3)

  double precision,parameter :: ppfd_to_par = 4.6d0  ! PAR energy convertion ratio with SW radiation (umol.J-1)
  double precision,parameter :: cp_air = 1004.6d0 ! Specific heat capacity of air; used in energy balance J.kg-1.K-1
  double precision,parameter :: gas_constant_d = 287.04d0 ! gas constant for dry air (J.K-1.mol-1)

  double precision,dimension(:),allocatable :: &
                    gbw, & ! boundary layer conductance for each canopy layer (m.s-1)
  leaf_heat_conductance, & ! boundary layer conductance for heat for each canopy layer
                  wetev    ! Wet evaporation from canopy (mm.t-1) per time step

end module gv_clim
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module  gv_daily_averages

  implicit none

  ! special variables needed for ACM recalibration runs are placed in here

  double precision, dimension(:), allocatable :: daily_ground_heat_flux &
                                                ,daily_canopy_temperature & 
                                                ,daily_shade_canopy_temperature &
                                                ,daily_sun_canopy_temperature &
                                                ,daily_weighted_SWP &
                                                ,daily_weighted_soilR &
                                                ,daily_soil_conductance  

end module  gv_daily_averages

!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_hourscale

  implicit none

  integer :: hour            ! ?
  double precision :: canopy_store = 0d0, & ! water retained on vegetation surface, mm
                discharge, & ! drainage from soil profile to groundwater, mm per timestep
               evap_store, & ! water evaporated from vegetation surface, mm per timestep
                      gaw, & ! boundary layer conductance (m.s-1), calculated each time step 
                      gws, & ! soil conductance to water vapour (m.s-1)
                  hourppt, & ! Met -> precip (mm per timestep)
                hourpress, & ! Met -> surface pressure (Pa)
                  hourrad, & ! Soil net radiation balance including both long and shortwave (W.m-2),
                             !  used in long wave determination.
                 hourrnet, & ! Met -> sw rad (W.m-2)
                 hourtemp, & ! Met -> temperature (oC) converted to (K)
                 hourtime, & ! Decimal time for reference
                   hourts, & ! soil surface temperature (K)
                  hourvpd, & ! Met -> VPD (kPa)
                 hourwind, & ! Met -> wind speed (m.s-1)
                 overflow, & ! over flow from water input into the soil layer. Set as a fixed 
                             !  proportion of surface_watermm (mm)
                 Qc = 0d0, & ! Soil heat flux (W.m-2) +ve heat moving upwards
                 Qe = 0d0, & ! Latent energy flux (soil) within the model structure (W.m-2); -ve evaporation occuring
                 Qh = 0d0, & ! Sensible heat flux from soil (W.m-2); -ve heat being released
                 Qn = 0d0, & ! Net LW emissions (W.m-2) - based upon the emissivity and the total radation penetration
                 Qm = 0d0, & ! Snowmelt heat flux (W.m-2) - same sign convention as Qe such that melting is negative
                 Qs = 0d0, & ! Latent energy flux from sublimation of snow (W.m-2)
                   runoff, & ! overland flow from non-infiltrated precip, mm per timestep
          surface_watermm, & ! timestep calculated surface water. i.e. the water that cannot be
                             !  infiltrated within the given timestep or surface layer saturation (mm) 
                    totet, & ! total evapotranspiration, m per timestep
                underflow, & ! water draining from soil profile to groundwater, mm per timestesp
            unintercepted    ! water (liquid) at soil surface, mm

  double precision, parameter :: freeze = 273.15d0 ! freezing point of water in Kelvin

end module gv_hourscale
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_Hydrol

  implicit none

  double precision,dimension(:),allocatable :: &
                soil_frac_clay, & ! Percentage of soil that is clay.
                soil_frac_sand    ! Percentage of soil that is sand.

end module gv_Hydrol
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_Irradiance_Sunshade

  implicit none

  double precision :: daylength

  double precision :: check,  & !
         skyabsl, & ! Total sky absorbance for long wave radiation (W.m-2); i.e. lw emitted back into the sky
          soilnet   ! absorbed radiation by soil (W.m-2); no long emissions, isothermal or otherwise

end module gv_Irradiance_Sunshade
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_metab

  implicit none

  double precision :: an, & ! GPP for given canopy layer in the shade or light leaf area loops (umolC.m-2.timestep-1)
          ci, & ! intra leaf CO2 concentration (ppm)
          ht, & ! Height of specific canopy layer, used in canopy loop (m)
 layer_capac, & ! Canopy layer specific canopy capacitance, based on canopy area
          rn, & ! Respiration constant for leaf area under photosynthetic analysis (umol CO2.gN.m-2 at 10oC)
      rplant, & ! Plant hydrolic resistance for each canopy layer (MPa.s-1.mmol-1)
       rsoil, & ! Root and soil resistence for a given canopoy layer (MPa.s-1.mmol-1)
         vcm, & ! Maximum rate of carboxylation (umol CO2.gN.m-2.s-1); based on kappaV coefficient
         vjm    ! Maximum rate of electrob transport (umol.m-2.s-1)

  ! metabolic temperature optimum (default = 30C, arctic = ~13C, tropics = 40C)
  double precision, parameter :: metabolic_opt_temp = 30d0, & ! metabolic temperature optimum
                                    vcmax_max_temp = 65.03d0, & ! max tolerated temperature for carboxylation
                                     jmax_max_temp = 57.05d0    ! max tolerated temperature for electron transport

end module gv_metab
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_meteo

  implicit none

  double precision ::   gbb, & ! boundary layer conductance of specific canopy layer in spa_canopy.F loop (m3.s-1)
            gbh, & ! boundary layer conductance for heat for a given canopy layer (m.s-1)
         roughl, & ! roughness length of the land surface (m)
       dew_evap, & ! dew or canopy evaporation for whole canopy (W.m-2)
      live_frac, & ! fraction of living canopy (crops only)
             la, & ! Leaf area in current photosynthetic pass, i.e. sun or shade leaf areas (m.canopy_layer-1))
            nit, & ! Canopy layer specific nitrogen content (gN.m-2 leaf area)
      layer_LCA, & ! Layer specific lca value used in leaf.f90 (gC.m-2)
            par, & ! PAR penetrating given canopy layer in spa_canopy.F loop (umol.m-2.s-1)
           psil, & ! canopy layer specific leaf water potential in spa_canopy.F loop (MPa)
         lwp_pd, & ! predawn leaf water potential in current canopy layer (MPa)
           psis, & ! Weighted soil water potential, by total evaporation estimate (MPa)
            rad, & ! net radiation penetration to canopy layer (kW.m-2)
           temp, & ! temperature of a given canopy layer (oC); see spa_canopy.F
           wdef, & ! air water content deficit (g m-3)
   displacement, & ! zero plane displacement height of the land surface (m)
 air_density_kg, & ! density of air (kg m-3)
   abs_pot_conv, & ! convertion of temperature from absolute to potential value (K -> K)
  abs_virt_conv    ! convertion of temperature from absolute to virtual value (K -> K)

  integer :: rad_pass ! looping variable for iteration of the canopy
  double precision,parameter :: gi = 1d0, & ! mesophyll conductance (m s-1) set high, and thus ignored, /0.0025/ old value
                                            ! v large mesophyll conductance - ACi curves have not been Epron adjusted
                       head = 0.009807d0, & ! head of pressure  (MPa/m)
                       Rcon = 8.3144d0      ! universal gas constant (J.mol-1.K-1)

end module gv_meteo
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_scale_declarations

  !! These declarations define the size of the soil  !!
  !!  and canopy profiles, & the time resolution.    !!
  !! They also declare constructs that are needed at !!
  !!  all levels, such as the user's config options, !!
  !!  the time holder, and the met-drivers.          !!
  !! Physical parameters, such as pi, are also here. !!

  implicit none

  ! Internal model lengths (fixed)..
  integer,parameter :: fname_length = 200, &  ! length of filename variables
                       max_nos_iterations = 2 ! number of iterations in math-loops

  ! Model grid sizing info
  type grid_holder
    integer :: canopy    = 10,    & ! number of canopy layers
               core      = 21,    & ! number of soil layers + 1
               soil      = 20,    & ! number of soil layers
               snow      =  5,    & ! number of snow layers 
               wetting   = 10       ! number of layers to use for wetting calcs
    double precision    :: latitude  = 50d0, & ! +ve == nth, -ve == sth
                           longitude = 0d0     ! 0--360. -ve not recognised.
  end type
  type(grid_holder),save :: grid

  ! All the time information in one place..
  type time_holder
    character(len=8) :: time_of_day ! purely diagnostic, just tells us what the time-of-the-day is.
    character(128) :: start_date = ""      ! starte date of simulation (fmt = '2002-01-01_00:00:00')
    integer, dimension(:), allocatable :: days_in_year ! number of whole days in each year simulated
    integer :: nos_of_years = 1,    & ! number of years to simulate
               run_day = 0,         & ! current day of the analysis
               steps_per_day = 24,  & ! number of timesteps per day
               year = 0,            & ! Current year
               day  = 0,            & ! Current day
               step = 0,            & ! Current step
               steps_count = 0              ! count of number of steps completed so far.
    double precision :: seconds_per_step = 3600d0, & ! timesteps (3600=>fixed at 24 steps per day)
                         seconds_per_day = 86400d0, &
                                 daytime = 0d0       ! date+time in days, e.g. noon on day 30 of yr 2
                                                     ! (if yr=365days) == 365 + 30 + 0.5 = 395.5
  end type time_holder
  type(time_holder),save :: time


  ! meteorological inputs/drivers
  type met_drivers
    double precision,allocatable,dimension(:) :: ambient_co2  ! (ppm) Ambient Carbon Dioxide concentration
    double precision,allocatable,dimension(:) :: ppfd         ! (umol.m-2.s-1) Photosynthetically active radiation
    double precision,allocatable,dimension(:) :: precip       ! (kg.m-2.s-1 or mm.t-1) precipitation
    double precision,allocatable,dimension(:) :: snowfall     ! (kg.m-2.s-1 or mm.t-1) snowfall
    double precision,allocatable,dimension(:) :: sfc_pressure ! (Pa)  Atmospheric surface pressure
    double precision,allocatable,dimension(:) :: sw_rad       ! (W.m-2) surface downward short-wave radiation
    double precision,allocatable,dimension(:) :: sw_diffuse   ! (W.m-2) Diffuse short wave radiation
    double precision,allocatable,dimension(:) :: lw_rad       ! (W.m-2) long wave incoming
    double precision,allocatable,dimension(:) :: temp         ! (C) temperature
    double precision,allocatable,dimension(:) :: vpd          ! (kPa) vapour pressure deficit
    double precision,allocatable,dimension(:) :: wind_spd     ! (m.s-1) wind strength
    double precision,allocatable,dimension(:) :: clearance_fraction   ! fraction of C pools removed by mechanical disturbance 
    double precision,allocatable,dimension(:) :: fire_fraction! fraction of C pools removed by fire disturbance
    double precision,allocatable,dimension(:) :: swc10 ! soil water content 10cm
    double precision,allocatable,dimension(:) :: swc20 ! soil water content 20cm
    double precision,allocatable,dimension(:) :: swc30 !
    double precision,allocatable,dimension(:) :: swc40 !
    double precision,allocatable,dimension(:) :: swc50 !
    double precision,allocatable,dimension(:) :: swc60 !
    double precision,allocatable,dimension(:) :: swc70 !
    double precision,allocatable,dimension(:) :: swc80 !
    double precision,allocatable,dimension(:) :: swc90 !
    double precision,allocatable,dimension(:) :: swc100 !
    double precision,allocatable,dimension(:) :: swc110 !
    double precision,allocatable,dimension(:) :: swc120 !
    double precision,allocatable,dimension(:) :: swc130 !
    double precision,allocatable,dimension(:) :: swc140 !
    double precision,allocatable,dimension(:) :: swc150 !
    integer,allocatable,dimension(:) :: clearance_management ! post clearance management type (see carbon_model.f90)
    ! If user does not provide variables, but still wants to alter co2/pressure from defaults
    ! then they can supply these in the 'met parameters' section of the config file..
    double precision :: const_ambient_co2  = 380d0    ! (ppm) Ambient Carbon Dioxide concentration
    double precision :: const_sfc_pressure = 101325d0 ! (Pa)  Atmospheric surface pressure
  end type met_drivers
  type(met_drivers),save :: met

  ! phenological inputs/drivers
  type phenology_drivers
    double precision,allocatable,dimension(:) :: lai   ! leaf area index from phenological drivers (m2/m2)
    double precision,allocatable,dimension(:) :: rootC ! root carbon input (gC.m-2)
  end type phenology_drivers
  type(phenology_drivers),save :: pheno

  ! All of the configuration info for SPA...
  type user_config_holder
     ! Input files --------------------------------------------------------
     character(fname_length) :: crops_filename     = ''  ! Crops file name/path
     character(fname_length) :: met_filename       = ''  ! Met file name/path
     logical                 :: met_file_is_nc = .false. ! Is met file in netcdf format?
     character(fname_length) :: phenology_filename = ''  ! phenology file name/path
     character(fname_length) :: soils_filename     = ''  ! Soils file name/path
     character(fname_length) :: veg_filename       = ''  ! Veg file name/path
     character(fname_length) :: carbon_filename    = ''  ! Carbon file name/path
     character(fname_length) :: restart_filename   = ''  ! file to use for restart
     logical                 :: load_restart = .false.   ! Should SPA start from a restart dump?
     ! Options ----------------------------------------------------------
     integer :: soil_hydrology_opts      = 1       ! Soil water potential parameterisation: Saxton = 1, Not current assigned = 2
     logical :: prescribed_swc           = .false.       ! false=prognostic soil water content, true=prescribe SWC in met forcing file
                                                   ! (Hawkins et al., )
     integer :: plant_func_type          = 1       ! Determines whether doing a deciduous/evergreen(1) & 
                                                   ! or arable crops(2), tubers (3)
     integer :: canopy_MDecay            = 0       ! canopy momentum decay (1 = Williams et al 1996, 2 = Smallman et al 2013)
     logical :: include_leap_years       = .false. ! should be include leap years, needed for long real data scenarios
     logical :: loop_met_driver          = .false. ! repeat the met-driver as needed?
     logical :: loop_pheno_driver        = .false. ! repeat phenology driver as needed
     integer :: iWUE                     = 1       ! use iWUE(1) or WUE(2) or Ball-Berry(3) or Medlyn (4) for gs optimisation
     logical :: dead_lai_energy_balance  = .false. ! include dead lai in the surface energy balance and 
                                                   ! radiative transfer (crop only)
     logical :: iterative_canopy         = .true.  ! use iterate long wave and canopy energy balance
     logical :: met_file_has_lw_rad      = .false. ! Does met file contain..Incoming Longwave (LIN) ?
     integer :: solve_canopy_temperature = 2       ! choice of how to solve canopy energy balance 
                                                   ! 1 = default steady state, 2 = uses an iterative approach
     logical :: temp_in_kelvin           = .false. ! whethe we are using kelvin or not for temperature input
     logical :: use_co2_from_met_file    = .false. ! Whether to use the CO2 values in the input met file
     logical :: use_ppfd_from_met_file   = .false. ! likewise for PAR.
     logical :: met_file_has_sfc_press   = .false. ! likewise for surface atmospheric pressure
     logical :: met_file_has_sw_diffuse  = .false. ! whether file has diffuse short wave radiation
     logical :: met_file_has_snowfall    = .false. ! does met file contain snowfall information
     logical :: met_file_has_disturbance = .false. ! does met file contain clearance or fire information
     logical :: precip_is_rate           = .false. ! Is precip in rate (mm/sec) or volume-per-timestep (mm)?
     logical :: use_dalec                = .true.  ! switch to use dalec phenology or not
     character(fname_length) :: vpd_or_specific_humidity = '' ! are we using vpd (kPa) or specific humidity (kg/kg)
     ! Output ----------------------------------------------------------
     integer :: print_msg_at_severity = 5 ! print-to-screen messages of this severity or higher
     character(fname_length) :: output_directory =  ''    ! directory to put output data in
     logical :: std_csv_output          = .true.  ! Set if you want to do calcs for deciduous veg
     logical :: Ecofluxes_csv_output    = .true.  ! Set if you want to do calcs for deciduous veg
     logical :: canopy_csv_output       = .true.  ! Set if you want to do calcs for deciduous veg
     logical :: Cstock_csv_output       = .true.  ! Set if you want to do calcs for deciduous veg
     logical :: soils_csv_output        = .true.  ! Set if you want to do calcs for deciduous veg
     logical :: make_restart_file       = .false. ! Should SPA make a restart dump?
     integer :: restart_write_frequency = 1000    ! frequency (in years) with which to produce restart dumps
  end type user_config_holder
  type(user_config_holder),save :: user_opts

  ! values used for minimum and maximum boundings
  double precision, parameter :: dble_two = 2d0, &
                                 dble_one = 1d0, & 
                                 dble_zero = 0d0

  ! General physical constants..
  double precision,parameter :: boltz = 5.670400d-8,    & ! (W m-2 K-4)
                                boltz_kW = boltz * 1d-3,& ! (kW.m-2.K-4)
                                pi = 3.14159265d0,  & ! (-)
                                deg_to_rad = pi / 180d0 ! ratio for converting degrees to radians

  ! molar conversion coefficients
  double precision, parameter :: &
                density_of_water = 998.2d0,    & ! density of water kg.m-3 
                 mol_to_g_carbon = 12d0,       & ! molecular mass of carbon (g)
                  mol_to_g_water = 18d0,       & ! molecular mass of water
                 g_to_mol_carbon = 1d0 / 12d0, & !
                  g_to_mol_water = 1d0 / 18d0, & ! molar ration for H2O = 18 g 
                umol_to_g_carbon = 1d-6 * mol_to_g_carbon, & ! molecular mass of carbon (g)
                 umol_to_g_water = 1d-6 * mol_to_g_water,  &
                g_to_umol_carbon = 1d6 * g_to_mol_carbon,  & !
                 g_to_umol_water = 1d6 * g_to_mol_water      ! molar ration for H2O = 18 g


  ! general timing constants      
  integer,dimension(12), parameter :: days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  double precision, parameter :: avg_days_in_year = 365.25

  ! Internal store of what version of the model this is,
  ! where (23,0) equates to "2.3.0"
  integer,parameter :: model_sci_version = 30, & ! 20=>2.0, 23=> 2.3, etc.
                       model_bug_version = 0     ! 0=>0, 1=>1, etc.  


end module gv_scale_declarations
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_soil_structure

  implicit none

  double precision, parameter :: abovebelow = 1d0 ! saxton water retention equation are off by default

  integer :: rooted_layers     ! number of root layers penetrated into the soil layers

  double precision :: drythick = 0.1d0, & ! Thickness of dry soil layer above water table (m); minimum level is set to 0.001 (m)
                             prevwater, & ! stores the initial and then previous total water content used in updating water balance
                             max_depth, & ! PARAMETER maximum rooting depth (m)
                           max_storage, & ! PARAMETER maximum canopy water storage (mm)
                  resp_rate_temp_coeff, & ! response of autotrophs to change in temperature
                          root_biomass, & ! root biomass (g Biomass.m-2) in total is determined as twice the root C (gC.m-2)
                            root_reach, & ! maximum depth of root system based on the available biomass. 
                                          ! i.e. does not have to be the max_rooting depth parameter
                                root_k, & ! PARAMETER mass of roots for reaching 50% maximum depth (g.m-2!)
                                  snow, & ! Used in calculation of change in soil profile wetting space
                          surf_biomass, & ! calculation of surface biomass (g), 
                                          ! assumes that 50 % of biomass is in top 25 % of the soil profile
                          through_fall, & ! PARAMETER fraction of precip as throughfall 
                           k_ice = 2.2, & ! thermal conductivity of ice, W/ m K
                          weighted_SWP    ! Soil water potential for whole soil profile, &
                                          ! weighted by the estimate of total evaporation. 
                                          ! This value is determined each time step for use in spa_canopy

  double precision :: thermal = 1.58d0    ! conversion rate of thermal conductivity from (W m-1 K-1) to (J m-1 K-1 h-1)
  double precision,dimension(4) :: swp_params=-0.4d0 ! swp params
  ! Key components of the Saxton soil water retention equations, these are re calculated on each timestep
  double precision,dimension(:),save,allocatable :: cond1, cond2, cond3, potA, potB

  double precision,dimension(:),allocatable :: &
                          conduc,  & ! Soil layer conductivity (m.s-1)
                 field_capacity,   & ! Field capacity of moisture for each layer (mm?), when soil water content at SWP = -10kPa
                fraction_uptake,   & ! fraction of evapotranspiration (i.e. root draw) from each layer. Reset at each timestep
                        iceprop,   & ! Soil ice proportion 
                    layer_depth,   & ! PARAMETER Soil layer depth (m)
                    mineralfrac,   & ! PARAMETER soil proportional inorganic content
                    organicfrac,   & ! PARAMETER soil proportional organic content
                       porosity,   & ! soil layer porosity
                        pptgain,   & ! soil water potential (MPa); used in determining vapour pressure (kPa) of soil air space
                                     !   for latent heat fluxes and water movement between between layers.
                    root_length,   & ! (m) based on root biomass and an exponential decline of mass with soil depth.
                                     !   Assumes root subroutine called each timestep.
                      root_mass,   & ! (g biomass) per soil layer are determined every timestep in root call.
                      soil_temp,   & ! soil temperature (K)
               soil_temp_nplus1,   & ! temporary soil temperature (K) vector used for determining the next time step temperature
                                     !   soil profile beginning with the new surface temperature ure and the constant soil core value.
                          soilR,   & ! soil root hydraulics resistance (MPa.s-1.m-2.mmol-1); calculated from soilR1 and SoilR2 
                         soilR1,   & ! soil root hydraulics component; represents some function of the root mass and volume
                         soilR2,   & ! soil root hydraulics component; represents some function of root mass and the root resistance
                                     !   and ratio of mass above and below ground
                            SWP,   & ! Soil water potential (MPa); used in determining vapour pressure (kPa) of soil air space for
                                     !   latent heat fluxes and water movement between between layers.
                      thickness,   & ! PARAMETER Soil layer thickness (m)
                      waterfrac,   & ! Extracts soil moisture by volumetric mixing ratio (m3.m-3) from its array  
                      watergain,   & ! Water gained (mm) from a given soil layer; reset with each timestep
                      watericemm,  & ! Total soil moisture content per layer (mm)
                      waterloss,   & ! Water lost (mm) from a given soil layer; reset with each timestep
                      wettingbot,  & ! Depth to bottom of wet soil layers (m)
                      wettingtop     ! Depth to top of wet soil layers (m)

end module gv_soil_structure
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_snow_info

  implicit none

  integer :: Nsnow         ! Number of snow layers

  double precision, parameter :: new_snownirref = 0.73d0, & ! NIR reflectance of new snow
                                 new_snowparref = 0.95d0    ! PAR reflectance of new snow

  double precision :: snowalb_nir, & ! NIR albedo of snow
                      snowalb_par, & ! PAR albedo of snow
                     snow_watermm, & ! Snow water equivalent (mm)
                            fsnow, & ! fraction of snow cover  
                       snowheight, & ! Snow depth (m)
                       snowweight    !

  double precision,dimension(:),allocatable ::  &
                       Dsfix,       & ! Fixed snow layer thickness (m)
                       Dsnow,       & ! Snow layer thicknesses (m)
                       Sice,        & ! Snow layer ice mass (kg/m2)
                       Sliq,        & ! Snow layer liquid water mass (kg/m2)
                       Tsnow          ! Snow layer temperatures (K) 

end module gv_Snow_Info
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_veg

  implicit none

  logical,allocatable :: c3(:) ! whether to perform c3 or c4 calcs for each layer

  integer :: nlink        = 1, &  !
             conductivity = 1     ! 1=>conductivity set, 0=>conductance set

  ! as parameter for the moment
  double precision, parameter :: emiss = 0.96d0

  ! some parameters
  double precision, parameter :: root_radius = 0.00029d0, & ! root radius (m) Bonen et al 2014 = 0.00029
                                root_density = 0.31d6    ! root density (g biomass m-3 root) 
                                                         ! 0.5e6 Williams et al 1996                                       
                                                         ! 0.31e6 Bonan et al 2014

  double precision :: rootresist = 25d0                  ! fine root hydraulic resistivity (MPa s g mmol-1 H2O)

  integer :: canopy_level      ! level at which canopy base is
  double precision,dimension(2) ::  dimen = 0.08d0 ! characteristic leaf dimension (m); 
                                                   ! deciduous leaf dimen(1) & dimen (2) = leaf width; 
                                                   ! needle leaves require dimen(2) = cone diameter
  double precision :: altitude, & !
                       avN,  & ! average foliar N (g N m-2 leaf area) = (total N)/LAI
               stock_roots,  & ! roots persists here because of need to work in mode without DALEC
             canopy_height,  & ! canopy height (m) is land cover type specific, generated from top canopy layer
               canopy_base,  & ! height above ground (m) of base of canopy
                     capac,  & ! capacitance (mmol m-2 LA MPa-1)
                    co2amb,  & ! ambient atmospheric CO2 (umol.mol-1)
                    gplant,  & ! depends on conductivity switch. if 0 => plant hydraulic conductance  (mmol m-2 s-1 MPa-1),
                               !                                 if 1 => plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      iWUE,  & ! stomatal efficiency (umolCO2.mmol-1H2O)
                               ! either threshold GPP increase per increase in stomatal openning
                       WUE,  & ! stomatal efficiency parameter (umolCO2.mmol-1H2O).
                               ! Threshold for GPP increase for given increase
                               ! in transpiration (See Bonan et al., 2014)
                cond_slope,  & ! slope parameter in Ball-Berry or Medlyn
                               ! See Hawkins et al.,
              cond_slope_b,  & ! Slope parameter sensitivity to LWP b parameter 
                               ! (Hawkins et al., )
              cond_slope_c,  & ! Slope parameter sensitivity to LWP c parameter
                               ! (Hawkins et al., )
                    kappac,  & ! Vcmax-N relationship (umol (gN)-1 s-1)
                    kappaj,  & ! Vjmax-N relationship (umol (gN)-1 s-1)
          latitude_radians,  & ! 
                    minlwp,  & ! minimum leaf water potential (MPa)
                   modRnet,  & ! modelled net radiationn (W.m-2)
                       tsk,  & ! effective skin temperature of surface (K)
                canopy_rad,  & ! net canopy radiation (W.m-2)
               dew_to_soil,  & ! dew formed on canopy (mm or kg.m-2.t-1) which is assumed to immediately run onto the soil
                 dew = 0d0,  & ! dew formed on canopy (mm or kg.m-2.t-1); added to canopy in next step
               prevC = 0d0,  & ! stored variable for us in difference calculation
                     sensh,  & ! NOT USED
                    totass,  & ! NOT USED
                   totevap,  & ! total evapotranspiration (leaf) (mm.t-1)
                     totla,  & ! total leaf area index
                      totN,  & ! total foliar N
!                 nfrac_top,  & !
                    totres,  & ! NOT USED
              tower_height     !

  ! Arrays..
  double precision,dimension(:),allocatable ::  &
                canopy_soil_resistance, & ! belowground hydraulic resistance of a canopy layer (MPa m2 s mmol-1)
                leaf_temp_intermediate, & ! leaf temperature used in intermediate steps of canopy iteration (oC) 
                             leaf_temp, & ! leaf temperature final value (oC)
                    can_sensible_layer, & ! canopy layer specific sensible heat flux (W.m-2)
                                   ess, & ! soil surface evaporation (W.m-2)
                            netrad_day, & ! net ecosystem radiation balance (W.m-2)
                system_energy_balance,  & ! system energy balance residual (W.m-2)
                                          ! +ve energy not used, -ve energy
                                          ! created
                                  gppt, & ! photosynthesis, micromoles/m2/s
                                lafrac, & ! leaf area fraction in canopy layer
                           lafrac_dead, & ! dead leaf area fraction in canopy (crop only)
                                   LCA, & ! leaf C per leaf area (gC.m-2); 
                                          ! NOTE: local declaration of LCA in carbon_model_crop
                                   lai, & ! leaf area index m2/m2
                               psil_pd, & ! predawn leaf water potential (MPa)
                              dead_lai, & ! lai which is dead but still standing; used in crops
!                          Cfol_profile, & ! profile tracking of LCA in case of future varying of LCA
                          layer_height, & ! m
                              LWPstore, & ! store LWP
                           LWPprevious, &
                                 nfrac, & ! fraction of foliar N in canopy layer
                                   Nla, & ! average foliar N per canopy layer (gN.m-2)
                                 respt, & ! autotrophic respiration, micromoles/m2/s
                              soiletmm, & ! soil evaporation, mm.t-1
                                 lemod, & ! total latent energy (W.m-2)
                                 mmmod, & ! total water flux (mm.step)
                                neemod, & ! net ecosystem C (umol m-2 s-1)
                              timeflux, & ! water flux canpy total (mmol m-2 GA s-1)
                                 wetle, & ! wet surface latent energy (W.m-2)
                 potential_evaporation, & ! potential evaporation (kg.m-2.t-1)
                         sensible_heat, & ! sensible heat (W.m-2)
                           transt_kg_s, & ! transpiration kg.m-2.s-1
                                transt    ! transpiration, W.m-2

  double precision,allocatable :: flux(:,:) !

end module gv_veg
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
module gv_carbon_model
  
  implicit none

  ! management / disturbance drivers updated each time step
  ! see
  double precision :: clearance_fraction   & ! 
         ,fire_fraction        & !
         ,clearance_management   !

  ! parameters for crop and default DALEC models
  integer :: nopars,nofluxes,nopools
  integer, parameter :: nopars_default   = 38, & !
                        nofluxes_default = 20, & !
                        nopools_default  = 7,  & !
                        nopars_crops     = 35, & !
                        nofluxes_crops   = 17, & !
                        nopools_crops    = 8     !

  ! Carbon model parameter vector
  double precision, dimension(:), allocatable :: pars
  ! Carbon pools (gC.m-2)...
  double precision, dimension(:,:), allocatable :: POOLS ! array of carbon pools

  ! Fluxes to/from carbon pools..
  double precision :: leaf_growth = 0d0, leaf_death = 0d0 ! growth and death fluxes (gC.m-2.day-1)
  double precision, dimension(:,:), allocatable :: FLUXES ! array of carbon fluxes (gC.m-2.day-1)
  double precision, dimension(:), allocatable :: gpp_lai_marginal & ! marginal return of gpp due to increased lai (gC.m-2.t-1)
                                    ,gpp_lai_marginal_reference & ! 
                                    ,daily_fluxes       ! daily sum of C fluxes in gC.m-2

  double precision, dimension(:), allocatable :: avg_min_airt_store,avg_dayl_store,avg_vpd_Pa_store &
                                    ,avg_min_airt,avg_dayl,avg_vpd_Pa ! GSI drivers averaged over 21 day period
                                                                      ! oC,seconds,Pa

  double precision :: GPP = 0d0, & ! carbon fluxes per given step (gC.m-2.day-1)
                      NEE = 0d0, & 
                resp_auto = 0d0, & ! Respiration - autotrophic
            resp_h_litter = 0d0, & ! Respiration - heterotrophic, in litter pool
     resp_h_soilOrgMatter = 0d0    ! Respiration - heterotrophic, in soil organic matter pool

  ! Reich model for maintenace respiration 
  double precision :: twq ! mean air temperature of warmest quarter (oC)

  ! Extra for Crops...
  double precision :: dleaf  !, SCA

  ! standard io names for carbon model. Note that I have not found a way of
  ! assigning the pool, par and flux descriptions that does not involve the
  ! character lengths having to be the same for a simultaneous assignment. Will
  ! try and fix this later
  character(200),dimension(:), allocatable :: fluxes_names,pools_names,pars_names
  character(200),parameter :: pars_names_default(nopars_default) = (/"lit->som_coef_day  ","GPP->Ra_fraction   " &
                                                                    ,"GSI_leaf_growth    ","GPP->root_coef     " &
                                                                    ,"leaf_turnover_day  ","wood_turnover_day  " &
                                                                    ,"root_turnover_day  ","Rhet_lit_coef_day  " &
                                                                    ,"Rhet_som_coef_day  ","resp_rate_coef(T)  " &
                                                                    ,"avgN_log10(gN_m2)  ","lab_turnover_day   " &
                                                                    ,"GPP->lab_coef      ","GSI_min_temp_coef  " &
                                                                    ,"GSI_max_temp_coef  ","GSI_min_dayl_coef  " &
                                                                    ,"LCA                ","initial_labile     " &
                                                                    ,"initial_foliage    ","initial_root       " &
                                                                    ,"initial_wood       ","initial_litter     " &
                                                                    ,"initial_som        ","GSI_max_dayl_coef  " &
                                                                    ,"GSI_min_vpd_coef   ","GSI_max_vpd_coef   " &
                                                                    ,"lai_gpp_return_coef","wood_is_branch_coef" &
                                                                    ,"wood_is_coarse_coef","Replant_labile     " &
                                                                    ,"Replant_foliage    ","Replant_root       " & 
                                                                    ,"Replant_wood       ","GSI_leaf_senescence" &
                                                                    ,"GSI_status_<_1_>   ","initial_GSI_(1-2)  " &
                                                                    ,"initial_cwd        ","cwd_turnover_day   "/) &
                             ,pars_names_crops(nopars_crops) = (/"lit->som_coef_day  ","GPP->Ra_fraction   " &
                                                                ,"max_DR_(DS_0-1)_day","max_DR_(DS_1-2)_day" &
                                                                ,"leaf_turnover_day  ","wood_turnover_day  " &
                                                                ,"leaf_shade_turnover","Vernalisation_days " &
                                                                ,"Rhet_lit_coef_day  ","Rhet_som_coef_day  " &
                                                                ,"avgN_log10(gN_m2)  ","sowing_day         " &
                                                                ,"resp_cost_lab->NPP ","PHU_for_emergence  " &
                                                                ,"harvest_day        ","plough_day         " &
                                                                ,"LCA                ","initial_labile     " &
                                                                ,"initial_foliage    ","initial_root       " &
                                                                ,"initial_wood       ","initial_litter     " &
                                                                ,"initial_som        ","initial_auto       " &
                                                                ,"initial_storage    ","minT_development   " &
                                                                ,"maxT_development   ","optT_development   " &
                                                                ,"minT_vernalisation ","maxT_vernalisation " &
                                                                ,"optT_vernalisation ","daylength_critial  " &
                                                                ,"daylength_coef     ","lab_turnover_day   " &
                                                                ,"auto_turnover_day  "/) &
                             ,pools_names_default(nopools_default) = (/"labile_gC.m-2 ","foliage_gC.m-2" &
                                                                      ,"roots_gC.m-2  ","wood_gC.m-2   " &
                                                                      ,"litter_gC.m-2 ","som_gC.m-2    " &
                                                                      ,"cwd_gC.m-2    "/) &
                             ,pools_names_crops(nopools_crops) = (/"labile_gC.m-2 ","foliage_gC.m-2" &
                                                                  ,"roots_gC.m-2  ","wood_gC.m-2   " &
                                                                  ,"litter_gC.m-2 ","som_gC.m-2    " &
                                                                  ,"rauto_gC.m-2  ","storage_gC.m-2"/)&
                             ,fluxes_names_default(nofluxes_default) = (/"GPP_gC.m-2.day-1      ","resp_rate(T)          " &
                                                                        ,"Rauto_gC.m-2.day-1    ","NPP->leaf_gC.m-2.day-1" &
                                                                        ,"NPP->lab_gC.m-2.day-1 ","NPP->root_gC.m-2.day-1" &
                                                                        ,"NPP->wood_gC.m-2.day-1","lab->leaf_gC.m-2.day-1" &
                                                                        ,"leaf->lit_fraction    ","leaf->lit_gC.m-2.day-1" &
                                                                        ,"wood->cwd_gC.m-2.day-1","root->lit_gC.m-2.day-1" &
                                                                        ,"Rhet_lit_gC.m-2.day-1 ","Rhet_som_gC.m-2.day-1 " &
                                                                        ,"lit->som_gC.m-2.day-1 ","lab->leaf_fraction    " &
                                                                        ,"Fire_gC.m-2.day-1     ","GSI_(0-1)             " &
                                                                        ,"cwd->lit_gC.m-2.day-1 ","Harvest_C_gC.m-2.day-1"/) &
                             ,fluxes_names_crops(nofluxes_crops) = (/"GPP_gC.m-2.day-1      ","resp_rate(T)          " &
                                                                    ,"Rauto_gC.m-2.day-1    ","NPP->leaf_gC.m-2.day-1" &
                                                                    ,"NPP->lab_gC.m-2.day-1 ","NPP->root_gC.m-2.day-1" &
                                                                    ,"NPP->wood_gC.m-2.day-1","lab->NPP_gC.m-2.day-1 " &
                                                                    ,"NPP->stor_gC.m-2.day-1","leaf->_gC.m-2.day-1   " &
                                                                    ,"wood->_gC.m-2.day-1   ","root->lit_gC.m-2.day-1" &
                                                                    ,"Rhet_lit_gC.m-2.day-1 ","Rhet_som_gC.m-2.day-1 " &
                                                                    ,"lit->som_gC.m-2.day-1 ","GPP->auto_gC.m-2.day-1" &
                                                                    ,"Harvest_C_gC.m-2.day-1"/)

  save

end module gv_carbon_model
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
