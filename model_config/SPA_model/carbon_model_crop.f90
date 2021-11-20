! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module carbon_model_crop_mod

  use gv_veg, only: stock_roots,dead_lai

  !!                                                      !!
  !! This module provides procedures to model the effect  !!
  !! of cropland management.                              !!
  !! It was developed by Oliver Sus' in May 2010.         !!
  !! Subsequent bug corrections and integration into SPA  !!
  !! by Luke Smallman January 2017                        !!
  !!                                                      !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: carbon_model_crop,stock_dead_foliage

  ! variables shared across spa..
  public :: BM_EX, Cshoot, DR, DR_P, DR_pre, DR_post, DR_T_PRA, DR_T_POA, DRAO_P, DRAO_T_POA, DRAO_T_PRA, &
            DS, DS_LRLV, DS_LRRT, DS_root, DS_shoot, emerged, fol_frac, fol_frac_intpol, FP, FT, FV, harvest_day, &
            LAICR, LCA_DS, LCA_ratio, LRLV, LRRT, lv_res, HI, PHCR, PHSC, PHUem, plough_day, remob, RDRSHMAX, root_frac, &
            root_frac_intpol, shoot_frac_intpol, stem_frac_intpol, sow_day, sown, stem_frac, st_res,    &
            tmax, tmax_v, tmin, tmin_v, topt, topt_v, use_seed_labile, VD, VDh, vernal_calcs, yield, &
            stock_seed_labile, raso, max_raso, stock_surface_litter

  ! variables local to this module..
  integer :: DLCOMP, grindex, LCADYN, nseasons, rank_DR_ratio, rank_DR_temp, y, year

  integer :: plough_day  = 170      ! day-of-year when field is ploughed  (default)
  integer :: sow_day     = 180      ! day-of-year when field is sown      (default)
  integer :: harvest_day = 300      ! day-of-year when field is harvested (default)
  ! crop management flags
  logical :: ploughed, vernal_calcs = .true., sown, use_seed_labile, emerged

  double precision,dimension(:),allocatable :: DR_P, DR_photo, DR_ratio, DR_ratio_photo, DR_ratio_post, DR_T_POA, &
                                   DR_T_PRA, DR_temp, DR_temp_post, DRAO_P, DRAO_T_POA, DRAO_T_PRA,   &
                                   DS_LRLV, DS_LRRT, DS_root, DS_shoot, fol_frac, LCA_DS, LCA_ratio,  &
                                   LRLV, LRRT, root_frac, stem_frac

  integer :: stmob, & ! remoblise stem C to labile (1 = on)
             turnover_labile_switch    ! begin turnover of labile C

double precision ::           ts_length, & ! time step length in hours
                      stock_seed_labile, & ! size of seed planted at sowning date
                                gpp_acm, & ! local GPP value as gCm-2.t-1
                                 deltat, & ! decimal day time step
                   stock_surface_litter, & ! used in radiative transfer
                    stock_storage_organ, & ! storage organ C pool, i.e. the desired crop (gC.m--2)
                     stock_dead_foliage, & ! dead but still standing foliage (gC.m--2)
                        stock_resp_auto, & ! autotrophic respiration pool (gC.m--2)
                           stock_labile, & ! labile C pool (gC.m--2)
                          stock_foliage, & ! foliage C pool (gC.m--2)
                             stock_stem, & ! stem C pool (gC.m--2)
!                            stock_roots, & ! roots C pool (gC.m--2)
                           stock_litter, & ! litter C pool (gC.m--2)
                    stock_soilOrgMatter, & ! SOM C pool (gC.m--2)
                              resp_auto, & ! autotrophic respiration (gC.m-2.t-1)
                          resp_h_litter, & ! litter heterotrophic respiration (gC.m-2.t-1)
                   resp_h_soilOrgMatter, & ! SOM heterotrophic respiration (gC.m-2)
                                    npp, & ! net primary productivity (gC.m-2.t-1)
                              nee_dalec, & ! net ecosystem exchange (gC.m-2.t-1)
                              lai_local, & !
                                     DS, & ! Developmental state and initial condition
                                   cLCA, & ! leaf mass area (gC.m-2)
            mean_alloc_to_storage_organ, & ! rolling average allocation of GPP to storage organ (gC.m-2)
        mean_alloc_to_storage_organ_old, & ! ...same but previous value...
                     decomposition_rate, & ! decomposition rate (frac / hr)
                     frac_GPP_resp_auto, & ! fraction of GPP allocated to autotrophic carbon pool
                  turnover_rate_foliage, & ! turnover rate of foliage (frac/hr)
                     turnover_rate_stem, & ! same for stem
                   turnover_rate_labile, & ! same for labile 
                turnover_rate_resp_auto, & ! same for autotrophic C pool
                 resp_cost_labile_trans, & ! labile lost to respiration per gC labile to GPP
             mineralisation_rate_litter, & ! mineralisation rate of litter
      mineralisation_rate_soilOrgMatter, & ! mineralisation rate of SOM
                                  PHUem, & ! emergance value for phenological heat units
                                    PHU, & ! phenological heat units
                                 DR_pre, & ! development rate coefficient DS 0->1
                                DR_post, & ! development rate coefficient DS 1->2
                                   tmin, & ! min temperature for development
                                   tmax, & ! max temperature for development
                                   topt, & ! optimum temperature for development
                                 tmin_v, & ! min temperature for vernalisation
                                 tmax_v, & ! max temperature for vernalisation
                                 topt_v, & ! optimim temperature for vernalisation
                                    VDh, & ! effective vernalisation days when plants are 50 % vernalised 
                                     VD, & ! count of vernalisation days
                               RDRSHMAX, & ! maximum rate of self shading turnover
                                   PHCR, & ! critical value of photoperiod for development
                                   PHSC, & ! photoperiod sensitivity
                                   raso, & ! rolling average for alloc to storage organ
                               max_raso, & ! maximum value for rolling average alloc to storage organ
                                  BM_EX, & ! 
                                     HI, & !
                                  yield, & ! crop yield (gC.m-2)
                     alloc_to_resp_auto, & ! amount of carbon to allocate to autotrophic respiration pool
                    turnover_rate_roots, & ! turnover over rate of roots interpolated each time step
                                gso_max, & !
                           max_raso_old, & !
                               raso_old, & !
            resp_cost_labile_to_foliage, & ! respiratory cost of moving carbon..from labile to foliage pools
            resp_cost_foliage_to_labile, & ! ..from foliage to labile pools
                              resp_rate, & ! rate of respiration at given temperature
                                 Cshoot, & !
                                     DR, & !
                        fol_frac_intpol, & !
                       stem_frac_intpol, & !
                               fP,fT,fV, & !                      
                                  remob, & !
                       root_frac_intpol, & !
                      shoot_frac_intpol, & !
                                 avtemp, & !
                 alloc_to_storage_organ, & !
                     litterfall_foliage, & !
                        litterfall_stem, & !
                       litterfall_roots, & !
                          decomposition, & !
                              npp_shoot, & !
                      alloc_from_labile, & !
                        alloc_to_labile, & !
                         alloc_to_roots, & !
                       alloc_to_foliage, & !
                          alloc_to_stem, & !
                                raremob, & !
                                  RDRSH, & !
                                  RDRDV, & !
                                    RDR


  ! 
  ! some hardcoded crop parameters
  ! 

  ! defines Q10 = 2 in exponential temperature response for heterotrophic
  ! respiration
  double precision, parameter :: resp_rate_temp_coeff = 0.0693d0
  ! residue fraction of leaves left post harvest
  double precision, parameter :: lv_res = 0.1d0
  ! residue fraction of stem left post harvest
  double precision, parameter :: st_res = 0.1d0
  ! LAI above which self shading turnover occurs
  double precision, parameter :: LAICR = 4.0d0
  ! allocation to storage organ relative to GPP 
  double precision, parameter :: rel_gso_max = 0.35d0
  ! minimum difference between DS and maturity for fol turnover
  double precision, parameter :: min_senescence_coef = 0.1d0

  ! 'save' indicates all these module variables should be held, even if the module itself goes out of scope.
  save 

contains
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! public procedures first.
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !

  subroutine CARBON_MODEL_CROP

    use gv_veg, only: lai
    use gv_clim, only: temp_top
    use gv_scale_declarations, only: time, avg_days_in_year, dble_one, dble_zero
    use gv_carbon_model, only: FLUXES,POOLS,GPP,NEE,pars &
                              ,leaf_growth,leaf_death

    ! The Data Assimilation Linked Ecosystem Carbon - Combined Deciduous
    ! Evergreen Analytical (DALEC_CDEA) model. The subroutine calls the
    ! Aggregated Canopy Model to simulate GPP and partitions between various
    ! ecosystem carbon pools. These pools are subject to turnovers /
    ! decompostion resulting in ecosystem phenology and fluxes of CO2

    implicit none

    ! declare local variables
!    double precision :: airt_weighting(3)
!    integer :: n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY

    ! POOLS are:
    ! 1 = labile
    ! 2 = foliar
    ! 3 = root
    ! 4 = wood
    ! 5 = litter
    ! 6 = som
    ! 7 = autotrophic
    ! 8 = storage organ C

    ! FLUXES are: 
    ! 1 = GPP
    ! 2 = temprate
    ! 3 = respiration_auto
    ! 4 = leaf production
    ! 5 = labile production
    ! 6 = root production
    ! 7 = wood production
    ! 8 = labile release
    ! 9 = alloc to storage
    ! 10 = leaf litter production
    ! 11 = woodlitter production
    ! 12 = rootlitter production
    ! 13 = respiration het litter
    ! 14 = respiration het som
    ! 15 = litter2som (decomposition)
    ! 16 = alloc to autotrophic pool

    ! PARAMETERS
    ! 16 values

    ! p(1) decomposition rate (frac/hr)
    ! p(2) Fraction of GPP allocated to autotrophic C pool
    ! p(3) DR coef for DS (0->1)
    ! p(4) DR coef for DS (1->2)
    ! p(5) turnover rate of foliage (frac/hr)
    ! p(6) Turnover rate of wood/stem (frac/hr)
    ! p(7) maximum rate of foliar turnover due to self shading
    ! p(8) effective vernalisation days when plant is 50 % vernalised
    ! p(9) mineralisation rate of som
    ! p(10) mineralisation rate of litter
    ! p(11) = log10(avgN) 
    ! p(12) = sow day
    ! p(13) = labile lost to respiration per gC labile top GPP
    ! p(14) = phenological heat units needed for emergence
    ! p(15) ! harvest day (doy)
    ! p(16) ! plough day (doy)
    ! p(17) ! leaf mass area (gC.m-2)
    ! p18,p19,p20,p21,p22,p23,p24,p25 = labile, foliar, roots, stem, litter,
    ! som,
    ! autotrophic and storage organ pools respectively
    ! p(26) ! min temperature for development
    ! p(27) ! max temperature for development
    ! p(28) ! optimum temperature for development
    ! p(29) ! min temperature for vernalisation
    ! p(30) ! max temperature for vernalisation
    ! p(31) ! optimim temperature for vernalisation
    ! p(32) ! critical value of photoperiod for development
    ! p(33) ! photoperiod sensitivity
    ! p(34) ! turnover rate of labile C
    ! p(35) ! turnover rate of autotrophic C

    ! only do initialisation once!
    if (time%steps_count == 1) then
       deltat = 1!/time%steps_per_day

       ! assigning initial conditions
       POOLS(1,1)=pars(18) ! Clabile
       POOLS(1,2)=pars(19) ! Cfoliar
       POOLS(1,3)=pars(20) ! Croots
       POOLS(1,4)=pars(21) ! Cstructural
       POOLS(1,5)=pars(22) ! Clitter
       POOLS(1,6)=pars(23) ! Csom
       POOLS(1,7)=pars(24) ! Cauto
       POOLS(1,8)=pars(25) ! Cstorage

       ! length of time step in hours..
       ts_length = time%seconds_per_step / 3600d0
   
       ! pair incoming variables to local module levels
       ! finally set some initial conditions
       ploughed        = .false.
       use_seed_labile = .false. ! whether to use seed labile for growth..
       sown            = .false. ! has farmer sown crop yet?
       emerged         = .false. ! has crop emerged yet?
       ! reset variables
       avtemp = 0d0
       yield = 0d0
       DS = -1d0
       DR = 0d0 
       fV = 0d0 ; fT = 0d0 ; fP = 0d0
       mean_alloc_to_storage_organ = 0d0
       mean_alloc_to_storage_organ_old = 0d0
       PHU = 0d0
       VD = 0d0
       BM_EX = 0d0
       HI = 0d0
       stock_dead_foliage = 0d0
       alloc_to_labile = 0d0
       stmob = 0
       max_raso = 0d0
       raso = 0d0
       RDRDV = 0d0 
       max_raso_old = 0d0
       raso_old  = 0d0

       ! parameters from file
       decomposition_rate                = pars(1) / 24d0  ! decomposition rate (day->hr)
       frac_GPP_resp_auto                = pars(2)  ! fraction of GPP allocated to autotrophic carbon pool
       DR_pre                            = pars(3)  ! development rate coefficient DS (0->1)
       DR_post                           = pars(4)  ! development rate coefficient DS (1->2)
       turnover_rate_foliage             = pars(6) / 24d0 ! pars(5)  ! turnover_rate of foliage (day->hr)
       turnover_rate_stem                = pars(6) / 24d0 ! turnover rate of stem (day->hr)
       RDRSHMAX                          = pars(7) / 24d0 ! maximum rate of foliar turnover due to self shading (day->hr)
       VDh                               = pars(8)  ! effective vernalisation days when plants are 50 % vernalised 
       mineralisation_rate_litter        = pars(9) / 24d0 ! mineralisation rate litter (day->hr)
       mineralisation_rate_soilOrgMatter = pars(10)/ 24d0 ! mineralisation rate som (day->hr)
       sow_day                           = nint(mod(pars(12),avg_days_in_year)) ! sow day (doy)
       resp_cost_labile_trans            = pars(13) ! labile lost to respiration per gC labile to GPP
       PHUem                             = pars(14) ! phenological heat units required for emergence
       harvest_day                       = nint(mod(pars(15),avg_days_in_year)) ! nint(mod(pars(15),365.25)) ! harvest day (doy)
       plough_day                        = nint(mod(pars(16),avg_days_in_year)) !nint(mod(pars(12)-2.0,365.25)) ! plough day (doy)
       cLCA                              = pars(17) ! leaf mass area (gC.m-2)
       stock_labile                      = pars(18) ! labile C
       stock_foliage                     = pars(19) ! foliar C
       stock_roots                       = pars(20) ! root C
       stock_stem                        = pars(21) ! stem / wood C
       stock_litter                      = pars(22) ! litter C
       stock_soilOrgMatter               = pars(23) ! som C
       stock_resp_auto                   = pars(24) ! autotrophic resp pool
       stock_storage_organ               = pars(25) ! storage organ (i.e. desired crop)
       tmin                              = pars(26)-273.15d0 ! min temperature for development
       tmax                              = pars(27)-273.15d0 ! max temperature for development
       topt                              = pars(28)-273.15d0 ! optimum temperature for development
       tmin_v                            = pars(29)-273.15d0 ! min temperature for vernalisation
       tmax_v                            = pars(30)-273.15d0 ! max temperature for vernalisation
       topt_v                            = pars(31)-273.15d0 ! optimim temperature for vernalisation
       PHCR                              = pars(32) ! critical value of photoperiod for development
       PHSC                              = pars(33) ! photoperiod sensitivity
       turnover_rate_labile              = pars(34)/ 24d0 ! turnover rate labile C (day->hr)
       turnover_rate_resp_auto           = pars(35)/ 24d0 ! turnover rate of autotrophic carbon for respiration (day->hr)

       ! assume that 50 % of litter is present on the surface...questionable...
       stock_surface_litter = stock_litter * 0.5d0

    endif ! if POOLS allocated

    ! 
    ! Begin looping through each time step
    ! 

    ! update local LAI value
    lai_local = sum(lai)

    ! load GPP value into local memory
    ! GPP (gC.m-2.day-1 -> gCm-2.t-1)
    gpp_acm = GPP / time%steps_per_day

    ! pass relevant variables into crop module memory
    if (time%step == 1) then
        avtemp = 0d0
        ! daily average of allocation to storage organ (needed to determine max.
        ! storage organ growth rate)
        mean_alloc_to_storage_organ_old = mean_alloc_to_storage_organ
        mean_alloc_to_storage_organ     = 0d0
    endif
    avtemp = avtemp + (temp_top/time%steps_per_day)

    ! Heterotrophic respiration rate (Q10):  doubles with 
    ! 10 degree temperature rise resprate from soil file = 0.0693
    resp_rate = 0.5d0 * exp( resp_rate_temp_coeff * temp_top )

    ! determine development stage (DS)
    if (time%step == time%steps_per_day) call development_stage(deltat)
    ! determine the carbon partitioning based on development stage
    call carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)
    ! begin carbon allocation for crops
    call calc_pools_crops(DS_LRRT,LRRT)
    ! conduct management updates at the end of the day
    if (time%step == time%steps_per_day) call management_dates

    ! calculate the NEE 
    NEE = nee_dalec * time%steps_per_day

    ! GPP (gC.m-2.day-1)
    FLUXES(time%step,1) = GPP
    ! temprate (i.e. temperature modified rate of metabolic activity))
    FLUXES(time%step,2) = resp_rate
    ! autotrophic respiration (gC.m-2.day-1)
    FLUXES(time%step,3) = resp_auto * time%steps_per_day
    ! leaf production rate (gC.m-2.day-1)
    FLUXES(time%step,4) = alloc_to_foliage* time%steps_per_day
    ! labile production (gC.m-2.day-1)
    FLUXES(time%step,5) = (alloc_to_labile + remob)* time%steps_per_day
    ! root production (gC.m-2.day-1)
    FLUXES(time%step,6) = alloc_to_roots * time%steps_per_day
    ! wood production 
    FLUXES(time%step,7) = alloc_to_stem * time%steps_per_day
    ! labile production
    FLUXES(time%step,8) = (alloc_from_labile + resp_cost_labile_to_foliage)* time%steps_per_day
    ! alloc to storage organ
    FLUXES(time%step,9) = alloc_to_storage_organ * time%steps_per_day
    ! total leaf litter production
    FLUXES(time%step,10) = litterfall_foliage * time%steps_per_day
    ! wood litter production
    FLUXES(time%step,11) = litterfall_stem * time%steps_per_day
    ! total root litter production
    FLUXES(time%step,12) = litterfall_roots * time%steps_per_day
    ! respiration heterotrophic litter
    FLUXES(time%step,13) = resp_h_litter * time%steps_per_day
    ! respiration heterotrophic som
    FLUXES(time%step,14) = resp_h_soilOrgMatter * time%steps_per_day
    ! litter to som
    FLUXES(time%step,15) = decomposition * time%steps_per_day
    ! alloc to autotrophic pool
    FLUXES(time%step,16) = (frac_GPP_resp_auto * gpp_acm) * time%steps_per_day
    ! extracted C from management
    FLUXES(time%step,17) = yield * time%steps_per_day
 
    ! values needed for SPA LAI update
    leaf_growth = alloc_to_foliage ; leaf_death = litterfall_foliage

    ! labile pool
    POOLS(1,1) = stock_labile
    ! foliar pool
    POOLS(1,2) = stock_foliage
    ! wood pool
    POOLS(1,4) = stock_stem
    ! root pool
    POOLS(1,3) = stock_roots
    ! litter pool
    POOLS(1,5) = stock_litter
    ! som pool
    POOLS(1,6) = stock_soilOrgMatter
    ! autotrophic pool 
    POOLS(1,7) = stock_resp_auto
    ! storage organ pool
    POOLS(1,8) = stock_storage_organ

  end subroutine CARBON_MODEL_CROP
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine calc_pools_crops(DS_LRRT,LRRT)

    use gv_scale_declarations, only: time, dble_zero, dble_one
    use gv_carbon_model, only: dleaf

    ! Allocated GPP to NPP and various carbon pools. Based !
    ! this on physiological responses to temperature       !
    ! vernalisation, and photoperiod.                      !

    implicit none

    ! arguments
    double precision, dimension(:), intent(inout) :: DS_LRRT, & !
                                                        LRRT    !

    ! local variables
    double precision :: decomp_efficency,  decomposition_surface_litter

    ! turnover rate of fine roots is now equal to the
    ! loss rate of roots (Penning de Vries, 1989)..
    turnover_rate_roots = interpolate( DS , DS_LRRT , LRRT , 5 ) / 24d0

    ! if sown turn on labile / seed turnover for growth
    if (sown) then
        ! turnover on
        turnover_labile_switch = 1
    else
        ! turnover off
        turnover_labile_switch = 0
    endif

    ! Initialise..
    resp_cost_foliage_to_labile = 0d0 ; yield = 0d0

    ! respiratory cost of C transfer from labile pool to short-term pool (NPP)
    ! (gC.m-2.t-1)
    resp_cost_labile_to_foliage = turnover_rate_labile * resp_cost_labile_trans * resp_rate &
                                * ts_length * dble(turnover_labile_switch)
    resp_cost_labile_to_foliage = stock_labile * min(dble_one,resp_cost_labile_to_foliage)

    ! allocation flux from labile C pool to NPP (gC.m-2.t-1)
    alloc_from_labile = turnover_rate_labile * ( 1d0 - resp_cost_labile_trans ) * resp_rate &
                      * ts_length * dble(turnover_labile_switch)
    alloc_from_labile = stock_labile * min(dble_one,alloc_from_labile)

    ! When GPP is higher than seed C content, remaining seed carbon enters
    ! litter
    ! C pool, as seedlings do not fully exhaust their seed (P. de Vries p 48)
    if ( ( gpp_acm .gt. alloc_from_labile ) .and. ( use_seed_labile ) ) then
        stock_litter = stock_litter + stock_labile
        stock_labile = 0d0
        use_seed_labile = .false.
    endif

    ! NPP as a fraction of GPP (1-.32=.68 or 68%) + allocation..
    npp = ( 1d0 - frac_GPP_resp_auto ) * gpp_acm + alloc_from_labile
    ! from labile pool; = SHORT-TERM POOL

    root_frac_intpol  = max(dble_zero,min(dble_one,root_frac_intpol))
    alloc_to_roots    = root_frac_intpol * npp         !
    shoot_frac_intpol = 1d0 - root_frac_intpol          ! 
    npp_shoot         = npp - alloc_to_roots           ! NPP remaining after root growth==SHOOT fraction
    alloc_to_foliage  = fol_frac_intpol  * npp_shoot   !
    alloc_to_stem     = stem_frac_intpol * npp_shoot   !
    alloc_to_storage_organ = max(dble_zero,npp_shoot - alloc_to_foliage - alloc_to_stem)
    if ( alloc_to_storage_organ > 0d0 ) then  ! allocation flux to storage organ limited by maximum growth rate
        gso_max  = ( stock_storage_organ + 0.5d0 ) * rel_gso_max / time%steps_per_day
        alloc_to_storage_organ = min( alloc_to_storage_organ , gso_max )
        if (sown) then
            alloc_to_labile = ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                            * ( 1d0 - resp_cost_labile_trans )
            resp_cost_foliage_to_labile = ( npp_shoot - alloc_to_foliage - alloc_to_stem - alloc_to_storage_organ ) &
                                        * resp_cost_labile_trans
        else
            alloc_to_labile             = 0d0
            resp_cost_foliage_to_labile = 0d0
        endif
    endif
    mean_alloc_to_storage_organ = mean_alloc_to_storage_organ + alloc_to_storage_organ

    ! set switches to (de)activate leaf, root and stem remobliization
    if (time%step == time%steps_per_day) then
        mean_alloc_to_storage_organ = mean_alloc_to_storage_organ / time%steps_per_day
        raso_old = raso
        ! running average of growth rate of storage organ..
        raso = ( mean_alloc_to_storage_organ + mean_alloc_to_storage_organ_old ) * 0.5d0
        max_raso_old = max_raso
        max_raso = max( raso , max_raso_old )
        ! Stem remobilisation triggered once running average of storage organ
        ! growth declines
        ! Second part prevents premature remobilisation
        if ( ( raso < raso_old ) .and. &
              ( mean_alloc_to_storage_organ .gt. ( mean_alloc_to_storage_organ_old + 0.5d0 ) / time%steps_per_day ) ) then
            stmob = 1
        else 
            stmob = 0
        endif
    endif

    ! Code for calculating relative death rate of leaves (RDR) as a 
    !  function of shading (RDRSH) or developmental stage (RDRT).

    ! GT 0 if LAI GT 4; 0. < RDRSH < RDRSHMAX (usually ~0.03)
    RDRSH = min( RDRSHMAX , max( dble_zero , RDRSHMAX * ( lai_local - LAICR ) / LAICR ) )
    if ( DS < 1d0 ) then
        RDRDV = 0d0
    else
       ! RDRDV dependant on DR and DS, values range typically between 0.02 <
       ! RDRDV < 0.25
!print*,"!! What RDRDV to use? !!"
!!$      RDRDV = DR /( max( 0.1 , 2. - DS ) )
!!$      RDRDV = RDRDV / 24. ! to get hourly senescence rate
       RDRDV = turnover_rate_foliage * ( 1d0 / ( ( max( 2d0 - DS , min_senescence_coef ) ) * 8d0 ) ) ** 2
    endif

    ! relative leaf death rate is the maximum value of the arguments RDRSH and
    ! RDRDV
    RDR = max( RDRSH , RDRDV )

    ! remobilization of foliar C and allocation to dead leaves pool (gC.m-2.t-1)
    litterfall_foliage = stock_foliage * min(dble_one,ts_length * RDR)
    litterfall_stem    = stock_stem    * min(dble_one,ts_length * DR * turnover_rate_stem * dble(stmob)) ! remobstem
    litterfall_roots   = stock_roots   * min(dble_one,ts_length * turnover_rate_roots)

    ! remobilized C to NPP (from both leaves and stems) (gC.m-2.t-1)
    remob   = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * ( 1d0 - resp_cost_labile_trans )
    ! respiratory cost of C transfer (conversion from starch to photosynthates)
    ! (gC.m-2.t-1)
    Raremob = ( litterfall_foliage * 0.5d0 + litterfall_stem ) * resp_cost_labile_trans

    ! for mass balance calculate the decompostion efficency
    decomp_efficency = decomposition_rate & 
                     / (decomposition_rate+mineralisation_rate_litter)

    ! total litter decomposition
    decomposition = stock_litter * (decomposition_rate+mineralisation_rate_litter) * resp_rate * ts_length
    ! heterotrophic respiration component 1: mineralisation of litter C pool
    ! (gC.m-2.t-1)
    resp_h_litter = decomposition * (1d0 - decomp_efficency)
    ! heterotrophic respiration component 2:  mineralisation of organic matter C
    ! pool (gC.m-2.t-1)
    resp_h_soilOrgMatter = stock_soilOrgMatter * min(dble_one,mineralisation_rate_soilOrgMatter * resp_rate * ts_length)
    ! decomposition of litter to soil organic matter (gC.m-2.t-1)
    decomposition = decomposition - resp_h_litter

    ! explicit calculation of the surface litter carbon pools for use in the
    ! radiative transfer scheme. This currently does not interfere with the
    ! normal stock_litter pool which still contains both surface and below
    ! surface litter
    ! total litter decomposition
    decomposition_surface_litter = stock_surface_litter * (decomposition_rate+mineralisation_rate_litter) * resp_rate * ts_length

    ! Recalculate Carbon Pools...

    stock_foliage       = max(dble_zero, stock_foliage + alloc_to_foliage - litterfall_foliage)
    dleaf               = alloc_to_foliage - litterfall_foliage ! relative leaf change (growth or senescence?)
    stock_stem          = max(dble_zero, stock_stem + alloc_to_stem - litterfall_stem)
    stock_storage_organ = max(dble_zero, stock_storage_organ + alloc_to_storage_organ)
    stock_roots         = max(dble_zero, stock_roots         + alloc_to_roots   - litterfall_roots)
    stock_litter        = max(dble_zero, stock_litter + litterfall_roots - resp_h_litter - decomposition)
    stock_soilOrgMatter = max(dble_zero, stock_soilOrgMatter + decomposition    - resp_h_soilOrgMatter)
    stock_dead_foliage  = max(dble_zero, stock_dead_foliage  + litterfall_foliage * 0.5d0)
    stock_labile        = max(dble_zero, stock_labile + alloc_to_labile  - alloc_from_labile - resp_cost_labile_to_foliage + remob)
    ! surface litter is a sub-component of the litter pool as such it is already
    ! accounted for the in mass balance via respiration so the lack of one here
    ! is not a mistake. surface litter only impacts the radiative transfer
    stock_surface_litter = stock_surface_litter - decomposition_surface_litter

    ! respiratory pool: new photosynthates are added (gC.m-2.t-1)
    stock_resp_auto = stock_resp_auto + frac_GPP_resp_auto * gpp_acm
    ! autotrophic respiration; Ra (7% of respiratory pool) (gC.m-2.t-1)
    resp_auto = stock_resp_auto * min(dble_one,turnover_rate_resp_auto * ts_length)
    ! respiratory pool reduced by Ra (amount of C respired by plant)
    stock_resp_auto = max(dble_zero, stock_resp_auto - resp_auto) 
    ! respiratory cost of C transfer from labile pool to short-term pool added
    ! to
    ! yield total autotrophic respiration (gC.m-2.t-1)
    resp_auto = resp_auto + resp_cost_labile_to_foliage + resp_cost_foliage_to_labile + Raremob

    ! nee (gC.m-2.t-1)
    nee_dalec = (resp_auto + resp_h_litter + resp_h_soilOrgMatter) - gpp_acm

  end subroutine calc_pools_crops
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine carbon_alloc_fractions(DS_shoot,DS_root,fol_frac,stem_frac,root_frac)

    use gv_scale_declarations, only: user_opts

    ! Determines carbon allocation fractions as a function !
    ! of developmental stage (DS).  Allocation fractions   !
    ! are from tables published in Penning de Vries (1989) !

    implicit none

    double precision, dimension(:), intent(inout) ::   DS_shoot, & !
                                                        DS_root, & !
                                                       fol_frac, & !
                                                      stem_frac, & !
                                                      root_frac    !

    ! local variables..
    double precision,dimension(:),allocatable :: frac_shoot, frac_root

    if ( sown ) then ! after sowing

       allocate( frac_shoot(size(DS_shoot)) , frac_root(size(DS_root)) )

       if (DS < 0d0 .and. user_opts%plant_func_type == 3) then

           ! Tubers only restriction
           fol_frac_intpol = 0d0
           stem_frac_intpol = 0d0
           shoot_frac_intpol = 0d0
           root_frac_intpol = 0d0
   
       else 

          ! loop over three crop "organs": 1) foliage 2) stems 3) root
          ! not necessary for storage organs, as all remaining C is allocated to
          ! these

          ! use different input for foliage and stem fractions, as they are
          ! relative to 
          ! the total shoot (or aboveground) allocation, root is relative to
          ! total plant 
          ! (above- and belowground) allocation..

          ! leaf development stages and corresponding fractions..

          frac_shoot = fol_frac
          ! interpolate between PdV allocation values with reference to 
          ! developmental stage (DS)..
          fol_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )
 
          ! stem DS and fracs..
          frac_shoot = stem_frac
          stem_frac_intpol = interpolate( DS , DS_shoot , frac_shoot , size(DS_shoot) )

          ! root DS and fracs..
          frac_root = root_frac
          root_frac_intpol = interpolate( DS , DS_root , frac_root , size(DS_root) )
  
       endif ! DS < 0.0 .and. user_opts%plant_func_type == 3

    endif ! sown or not

  end subroutine carbon_alloc_fractions
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine development_stage(days_in_step)

    ! Based on modified Wang & Engel model (Streck et al., 2003), !
    ! but with only 2 sub-phases, vegetative and reproductive     !
    ! (i.e. only two different DRmax).   O. Sus, May 2010.        !

    implicit none

    ! agruments
    double precision :: days_in_step

    ! local variables..
    double precision ::  doptmin, & ! Difference between optimum and minimum temperature
                         dmaxmin, & ! Difference between maximum and minimum temperature
                          dttmin, & ! Difference between daiy average and minimum temperatures
                       doptmin_v, & ! Difference between optimum and minimum vernalization temperatures
                       dmaxmin_v, & ! Difference between maximum and minimum vernalization temperatures
                        dttmin_v    ! Difference between daily average and minimum vernalization temperatures

    doptmin   = topt   - tmin   ! difference between optimal and minimum cardinal temperatures
    dmaxmin   = tmax   - tmin   ! difference between maximum and minimum cardinal temperatures
    dttmin    = avtemp - tmin   ! difference between daily average and minimum cardinal temperatures
    doptmin_v = topt_v - tmin_v ! same as above,
    dmaxmin_v = tmax_v - tmin_v !       but for vernalization 
    dttmin_v  = avtemp - tmin_v ! cardinal temperatures

    ! Calculation of developmental function values: vernalization (fV),
    ! temperature (fT) and
    ! photoperiod (fP) these values are multiplicative factors of DRmax (maximum
    ! developmental
    ! rate), each ranging between 0 (no development) and 1 (unrestricted
    ! development).

    ! Summation of vernalization days (VD), not before sowing and only if
    ! average temperature is within min and max cardinal temperatures..
    if ( ( avtemp .gt. tmin_v ) .and. ( avtemp .lt. tmax_v ) .and. sown ) then
        fV = vernalization( doptmin_v , dmaxmin_v , dttmin_v , days_in_step )
    endif

    ! Only calculate temperature coefficient if avtemp lies within (tmin,tmax)
    ! range.
    if ( (avtemp .gt. tmin ) .and. ( avtemp .lt. tmax ) ) then
        fT = temperature_impact( doptmin , dmaxmin , dttmin )
    else
        fT = 0d0
    endif

    fP = photoperiod_impact( PHCR , PHSC ) ! calculation of photoperiod coefficient

    if ( emerged .and. ( DS .lt. 2d0 ) ) then   ! sum up daily DR values between emergence and maturity (DS=2)

       if ( DS .lt. 1d0 ) then  ! in the vegetative phase (before flowering):

          DR = DR_pre * fT * fP   ! DR is affected by temperature, photoperiod...

          if ( vernal_calcs ) DR = DR * fV ! ...and vernalization (for winter cereals)

          DS = DS + (DR * days_in_step)    ! developmental stage (DS), calculated as the sum of daily developmental rates

       else    ! in the reproductive phase (after flowering):

          DR = DR_post * fT   ! DR is affected only by temperature

          DS = DS + (DR * days_in_step)

       endif

    endif ! emerged but not fully mature

  end subroutine development_stage
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine management_dates 

    use gv_scale_declarations, only: time, dble_zero, dble_one

    ! This routine should be called at the end of each day of a crops  !
    ! simulation.  It checks whether we should plough/sow/harvest, and !
    ! during the growing establishes when the crop will emerge after   !
    ! sowing, based on heat accumulation (Phenological Heat Units).    !

    implicit none

    ! local varibles
    double precision :: tmp
    logical :: plough_sanity,sow_sanity,harvest_sanity

    ! reset
    plough_sanity = .false. ; sow_sanity = .false. ; harvest_sanity = .false.

    ! spring crop
    if (sow_day < harvest_day .and. time%day < harvest_day) sow_sanity = .true.
    if (plough_day < harvest_day .and. time%day < harvest_day) plough_sanity = .true.
    if (harvest_day > sow_day) harvest_sanity = .true.
    ! winter crops
    if (sow_day > harvest_day) sow_sanity = .true.
    if (plough_day > harvest_day) plough_sanity = .true.
    if (harvest_day < plough_day .and. time%day < plough_day) harvest_sanity = .true.

    if ( .not. sown ) then

       ! fresh field...

       if ( plough_sanity .and. .not.ploughed .and. time%day >= plough_day ) then
          ! the field needs ploughing..
          call plough 

       elseif ( sow_sanity .and. time%day >= sow_day ) then

          ! ensure that the field has indeed been ploughed
          if (.not.ploughed) call plough

          ! the field needs sowing..
          sown = .true.
          ! this switch controls whether the labile carbon within the seed is used
          ! for growth
          use_seed_labile = .true.
          stock_labile = stock_seed_labile

       endif ! plough or sow?

    else ! sown or not

       ! crop in field..

       ! calculate when crop emerges..
       if ( .not. emerged ) then

          ! estimate emergence date based on the accumulated phenological heat
          ! units (PHU)
          ! where PHU is the (positive) heat over tmin..
          tmp = max( avtemp - tmin , dble_zero )!*days_per_step
          PHU = PHU + tmp

          ! set the development stage and emergence..
          if ( PHU >= PHUem ) then
              emerged = .true.
              DS = 0d0
          else
              emerged = .false.
              DS = -1d0
          endif

       endif ! emerged or not

       ! note that in this case harvest day has been fixed relative to the sow
       ! day
       if ( harvest_sanity .and. time%day >= harvest_day) then
          ! the field needs harvesting..
          call harvest
       endif ! 

    endif ! sown or not

  end subroutine management_dates
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! PROCEDURES BELOW ARE PRIVATE, IE THEIR USE IS LIMITED TO THIS MODULE
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine harvest

    use gv_veg, only: lai
    use gv_scale_declarations, only: time

    implicit none

    ! shoot biomass..
    Cshoot = stock_foliage + stock_stem + stock_storage_organ + stock_labile

    ! determine harvest index..
    HI = stock_storage_organ / Cshoot

    ! the stuff we actually want from the harvest...
    yield = stock_storage_organ

    ! the biomass that is harvested in addition to the storage-organ..
    BM_EX  = stock_foliage * ( 1d0 - lv_res )          &
              + stock_stem * ( 1d0 - st_res )          &
               + stock_dead_foliage * ( 1d0 - lv_res ) &
                + stock_labile

    ! what's left (will fall to the ground)..
    stock_litter  = stock_litter                     &
                    + stock_resp_auto                &
                     + stock_foliage * lv_res        &
                      + stock_stem * st_res          &
                       + stock_dead_foliage * lv_res

    ! keep track of surface litter component explicitly too
    stock_surface_litter = stock_foliage * lv_res        &
                          + stock_stem * st_res          &
                          + stock_dead_foliage * lv_res

    ! empty the plant stocks..
    stock_storage_organ = 0d0
    stock_foliage       = 0d0
    stock_stem          = 0d0
    stock_dead_foliage  = 0d0
    stock_labile        = 0d0
    stock_resp_auto     = 0d0
    lai                 = 0d0
    dead_lai            = 0d0

    ! roots stay in ground and slowly decompose (until/unless the field is
    ! ploughed)

    ! reset logical variables..
    sown    = .false.
    emerged = .false.
    ploughed = .false.
    DS = -1d0 ; fV = 0d0 ; fT = 0d0 ; fP = 0d0 ; VD = 0d0 

  end subroutine harvest
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function photoperiod_impact( PH_crit , PH_sens )

    use gv_Irradiance_Sunshade, only: daylength

    ! Function to determine the coefficient for !
    ! photoperiod impact on developmental rate. !
    ! From Streck et al., 2003                  !

    implicit none

    ! arguments..
    double precision,intent(in) :: PH_crit, & ! critical photoperiod below which no development occurs
                                   PH_sens    ! photoperiod sensitivity

    photoperiod_impact = max(0d0, 1d0 - exp ( - PH_Sens * ( daylength - PH_crit ) ))

  end function photoperiod_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine plough

    use gv_scale_declarations, only: time
    ! this s/r will reset various carbon pools, to mimic the effect of the
    ! farmer ploughing. !

    implicit none

    ! Move all plant stocks into the litter pool.
    ! ( many of these should already be empty after the harvest, )
    ! ( e.g. the stocks for labile, foliage, storage-organ stem. )
    stock_litter        = stock_litter + stock_dead_foliage &
                          + stock_foliage + stock_labile    &
                           + stock_roots + stock_stem       &
                            + stock_storage_organ

    stock_dead_foliage  = 0d0
    stock_foliage       = 0d0
    stock_labile        = 0d0
    stock_roots         = 0d0
    stock_stem          = 0d0
    stock_storage_organ = 0d0
    ! at the moment assume surface litter is ploughed into soil. This assumption
    ! is false in no-til management
    stock_surface_litter = 0d0

    ! Reset the development stage & phenological heat units..
    ploughed = .true. ; DS = -1d0 ; PHU = 0d0 
    max_raso = 0d0 ; raso = 0d0 ; max_raso_old = 0d0 ; raso_old = 0d0
    mean_alloc_to_storage_organ_old = 0d0 ; mean_alloc_to_storage_organ = 0d0

  end subroutine plough
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function temperature_impact( doptmin , dmaxmin , dttmin )

    ! Function to determine the coefficent for  !
    ! temperature impact on developmental rate. !
    ! From Streck et al., 2003.                 !

    implicit none

    ! arguments..
    double precision,intent(in) :: doptmin , dmaxmin , dttmin   ! temperature differences

    ! local variables..
    double precision :: a , nmr , dnr

    a   = log( 2d0 ) / ( log( ( dmaxmin ) / doptmin ) )

    nmr = 2d0 * ( ( dttmin ) ** a ) * ( doptmin ** a ) - ( ( dttmin ) ** ( 2d0 * a ) )

    dnr = doptmin ** ( 2d0 * a )

    temperature_impact = nmr / dnr

  end function temperature_impact
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function vernalization( doptmin_v , dmaxmin_v , dttmin_v , days_in_step )

    use gv_scale_declarations, only: dble_zero, dble_one

    ! Function to determine the coefficent for vernalization !
    ! impact on developmental rate. See Streck et al., 2003. !

    implicit none

    ! arguments..
    double precision,intent(in) :: dmaxmin_v , doptmin_v , dttmin_v & ! temperature differences
                                  ,days_in_step

    ! local variables..
    double precision :: a , dnr , fvn , nmr

    a   = log( 2d0 ) / ( log( ( dmaxmin_v ) / doptmin_v ) )
    nmr = 2d0 * ( ( dttmin_v ) ** a ) * ( doptmin_v ** a ) - ( ( dttmin_v ) ** (2d0 * a ) )
    dnr = doptmin_v ** ( 2d0 * a )
    fvn = nmr / dnr

    VD = VD + (fvn*days_in_step)

    ! final output value..
    vernalization = max( dble_zero , min( dble_one , ( VD ** 5 ) / ( ( VDh ** 5 ) + (VD ** 5 ) ) ) )

  end function vernalization
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function interpolate( x , reference_x , reference_y , row )

    ! Interpolation function.                    !
    ! x is input value, interpol is output value !
    ! reference_x/y are reference input data.    !

    implicit none

    ! arguments..
    integer,intent(in)             :: row
    double precision,intent(in)                :: x
    double precision,dimension(row),intent(in) :: reference_x , reference_y

    ! local variables..
    integer::i

    do i = 1 , row

       if ( x .le. reference_x(1) ) then
          interpolate = reference_y(1)
          exit
       endif

       ! cycling means growth rate remains constant between DS levels
       if ( ( x .gt. reference_x(i) ) .and. ( i .lt. row ) ) cycle

       if ( x .eq. reference_x(i) ) then
          interpolate = reference_y(i)
          exit
       endif

       if ( x .lt. reference_x(i) ) then
          interpolate = reference_y(i-1) + ( x - reference_x(i-1) ) &
                       * ( reference_y(i) - reference_y(i-1) )      &
                       / ( reference_x(i) - reference_x(i-1) )
          exit
       else
          interpolate = reference_y(row)
       endif

    enddo

  end function interpolate
!
!--------------------------------------------------------------------
!
end module carbon_model_crop_mod
