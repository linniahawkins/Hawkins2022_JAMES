! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module carbon_model_mod

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: CARBON_MODEL     &
           ,disturbance_residue_to_litter &
           ,disturbance_residue_to_som    &
           ,disturbance_loss_from_litter  &
           ,disturbance_loss_from_som     &
           ,itemp,ivpd,iphoto&
           ,extracted_C

  ! forest rotation specific info
  double precision, allocatable, dimension(:) :: extracted_C,itemp,ivpd,iphoto

  double precision, allocatable, dimension(:) :: disturbance_residue_to_litter, &
                                                 disturbance_residue_to_som,    &
                                                 disturbance_loss_from_litter,  &
                                                 disturbance_loss_from_som 
  integer :: gsi_lag, nodays
  double precision, allocatable, dimension(:) :: tmp_x, tmp_m
  double precision :: deltat,deltat_1 &
         ,Tfac,Photofac,VPDfac  ! oC, seconds, Pa

  ! local fire related variables
  double precision :: CFF(7) = 0d0, CFF_res(4) = 0d0    & ! combusted and non-combustion fluxes
         ,NCFF(7) = 0d0, NCFF_res(4) = 0d0  & ! with residue and non-residue seperates
         ,combust_eff(5)                & ! combustion efficiency
         ,rfac                            ! resilience factor

  ! local deforestation related variables
  double precision, dimension(4) :: post_harvest_burn   & ! how much burning to occur after
                       ,foliage_frac_res    &
                       ,roots_frac_res      &
                       ,rootcr_frac_res     &
                       ,stem_frac_res       &
                       ,branch_frac_res     &
                       ,Cbranch_part        &
                       ,Crootcr_part        &
                       ,soil_loss_frac

  double precision :: labile_loss,foliar_loss      &
         ,roots_loss,wood_loss         &
         ,labile_residue,foliar_residue&
         ,roots_residue,wood_residue   &
         ,C_total,labile_frac_res      &
         ,Cstem,Cbranch,Crootcr        &
         ,stem_residue,branch_residue  &
         ,coarse_root_residue          &
         ,soil_loss_with_roots

  integer :: reforest_day, harvest_management,restocking_lag

  ! local variables for GSI phenology model
  double precision :: tmp,gradient &
         ,deltaGPPfrac &
         ,fol_turn_crit,lab_turn_crit &
         ,gsi_history(22),just_grown

  contains
    !
    !----------------------------------------------------------------------
    !
    ! public procedures first.
    !
    !----------------------------------------------------------------------
    !
    subroutine CARBON_MODEL

      use gv_clim, only: temp_top
      use canopy_optimisation_mod, only: calculate_lai_marginal_return
      use gv_carbon_model, only: POOLS,FLUXES,GPP,NEE,pars &
                                ,gpp_lai_marginal,fire_fraction &
                                ,gpp_lai_marginal_reference &
                                ,clearance_fraction,clearance_management &
                                ,avg_min_airt,avg_dayl,avg_vpd_Pa &
                                ,leaf_growth,leaf_death
      use gv_scale_declarations, only: time

      ! The Data Assimilation Linked Ecosystem Carbon - Growing Season
      ! Index - Forest Rotation (DALEC_GSI_FR) model. 
      ! The subroutine calls the Aggregated Canopy Model to simulate GPP and 
      ! partitions between various ecosystem carbon pools. These pools are
      ! subject to turnovers / decompostion resulting in ecosystem phenology and fluxes of CO2

      implicit none

      integer :: i,f,n
    
      ! met drivers are:
      ! 1st run day
      ! 2nd min daily temp (oC)
      ! 3rd max daily temp (oC)
      ! 4th Radiation (MJ.m-2.day-1)
      ! 5th CO2 (ppm)
      ! 6th DOY
      ! 7th lagged precip
      ! 8th deforestation fraction
      ! 9th burnt area fraction
      ! 10th 21 day average min temperature (K)
      ! 11th 21 day average photoperiod (sec)
      ! 12th 21 day average VPD (Pa)
      ! 13th Forest management practice to accompany any clearing
      ! 14th Mean temperature

      ! POOLS are:
      ! 1 = labile (p18)
      ! 2 = foliar (p19)
      ! 3 = root   (p20)
      ! 4 = wood   (p21)
      ! 5 = litter (p22)
      ! 6 = som    (p23)

      ! p(30) = labile replanting
      ! p(31) = foliar replanting
      ! p(32) = fine root replanting
      ! p(33) = wood replanting

      ! FLUXES are: 
      ! 1 = GPP
      ! 2 = temprate
      ! 3 = Rauto
      ! 4 = NPP->leaf (disused)
      ! 5 = NPP->labile
      ! 6 = NPP->root
      ! 7 = NPP->wood
      ! 8 = labile->leaf
      ! 9 = leaffall factor
      ! 10 = leaflitter production
      ! 11 = wood->CWD
      ! 12 = root->litter
      ! 13 = Rhet_litter
      ! 14 = Rhet_som
      ! 15 = litter->som
      ! 16 = labile release factor
      ! 17 = carbon flux due to fire
      ! 18 = growing season index
      ! 19 = CWD->litter

      ! PARAMETERS
      ! 17+4(GSI) values

      ! p(1) Litter to SOM conversion rate
      ! p(2) Fraction of GPP respired
      ! p(3) GSI sensitivity for leaf growth
      ! p(4) Fraction of GPP-(Rauto+flab) allocated to roots 
      ! p(5) max leaf turnover rate (gC.m-2.day-1)
      ! p(6) Turnover rate of wood (fraction.day-1)
      ! p(7) Turnover rate of roots (fraction.day-1)
      ! p(8) Litter turnover rate 
      ! p(9) SOM turnover rate  
      ! p(10) Parameter in exponential term of temperature 
      ! p(11) Canopy efficiency parameter - C_eff (part of ACM) not used here
      ! p(12) = max labile turnover rate (gC.m-2.day-1)
      ! p(13) = Fraction allocated to Clab 
      ! p(14) = min temp threshold (GSI) 
      ! p(15) = max temp threshold (GSI)
      ! p(16) = min photoperiod threshold (GIS) 
      ! p(17) = LCA (gC.m-2) not used here
      ! p(24) = max photoperiod threshold (GSI)
      ! p(25) = min VPD threshold (GSI)
      ! p(26) = max VPD threshold (GSI)
      ! p(27) = minimum GPP benefit of increased LAI for labile allocation to be allowed
      ! p(28) = fraction of Cwood which is Cbranch
      ! p(29) = fraction of Cwood which is Ccoarseroot
      ! p(30) = 
      ! p(31) = 
      ! p(32) = 
      ! p(33) = 
      ! p(34) = GSI sensitivity for leaf scenescence
      ! p(35) = GSI have I just left a growing state (>1)
      ! p(36) = initial GSI value
      ! p(38) = turnover rate of CWD
  
      ! variables related to deforestation
      ! labile_loss = total loss from labile pool from deforestation
      ! foliar_loss = total loss form foliar pool from deforestation
      ! roots_loss = total loss from root pool from deforestation
      ! wood_loss = total loss from wood pool from deforestation
      ! labile_residue = harvested labile remaining in system as residue
      ! foliar_residue = harested foliar remaining in system as residue
      ! roots_residue = harvested roots remaining in system as residue
      ! wood_residue = harvested wood remaining in system as residue
      ! coarse_root_residue = expected coarse woody root left in system as residue

      ! parameters related to deforestation
      ! labile_frac_res = fraction of labile harvest left as residue
      ! foliage_frac_res = fraction of foliage harvest left as residue
      ! roots_frac_res = fraction of roots harvest left as residue
      ! wood_frac_res = fraction of wood harvest left as residue
      ! Crootcr_part = fraction of wood pool expected to be coarse root
      ! Crootcr_frac_res = fraction of coarse root left as residue
      ! soil_loss_frac = fraction determining Csom expected to be physically
      ! removed along with coarse roots

      ! in first model time step only
      if (.not.(allocated(gpp_lai_marginal))) then
          ! time variables to memory
          nodays = sum(time%days_in_year)
          deltat = 1d0/time%steps_per_day
          deltat_1 = time%steps_per_day

          ! allocate dimensions and zero
          allocate(gpp_lai_marginal(time%steps_per_day) &
                  ,gpp_lai_marginal_reference(time%steps_per_day))
          POOLS = 0d0 ; FLUXES = 0d0
          ! assign initial conditions to pools
          POOLS(1,1)=pars(18) 
          POOLS(1,2)=pars(19)
          POOLS(1,3)=pars(20)
          POOLS(1,4)=pars(21)
          POOLS(1,5)=pars(22)
          POOLS(1,6)=pars(23)
          POOLS(1,7)=pars(37)
          ! calculate initial GSI value
          FLUXES(:,18)=calculate_gsi(pars(14),pars(15),pars(16) &
                                    ,pars(24),pars(25),pars(26) &
                                    ,avg_min_airt(1),avg_dayl(1),avg_vpd_Pa(1))

          ! initial values for deforestation variables
          labile_loss = 0d0    ; foliar_loss = 0d0
          roots_loss = 0d0     ; wood_loss = 0d0
          labile_residue = 0d0 ; foliar_residue = 0d0
          roots_residue = 0d0  ; wood_residue = 0d0
          stem_residue = 0d0   ; branch_residue = 0d0
          reforest_day = 0
          soil_loss_with_roots = 0d0
          coarse_root_residue = 0d0
          post_harvest_burn = 0d0

          ! now load the hardcoded forest management parameters into their locations

          ! Parameter values for deforestation variables
          ! scenario 1
          ! harvest residue (fraction); 1 = all remains, 0 = all removed
          foliage_frac_res(1) = 1d0
          roots_frac_res(1)   = 1d0
          rootcr_frac_res(1) = 1d0
          branch_frac_res(1) = 1d0
          stem_frac_res(1)   = 0d0 ! 
          ! wood partitioning (fraction)
          Crootcr_part(1) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
          ! Black et al 2009; Morison et al 2012)
          Cbranch_part(1) =  0.20d0 ! (Ares & Brauers 2005)
          ! actually < 15 years branches = ~25 %
          !          > 15 years branches = ~15 %.
          ! Csom loss due to phyical removal with roots 
          ! Morison et al (2012) Forestry Commission Research Note
          soil_loss_frac(1) = 0.02d0 ! actually between 1-3 %
          ! was the forest burned after deforestation
          post_harvest_burn(1) = 1d0 

          !## scen 2
          ! harvest residue (fraction); 1 = all remains, 0 = all removed
          foliage_frac_res(2) = 1d0
          roots_frac_res(2)   = 1d0
          rootcr_frac_res(2) = 1d0
          branch_frac_res(2) = 1d0
          stem_frac_res(2)   = 0d0 ! 
          ! wood partitioning (fraction)
          Crootcr_part(2) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
          ! Black et al 2009; Morison et al 2012)
          Cbranch_part(2) =  0.20d0 ! (Ares & Brauers 2005)
          ! actually < 15 years branches = ~25 %
          !          > 15 years branches = ~15 %.
          ! Csom loss due to phyical removal with roots 
          ! Morison et al (2012) Forestry Commission Research Note
          soil_loss_frac(2) = 0.02d0 ! actually between 1-3 %
          ! was the forest burned after deforestation
          post_harvest_burn(2) = 0d0
    
          !## scen 3
          ! harvest residue (fraction); 1 = all remains, 0 = all removed
          foliage_frac_res(3) = 0.5d0
          roots_frac_res(3)   = 1d0
          rootcr_frac_res(3) = 1d0
          branch_frac_res(3) = 0d0
          stem_frac_res(3)   = 0d0 ! 
          ! wood partitioning (fraction)
          Crootcr_part(3) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
          ! Black et al 2009; Morison et al 2012)
          Cbranch_part(3) =  0.20d0 ! (Ares & Brauers 2005)
          ! actually < 15 years branches = ~25 %
          !          > 15 years branches = ~15 %.
          ! Csom loss due to phyical removal with roots 
          ! Morison et al (2012) Forestry Commission Research Note
          soil_loss_frac(3) = 0.02d0 ! actually between 1-3 %
          ! was the forest burned after deforestation
          post_harvest_burn(3) = 0d0

          !## scen 4
          ! harvest residue (fraction); 1 = all remains, 0 = all removed
          foliage_frac_res(4) = 0.5d0
          roots_frac_res(4)   = 1d0
          rootcr_frac_res(4) = 0d0
          branch_frac_res(4) = 0d0
          stem_frac_res(4)   = 0d0 ! 
          ! wood partitioning (fraction)
          Crootcr_part(4) = 0.32d0 ! Coarse roots (Adegbidi et al 2005;
          ! Black et al 2009; Morison et al 2012)
          Cbranch_part(4) =  0.20d0 ! (Ares & Brauers 2005)
          ! actually < 15 years branches = ~25 %
          !          > 15 years branches = ~15 %.
          ! Csom loss due to phyical removal with roots 
          ! Morison et al (2012) Forestry Commission Research Note
          soil_loss_frac(4) = 0.02 ! actually between 1-3 %
          ! was the forest burned after deforestation
          post_harvest_burn(4) = 0d0

          ! for the moment override all paritioning parameters with those coming from
          ! CARDAMOM
          Cbranch_part = pars(28)
          Crootcr_part = pars(29)

          ! declare fire constants (labile, foliar, roots, wood, litter)
          combust_eff(1) = 0.1d0 ; combust_eff(2) = 0.9d0
          combust_eff(3) = 0.1d0 ; combust_eff(4) = 0.5d0
          combust_eff(5) = 0.3d0 ; rfac = 0.5d0

          if (.not.allocated(disturbance_residue_to_som)) then
              allocate(disturbance_residue_to_litter(1), &
                       disturbance_residue_to_som(1),    &
                       disturbance_loss_from_litter(1),  &
                       disturbance_loss_from_som(1))
          endif
          disturbance_residue_to_litter = 0d0 ; disturbance_loss_from_litter = 0d0
          disturbance_residue_to_som = 0d0 ; disturbance_loss_from_som = 0d0

          ! assign climate sensitivities
          fol_turn_crit=pars(34)-1d0
          lab_turn_crit=pars(3)-1d0
          just_grown=pars(35)

          ! calculate some values once as these are invarient between DALEC runs
          if (.not.allocated(tmp_x)) then
              ! 21 days is the maximum potential so we will fill the maximum potential
              ! + 1 for safety
              allocate(tmp_x(22),tmp_m(nodays))
              do f = 1, 22
                 tmp_x(f) = f
              end do
              do n = 1, nodays
                 ! calculate the gradient / trend of GSI
                 if (n < 21) then
                     ! integrate over the available run period if less than 21 days
                     tmp_m(n) = n-1
                 else
                     ! or over the 21 day period. NOTE that this is a simplification
                     ! of the CARDAMOM code because time step here will never be > 1
                     ! day
                     tmp_m(n) = 20
                 endif ! for calculating gradient
              end do ! calc daily values once
              ! allocate GSI history dimension
              gsi_lag=nint(max(2d0,maxval(tmp_m)))
          end if ! .not.allocated(tmp_x)
          ! assign our starting value
          gsi_history = pars(36)-1d0

      endif ! first time here

      ! 
      ! Begin GPP allocation
      ! 
 
      ! GPP (gC.m-2.day-1); from SPA
      FLUXES(time%step,1) = GPP
      ! temprate (i.e. temperature modified rate of metabolic activity))
      FLUXES(time%step,2) = exp(pars(10)*temp_top)
      ! autotrophic respiration (gC.m-2.day-1)
      FLUXES(time%step,3) = pars(2)*FLUXES(time%step,1)
      ! leaf production rate (gC.m-2.day-1)
      FLUXES(time%step,4) = 0d0
      ! labile production (gC.m-2.day-1)
      FLUXES(time%step,5) = (FLUXES(time%step,1)-FLUXES(time%step,3)-FLUXES(time%step,4))*pars(13)
     ! root production (gC.m-2.day-1)
      FLUXES(time%step,6) = (FLUXES(time%step,1)-FLUXES(time%step,3)-FLUXES(time%step,4)-FLUXES(time%step,5))*pars(4)
      ! wood production 
      FLUXES(time%step,7) = FLUXES(time%step,1)-FLUXES(time%step,3)-FLUXES(time%step,4)-FLUXES(time%step,5)-FLUXES(time%step,6)

      ! accumulate GPP marginal return for lai increment
      tmp = pars(12)*FLUXES(time%steps_per_day,18) ; tmp = POOLS(1,1)*(1d0-(1d0-tmp))
      tmp = (POOLS(1,2)+tmp)/pars(17) ! new total LAI
      gpp_lai_marginal(time%step) = calculate_lai_marginal_return(tmp)
      gpp_lai_marginal_reference(time%step) = GPP !calculate_lai_marginal_return(POOLS(1,2)/pars(17))
      ! update phenology and make leaf allocation at the end of the day
      if (time%step == time%steps_per_day) then

          ! calculate fractional GPP return for LAI increment: NOTE conditional
          ! to protect against NaN
          if (sum(gpp_lai_marginal) == sum(gpp_lai_marginal_reference)) then
              deltaGPPfrac = 0d0
          else
              if (sum(gpp_lai_marginal_reference) > 0d0) then
                  deltaGPPfrac = ((sum(gpp_lai_marginal) - sum(gpp_lai_marginal_reference))/sum(gpp_lai_marginal_reference))
              else
                  deltaGPPfrac = sum(gpp_lai_marginal) - sum(gpp_lai_marginal_reference)
              endif
          endif
          

          ! GSI added to fortran version by TLS 24/11/2014
          !   25/09/14 - JFE
          ! Here we calculate the Growing Season Index based on 
          ! Jolly et al. A generalized, bioclimatic index to predict foliar 
          ! phenology in response to climate Global Change Biology, Volume 11, page 619-632,
          ! 2005 (doi: 10.1111/j.1365-2486.2005.00930.x) 
          ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
          ! Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res.,
          ! 116, G03020, doi:10.1029/2010JG001545.
        
          ! It is the product of 3 limiting factors for temperature, photoperiod and
          ! vapour pressure deficit that grow linearly from 0 to 1 between a calibrated 
          ! min and max value. Photoperiod, VPD and avgTmin are direct input
          FLUXES(:,18) = calculate_gsi(pars(14),pars(15),pars(16),pars(24),pars(25),pars(26) &
                                      ,avg_min_airt(1),avg_dayl(1),avg_vpd_Pa(1))

          ! we will load up some needed variables
!          m = int(tmp_m(time%run_day))
!          ! update gsi_history for the calculation
!          if (time%run_day <= m) then
!              ! in first step only we want to take the initial GSI value only
!              gsi_history(m) = FLUXES(time%step,18)
!          else
              ! now we have filled the vector time to loop through in
              ! replacement
              do i = 2,gsi_lag
                 gsi_history(i-1)=gsi_history(i)
              enddo
              gsi_history(gsi_lag) = FLUXES(time%step,18)
!          endif

          ! calculate gradient
          gradient = linear_model_gradient(tmp_x(1:(gsi_lag)),gsi_history(1:gsi_lag),gsi_lag)

          ! first assume that nothing is happening
          FLUXES(:,9) = 0d0  ! leaf turnover
          FLUXES(:,16) = 0d0 ! leaf growth

          ! now update foliage and labile conditions based on gradient calculations
          if (gradient < fol_turn_crit .or. FLUXES(time%step,18) == 0d0) then
              ! we are in a decending condition so foliar turnover
              FLUXES(:,9) = pars(5)*(1d0-FLUXES(:,18))
              just_grown = 0.5d0
          else if (gradient > lab_turn_crit) then
              ! we are in a assending condition so labile turnover
              FLUXES(:,16) = pars(12)*FLUXES(:,18)
              just_grown = 1.5d0
              ! determine if increase in LAI leads to an improvement in GPP greater
              ! than
              ! critical value, if not then no labile turnover allowed
              if ( deltaGPPfrac < pars(27) ) then
                  FLUXES(:,16) = 0d0
              endif
          else ! just grown?
              ! probaly we want nothing to happen, however if we are at the seasonal
              ! maximum we will consider further growth still
              if (just_grown >= 1d0) then
                 ! we are between so definitely not losing foliage and we have
                 ! previously been growing so maybe we still have a marginal return on
                 ! doing so again
                 FLUXES(:,16) = pars(12)*FLUXES(:,18)
                 ! determine if increase in LAI leads to an improvement in GPP
                 ! greater than critical value, if not then no labile turnover allowed
                 if ( deltaGPPfrac < pars(27) ) then
                     FLUXES(:,16) = 0d0
                 endif
              endif ! Just grown?
           endif ! gradient choice

          ! total labile release
          FLUXES(:,8)  = POOLS(1,1)*(1d0-(1d0-FLUXES(:,16))**deltat)*deltat_1
          ! total leaf litter production
          FLUXES(:,10) = POOLS(1,2)*(1d0-(1d0-FLUXES(:,9))**deltat)*deltat_1

      endif ! daily adjustment to phenology

      ! 
      ! turnovers which occur at each time step
      ! 

      ! total wood litter production
      FLUXES(time%step,11) = POOLS(1,4)*(1d0-(1d0-pars(6))**deltat)*deltat_1
      ! total root litter production
      FLUXES(time%step,12) = POOLS(1,3)*(1d0-(1d0-pars(7))**deltat)*deltat_1

      ! 
      ! those with temperature AND time dependancies
      ! 

      ! respiration heterotrophic litter
      FLUXES(time%step,13) = POOLS(1,5)*(1d0-(1d0-FLUXES(time%step,2)*pars(8))**deltat)*deltat_1
      ! respiration heterotrophic som
      FLUXES(time%step,14) = POOLS(1,6)*(1d0-(1d0-FLUXES(time%step,2)*pars(9))**deltat)*deltat_1
      ! litter to som
      FLUXES(time%step,15) = POOLS(1,5)*(1d0-(1d0-FLUXES(time%step,2)*pars(1))**deltat)*deltat_1
      ! CWD to litter
      FLUXES(time%step,19) = POOLS(1,7)*(1d0-(1d0-FLUXES(time%step,2)*pars(38))**deltat)*deltat_1
      ! calculate the NEE  (gC.m-2.day-1)
      NEE = (-FLUXES(time%step,1)+FLUXES(time%step,3)+FLUXES(time%step,13)+FLUXES(time%step,14))

      !
      ! update pools for next timestep
      ! 

      if (time%step == time%steps_per_day) then

          ! labile pool
          POOLS(1,1) = POOLS(1,1) + ( sum(FLUXES(:,5)) -sum(FLUXES(:,8)) )*deltat
          ! foliar pool
          POOLS(1,2) = POOLS(1,2) + ( sum(FLUXES(:,4)) -sum(FLUXES(:,10))+sum(FLUXES(:,8)) )*deltat
          ! values needed for outside lai allocation in SPA
          leaf_growth = sum(FLUXES(:,4)+FLUXES(:,8))*deltat ; leaf_death = sum(FLUXES(:,10))*deltat
          ! wood pool
          POOLS(1,4) = POOLS(1,4) + ( sum(FLUXES(:,7)) -sum(FLUXES(:,11)) )*deltat
          ! root pool
          POOLS(1,3) = POOLS(1,3) + ( sum(FLUXES(:,6)) -sum(FLUXES(:,12)) )*deltat
          ! litter pool
          POOLS(1,5) = POOLS(1,5) + ( sum(FLUXES(:,10))+sum(FLUXES(:,12))+sum(FLUXES(:,19)) &
                                     -sum(FLUXES(:,13))-sum(FLUXES(:,15)) )*deltat
          ! som pool
          POOLS(1,6) = POOLS(1,6) + ( sum(FLUXES(:,15))-sum(FLUXES(:,14)) )*deltat
          ! cwd pool
          POOLS(1,7) = POOLS(1,7) + ( sum(FLUXES(:,11))-sum(FLUXES(:,19)) )*deltat

      endif ! end of day pool

      ! 
      ! deal first with deforestation
      ! 

      if (time%run_day == reforest_day) then
          POOLS(1,1) = pars(30) 
          POOLS(1,2) = pars(31) 
          POOLS(1,3) = pars(32) 
          POOLS(1,4) = pars(33) 
      endif 

      if (clearance_fraction > 0d0) then
  
          ! pass harvest management to local integer
          harvest_management = int(clearance_management)
          ! assume that labile is proportionally distributed through the plant
          ! and therefore so is the residual fraction
          C_total = POOLS(1,2) + POOLS(1,3) + POOLS(1,4)
          ! partition wood into its components
          Cbranch = POOLS(1,4)*Cbranch_part(harvest_management)
          Crootcr = POOLS(1,4)*Crootcr_part(harvest_management)
          Cstem   = POOLS(1,4)-(Cbranch + Crootcr)
          ! now calculate the labile fraction of residue
          labile_frac_res = ( (POOLS(1,2)/C_total) * foliage_frac_res(harvest_management) ) & 
                          + ( (POOLS(1,3)/C_total) * roots_frac_res(harvest_management)   ) & 
                          + ( (Cbranch/C_total)      * branch_frac_res(harvest_management)  ) &
                          + ( (Cstem/C_total)        * stem_frac_res(harvest_management)    ) &
                          + ( (Crootcr/C_total)      * rootcr_frac_res(harvest_management)  ) 

          ! loss of carbon from each pools
          labile_loss = POOLS(1,1)*clearance_fraction
          foliar_loss = POOLS(1,2)*clearance_fraction
          roots_loss  = POOLS(1,3)*clearance_fraction
          wood_loss   = POOLS(1,4)*clearance_fraction
          ! transfer fraction of harvest waste to litter or som pools
          ! easy pools first
          labile_residue = POOLS(1,1)*clearance_fraction*labile_frac_res
          foliar_residue = POOLS(1,2)*clearance_fraction*foliage_frac_res(harvest_management)
          roots_residue  = POOLS(1,3)*clearance_fraction*roots_frac_res(harvest_management)
          ! explicit calculation of the residues from each fraction
          coarse_root_residue  = Crootcr*clearance_fraction*rootcr_frac_res(harvest_management)
          branch_residue = Cbranch*clearance_fraction*branch_frac_res(harvest_management)
          stem_residue = Cstem*clearance_fraction*stem_frac_res(harvest_management)
          ! now finally calculate the final wood residue
          wood_residue = stem_residue + branch_residue + coarse_root_residue 
          ! mechanical loss of Csom due to coarse root extraction                 
          soil_loss_with_roots = Crootcr*clearance_fraction*(1d0-rootcr_frac_res(harvest_management)) &
                               * soil_loss_frac(harvest_management)

          ! update living pools directly
          POOLS(1,1) = max(0d0,POOLS(1,1)-labile_loss)
          POOLS(1,2) = max(0d0,POOLS(1,2)-foliar_loss)
          POOLS(1,3) = max(0d0,POOLS(1,3)-roots_loss)
          POOLS(1,4) = max(0d0,POOLS(1,4)-wood_loss)
          ! then work out the adjustment due to burning if there is any
          if (post_harvest_burn(harvest_management) > 0d0) then
              ! first fluxes 
              ! LABILE 
              CFF(1) = POOLS(1,1)*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF(1) = POOLS(1,1)*post_harvest_burn(harvest_management)*(1d0-combust_eff(1))*(1d0-rfac)
              CFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*combust_eff(1)
              NCFF_res(1) = labile_residue*post_harvest_burn(harvest_management)*(1d0-combust_eff(1))*(1d0-rfac)
              ! foliar 
              CFF(2) = POOLS(1,2)*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF(2) = POOLS(1,2)*post_harvest_burn(harvest_management)*(1d0-combust_eff(2))*(1d0-rfac)
              CFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*combust_eff(2)
              NCFF_res(2) = foliar_residue*post_harvest_burn(harvest_management)*(1d0-combust_eff(2))*(1d0-rfac)
              ! root 
              CFF(3) = 0d0 !POOLS(1,3)*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF(3) = 0d0 !POOLS(1,3)*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
              CFF_res(3) = 0d0 !roots_residue*post_harvest_burn(harvest_management)*combust_eff(3)
              NCFF_res(3) = 0d0 !roots_residue*post_harvest_burn(harvest_management)*(1-combust_eff(3))*(1-rfac)
              ! wood 
              CFF(4) = POOLS(1,4)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(4) = POOLS(1,4)*post_harvest_burn(harvest_management)*(1d0-combust_eff(4))*(1d0-rfac)
              CFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF_res(4) = wood_residue*post_harvest_burn(harvest_management)*(1d0-combust_eff(4))*(1d0-rfac)
              ! litter 
              CFF(5) = POOLS(1,5)*post_harvest_burn(harvest_management)*combust_eff(5)
              NCFF(5) = POOLS(1,5)*post_harvest_burn(harvest_management)*(1d0-combust_eff(5))*(1d0-rfac)
              ! CWD  Using Combustion factors for wood
              CFF(7) = POOLS(1,7)*post_harvest_burn(harvest_management)*combust_eff(4)
              NCFF(7) = POOLS(1,7)*post_harvest_burn(harvest_management)*(1d0-combust_eff(4))*(1d0-rfac)
              ! fires as daily averages to comply with units 
              FLUXES(time%step,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5)+CFF(7) & 
                          +CFF_res(1)+CFF_res(2)+CFF_res(3)+CFF_res(4))/deltat

              ! update the residue terms
              labile_residue = labile_residue - CFF_res(1) - NCFF_res(1)
              foliar_residue = foliar_residue - CFF_res(2) - NCFF_res(2)
              roots_residue  = roots_residue  - CFF_res(3) - NCFF_res(3)
              wood_residue   = wood_residue   - CFF_res(4) - NCFF_res(4)
              ! now update NEE
              NEE=NEE+FLUXES(time%step,17)
          else
              FLUXES(time%step,17) = 0d0
              CFF = 0d0 ; NCFF = 0d0
              CFF_res = 0d0 ; NCFF_res = 0d0
          end if
          ! update all pools this time
          POOLS(1,1) = max(0d0, POOLS(1,1) - CFF(1) - NCFF(1) )
          POOLS(1,2) = max(0d0, POOLS(1,2) - CFF(2) - NCFF(2) )
          POOLS(1,3) = max(0d0, POOLS(1,3) - CFF(3) - NCFF(3) )
          POOLS(1,4) = max(0d0, POOLS(1,4) - CFF(4) - NCFF(4) )
          POOLS(1,5) = max(0d0, POOLS(1,5) + (labile_residue+foliar_residue+roots_residue) &
                                          + (NCFF(1)+NCFF(2)+NCFF(3))-CFF(5)-NCFF(5) )
          POOLS(1,6) = max(0d0, POOLS(1,6) + (soil_loss_with_roots) + (NCFF(4)+NCFF(5)+NCFF(7)))
          POOLS(1,7) = max(0d0, POOLS(1,7) + wood_residue - CFF(7) - NCFF(7) )
          ! some variable needed for the EDCs
          ! reallocation fluxes for the residues
          disturbance_residue_to_litter = (labile_residue+foliar_residue+roots_residue) & 
                                        + (NCFF(1)+NCFF(2)+NCFF(3))
          disturbance_loss_from_litter  = CFF(5)+NCFF(5)
          disturbance_residue_to_som    = NCFF(4)+NCFF(5)+NCFF(7)
          disturbance_loss_from_som     = soil_loss_with_roots
          ! harvested carbon from all pools
          FLUXES(time%step,20) = (wood_loss-(wood_residue+CFF_res(4)+NCFF_res(4))) &
                               + (labile_loss-(labile_residue+CFF_res(1)+NCFF_res(1))) &
                               + (foliar_loss-(foliar_residue+CFF_res(2)+NCFF_res(2))) &
                               + (roots_loss-(roots_residue+CFF_res(3)+NCFF_res(3)))
          ! total carbon loss from the system
          C_total = (labile_residue+foliar_residue+roots_residue+wood_residue+sum(NCFF)) &
                  - (labile_loss+foliar_loss+roots_loss+wood_loss+soil_loss_with_roots+sum(CFF))

          ! if total clearance occured then we need to ensure some minimum
          ! values and reforestation is assumed one year forward
          if (clearance_fraction > 0.99d0) then
              ! FC Forest Statistics 2015 lag between harvest and restocking ~ 2 year
              restocking_lag = 365*2
              reforest_day = min((time%run_day+restocking_lag), sum(time%days_in_year))
          endif ! if total clearance

      endif ! end deforestation info

      ! 
      ! then deal with fire
      ! 

      if (fire_fraction > 0d0) then
  
          ! first fluxes 
          ! LABILE 
          CFF(1) = POOLS(1,1)*fire_fraction*combust_eff(1)
          NCFF(1) = POOLS(1,1)*fire_fraction*(1d0-combust_eff(1))*(1d0-rfac)
          ! foliar 
          CFF(2) = POOLS(1,2)*fire_fraction*combust_eff(2)
          NCFF(2) = POOLS(1,2)*fire_fraction*(1d0-combust_eff(2))*(1d0-rfac)
          ! root 
          CFF(3) = 0d0 ! POOLS(1,3)*fire_fraction*combust_eff(3)
          NCFF(3) = 0d0 ! POOLS(1,3)*fire_fraction*(1-combust_eff(3))*(1-rfac)
          ! wood 
          CFF(4) = POOLS(1,4)*fire_fraction*combust_eff(4)
          NCFF(4) = POOLS(1,4)*fire_fraction*(1d0-combust_eff(4))*(1d0-rfac)
          ! litter 
          CFF(5) = POOLS(1,5)*fire_fraction*combust_eff(5)
          NCFF(5) = POOLS(1,5)*fire_fraction*(1d0-combust_eff(5))*(1d0-rfac)
          ! fires as daily averages to comply with units 
          FLUXES(time%step,17)=(CFF(1)+CFF(2)+CFF(3)+CFF(4)+CFF(5))/deltat
    
          ! all fluxes are at a daily timestep 
          NEE=NEE+FLUXES(time%step,17)
      
          !// update pools
          ! Adding all fire pool transfers here 
          POOLS(1,1)=POOLS(1,1)-CFF(1)-NCFF(1)
          POOLS(1,2)=POOLS(1,2)-CFF(2)-NCFF(2)
          POOLS(1,3)=POOLS(1,3)-CFF(3)-NCFF(3)
          POOLS(1,4)=POOLS(1,4)-CFF(4)-NCFF(4)
          POOLS(1,5)=POOLS(1,5)-CFF(5)-NCFF(5)+NCFF(1)+NCFF(2)+NCFF(3)
          POOLS(1,6)=POOLS(1,6)+NCFF(4)+NCFF(5)+NCFF(7)
          POOLS(1,7)=POOLS(1,7)-CFF(7)-NCFF(7)
          ! some variable needed for the EDCs
          ! reallocation fluxes for the residues
          disturbance_residue_to_litter = (NCFF(1)+NCFF(2)+NCFF(3))
          disturbance_residue_to_som = (NCFF(4)+NCFF(5)+NCFF(7))
          disturbance_loss_from_litter  = CFF(5)+NCFF(5)

        end if ! end burnst area issues
!
  end subroutine CARBON_MODEL
  !
  !----------------------------------------------------------------------
  !
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_gsi(min_airt,max_airt,min_dayl,max_dayl,min_vpd,max_vpd &
                                         ,min_airt_in,avg_dayl_in,avg_vpd_Pa_in)

    ! GSI added to fortran version by TLS 24/11/2014
    ! 25/09/14 - JFE Here we calculate the Growing Season Index based on 
    ! Jolly et al. A generalized, bioclimatic index to predict foliar 
    ! phenology in response to climate Global Change Biology, Volume 11,
    ! page 619-632, 2005 (doi: 10.1111/j.1365-2486.2005.00930.x) 
    ! Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S.
    ! Denning (2011), A global reanalysis of vegetation phenology, J.
    ! Geophys. Res., 116, G03020, doi:10.1029/2010JG001545.

    ! It is the product of 3 limiting factors for temperature, photoperiod
    ! and vapour pressure deficit that grow linearly from 0 to 1 between a
    ! calibrated min and max value. 
    ! Photoperiod, VPD and avgTmin are direct input averaged over previous 21 days

    implicit none

    double precision, intent(in) :: min_airt,max_airt & ! min / max parameters for air temperature 
                       ,min_dayl,max_dayl & !  " " for day length
                       ,min_vpd,max_vpd   & !  " " VPD
                       ,min_airt_in       & ! oC
                       ,avg_dayl_in       & ! seconds
                       ,avg_vpd_Pa_in       ! Pa

    ! temperature limitation, then restrict to 0-1; correction for k->oC
    Tfac = (min_airt_in-(min_airt-273.15d0)) / (max_airt-min_airt)
    Tfac = min(1d0,max(0d0,Tfac))
    ! photoperiod limitation
    Photofac = (avg_dayl_in-min_dayl) / (max_dayl-min_dayl)
    Photofac = min(1d0,max(0d0,Photofac))
    ! VPD limitation
    VPDfac = 1d0 - ( (avg_vpd_Pa_in-min_vpd) / (max_vpd-min_vpd) )
    VPDfac = min(1d0,max(0d0,VPDfac))

    ! calculate and store the GSI index (assumed here to go to the end of day
    ! position)
    calculate_gsi = Tfac*Photofac*VPDfac

    return

  end function calculate_gsi
  !
  !----------------------------------------------
  !
  double precision function linear_model_gradient(x,y,interval)

    ! Function to calculate the gradient of a linear model for a given depentent
    ! variable (y) based on predictive variable (x). The typical use of this
    ! function will in fact be to assume that x is time.

    implicit none

    ! declare input variables
    integer :: interval
    double precision, dimension(interval) :: x,y

    ! declare local variables
    double precision :: sum_x, sum_y, sumsq_x,sum_product_xy

    ! calculate the sum of x
    sum_x = sum(x)
    ! calculate the sum of y
    sum_y = sum(y)
    ! calculate the sum of squares of x
    sumsq_x = sum(x*x)
    ! calculate the sum of the product of xy
    sum_product_xy = sum(x*y)
    ! calculate the gradient
    linear_model_gradient = ( (interval*sum_product_xy) - (sum_x*sum_y) ) &
                          / ( (interval*sumsq_x) - (sum_x*sum_x) )

    ! for future reference here is how to calculate the intercept
!    intercept = ( (sum_y*sumsq_x) - (sum_x*sum_product_xy) ) &
!              / ( (interval*sumsq_x) - (sum_x*sum_x) )

    ! don't forget to return to the user
    return

  end function linear_model_gradient
  !
  !----------------------------------------------------------------------
  !
end module carbon_model_mod
!
!------------------------------------------------------------------------
!
