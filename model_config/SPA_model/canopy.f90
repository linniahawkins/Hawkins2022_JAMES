! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module canopy

  !! >this module requires a summary< !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: timestep_calcs

  integer :: max_pass
  double precision :: dew_potential
  double precision :: dead_frac & ! fraction of dead but still standing canopy (crops only)
              ,bulk_conductance & ! bulk surface conductance (m.s-1)
                     ,ustar     & ! friction velocity (m.s-1)
                     ,windsurf    ! wind speed under the canopy (m.s-1)

  ! parameters
  double precision, parameter :: min_lai = 1.5d0, & ! minimum LAI used in momentum decay equations 
                         min_throughfall = 0.2d0, & ! minimum fraction of
                                                  ! precipitation which is through fall 
                             min_storage = 0.2d0    ! minimum canopy water (surface) storage (mm)

  ! marginal return calculation play
!  logical :: update_N_profile = .false.
!  double precision :: marginal_N_LCA_extra = 0.0 &
!                     ,marginal_N_LCA_less  = 0.0 &
!                     ,marginal_N_LCA_reference = 0.0
!  double precision, dimension(:), allocatable :: Nla_possible, LCA_possible

contains
  !
  !----------------------------------------------------------------------
  !
  subroutine timestep_calcs( time )

    !  timestep_calcs (previously called hour) runs first by timestep, !
    !  and within each step then runs through each layer.              !
    !  (calculate daily carbon gain for one layer)                     !

    use gv_clim,               only: ppt, coa, wetev, cp_air, temp_top, leaf_heat_conductance, lambda_bulk, &
                                     wdtop, ppfd_to_par, sw_rad, wind_spd
    use gv_hourscale,          only: totet, Qh, Qc ,Qe, Qn, Qm, Qs
    use gv_meteo,              only: la, psil, psis, rad_pass, air_density_kg, abs_pot_conv, dew_evap, live_frac, &
                                     roughl, displacement
    use gv_scale_declarations, only: grid, time_holder, user_opts, umol_to_g_carbon, g_to_umol_carbon, &
                                     g_to_mol_carbon, dble_zero, dble_one
    use gv_soil_structure,     only: weighted_SWP, through_fall, max_storage
    use gv_veg,                only: co2amb, gppt, lafrac, lafrac_dead, lai, LWPstore, LWPprevious, &
                                     respt, stock_roots, totevap, transt, transt_kg_s, leaf_temp, &
                                     leaf_temp_intermediate, can_sensible_layer,dead_lai, canopy_rad, &
                                     dew, tsk, system_energy_balance, sensible_heat, wetle, lemod, mmmod, &
                                     neemod, timeflux, flux, ess, potential_evaporation, netrad_day, canopy_height
    use leaf,                  only: leaf_balance, leaf_temperature, gs2
    use light,                 only: solar, solar_update_for_iterative_canopy, &
                                     leaf_sun, checkpar, soilnir, soilpar, soilnir_coef, soilpar_coef, &
                                     parrefl_crop, nirrefl_crop, nirabs, parabs_shade, parabs_sun, longem
    use soil_functions,        only: soil_processes, soil_conductance, lhs
    use gv_daily_averages,     only: daily_canopy_temperature,daily_sun_canopy_temperature, daily_shade_canopy_temperature, &
                                     daily_soil_conductance
    use log_tools,             only: write_log, msg_warning, message
    ! DALEC related 
    use CARBON_MODEL_MOD,      only: CARBON_MODEL
    use CARBON_MODEL_CROP_MOD, only: CARBON_MODEL_CROP, stock_surface_litter
    use gv_carbon_model,       only: FLUXES, POOLS, GPP, NEE, pars, leaf_growth, leaf_death &
                                    ,resp_auto, resp_h_litter, resp_h_soilOrgMatter 
    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    integer :: i
    double precision    :: litter_frac, & ! fraction of soil surface covered by litter (crops only)
                            frac_shade, & ! part of leaf not hit by sunlight
                              frac_sun, & ! part of leaf hit by sunlight
                            psil_store, & ! temporary storage of psil value
                                la_sun, & ! area of leaf sunlit
                              la_shade, & ! area of leaf shaded
                           net_can_rad, & ! Net canopy radiation (kW.m-2)
                            shade_psil, & ! psil value for shaded leaf
                              sun_psil, & ! psil value for sunlit leaf
                              sum_cond, & ! sum of conductance weighted by area for all surfaces in system (m.s-1)
                          tot_est_evap    ! total estimated evaporation

     double precision   ::  dead_sensible(grid%canopy), & ! dead foliage sensible heat (W.m-2)
                             sun_sensible(grid%canopy), & ! sun area leaf sensible heat (W.m-2)
                           shade_sensible(grid%canopy), & ! shade area leaf sensible heat (W.m-2)
                           dead_leaf_temp(grid%canopy), & ! dead leaf area temperature (oC)
                            leaf_temp_sun(grid%canopy), & ! leaf temperature in sun (oC)
                          leaf_temp_shade(grid%canopy)    ! leaf temperature in shade (oC)

    ! some initial conditions for the first timestep
    if (time%steps_count == 1) then
        tot_est_evap = 0d0 ; totet = 0d0 
    endif
    if (time%step == 1) then
        gppt = 0d0  ;  respt = 0d0  ;  transt = 0d0 ; wetev = 0d0
        system_energy_balance = 0d0 ;  transt_kg_s = 0d0
        sensible_heat = 0d0
    end if

    ! Latent heat of vapourisation (J.kg-1)
    ! TLS: WARNING WARNING lambda is calculated multiple times based on
    ! different assumptions of air temperature above canopy, leaf layer
    ! specific temperature, soil temperature etc. This is a problem!
    if (temp_top > 0d0) then
        lambda_bulk = 2501000.0d0-2364.0d0*temp_top
    else 
        ! sublimation must be occuring
        lambda_bulk = lhs 
    endif

    ! leaf death must occur first as the current LCA has greater relevance
    ! to the old leaves and not the new leaves about to be growth (there
    ! are any).

    ! assume decay
    if (sum(lai) > 0d0 .and. leaf_death > 0d0) then  ! do this only during leaf senescence
        ! for senescence, average LCA is used for calculation of
        ! reduction in lai
        lai = lai - ( leaf_death * lafrac / pars(17) )
    end if ! leaf scenescence

    ! update LAI based on growth and death for each C cycle model
    if ( leaf_growth > 0d0 ) then    ! do this only during leaf growth

        ! LAI calculated through accumulation of leaf mass changes divided by 
        ! leaf carbon per area parameter (LCA or pars(17)).
        ! This calculation allows LCA to be dynamic, thus any changes in LCA
        ! would only affect added leaf C/area 
        ! NOTE: that pars(17) continues to be used as profile LCA is
        ! expected to maintain the canopy average
        lai = lai + (leaf_growth * lafrac / pars(17))
    endif ! leaf growth

    ! calculate leaf area index...
    select case( user_opts%plant_func_type )

      ! but allow for specific modifications based on ecosystem type
      case( 0 , 1 ) ! deciduous or evergreen forest..

        ! no special instructions as yet      

      case(2 , 3) ! cereal crops/C4 grasslands(2) or tubers (3)

        ! 50 % carbon loss from lai goes to remobilise and 50 % goes to
        ! dead foliage
        ! while only 50 % carbon goes to dead lai the area itself
        ! is maintained
        if (user_opts%dead_lai_energy_balance .and. leaf_death > 0d0) then
            dead_lai = dead_lai + ( leaf_death * lafrac_dead / pars(17) )
        end if

        ! equation relating litter C mass to area cover (% --> fractional) and
        ! reflectance are from Nagler et al., (2003)
        ! equation for wheat is used
        if (stock_surface_litter > 0d0) then
            litter_frac = (1d-2*((-0.0007d0*(stock_surface_litter*2d0)**2d0)+(0.5053d0*stock_surface_litter*2d0)+7.4017d0))
            litter_frac = min(1d0,max(0d0, litter_frac))
            soilnir_coef = (soilnir*(1d0-litter_frac))+(nirrefl_crop*litter_frac)
            soilpar_coef = (soilpar*(1d0-litter_frac))+(parrefl_crop*litter_frac)
        else
            litter_frac = 0d0
            soilnir_coef = 1d0
            soilpar_coef = 1d0
        end if

    end select ! for pdf

    ! zero values for next time step
    leaf_death = 0d0 ; leaf_growth = 0d0

!    ! calculate Nitrogen profile based on total N fractional allocation over
!    ! canopy layers
!    totN = avN * sum( lai ) ! total N
!    Nla  = nfrac * totN     ! N per layer
!    where (lai > 0.0) Nla = Nla / lai

    ! sometime decay can overflow in systems with management (e.g. grazed
    ! grassland or crops)
    if (sum(lai) < 0d0) lai = 0d0

    ! determine fraction of living lai as opposed to dead; as happens in crops
    if (sum(dead_lai) > 0d0) then
        live_frac = sum(lai)/(sum(lai)+sum(dead_lai))
    else
        live_frac = 1d0
    end if
    ! dead frac is obs the reverse
    dead_frac = 1d0-live_frac

    ! no point linking these to LAI unless we are iterating the canopy as the
    ! impact is only significant, and required, when wet evap is coupled to the
    ! canopy
    if (user_opts%iterative_canopy) then

       ! determine maximum canopy storage & through fall fraction (mm); based on
       ! LAI do not allow canopy dew or wet surface evaporation to occur on small
       ! canopies
       if ((sum(lai)+sum(dead_lai)) > 0d0) then
           through_fall = max(min_throughfall,exp(-(sum(dead_lai)+sum(lai))*0.5d0)*exp(wind_spd**2*0.10d0))
       else
           through_fall = 1d0
       end if
       ! maximum canopy storage (mm); minimum is applied to prevent errors in
       ! drainage calculation. Assume minimum capacity due to wood
       max_storage = max(min_storage,(0.1d0*(sum(lai)+sum(dead_lai))))

    end if ! iterative canopy

    ! calculate canopy level boundary conductance for heat and water vapour
    call boundary
    ! carry out radiative transfer for longwave, NIR and PAR
    call solar
    ! carry out soil surface energy balance and hydrology (and for some reason
    ! canopy water drainage)
    call soil_processes( time , tot_est_evap )

    ! save weighted soil water potential before entering the loop
    psis = weighted_SWP
    ! remember current leaf water potentials for use in marginal return
    ! calculations
    LWPprevious = LWPstore

    ! Cycle through canopy hydraulogical & carbon cycles; providing a linking
    ! into the energy cycle. 
    ! Due to isothermal net radiation being calculated in solar(), this loop is
    ! ran 3 times 
    ! 1) to allow for a long wave radiation balance to be updated making net
    ! radiation; 
    !    Net radiation used to calculate net to calculate expected dew or wet
    !    canopy evaporation.
    ! 2) the impact on canopy temperature of dew formation or canopy
    ! evaporation. Finish here
    ! Needs to be updated to cycle until residual less than 10 W.m-2 or max iter
    if (user_opts%iterative_canopy) then
        max_pass = 2 !3 !5 
    else
        max_pass = 1
    end if

    ! now iterate the canopy how ever many as requested
    do rad_pass = 1, max_pass

       if (rad_pass == 1) then
         dew = 0d0 ; dew_evap = 0d0
         ! only need to assign this value if the iterative canopy is not in use.
         ! If we are iterating the canopy this will be calculated below in
         ! conjunction with the canopy energy balance
         if ( .not. user_opts%iterative_canopy ) then
             dew_evap = lambda_bulk * wetev( time%step ) / time%seconds_per_step
         end if
       end if

       ! Reset GPP, respiration, transpiration and net radiation fluxes for each
       ! grid box
       gppt(time%step) = 0d0 ; respt(time%step) = 0d0
       transt(time%step) = 0d0 ; transt_kg_s(time%step) = 0d0
       sensible_heat(time%step) = 0d0
       canopy_rad = 0d0 ; net_can_rad = 0d0 ; can_sensible_layer = 0d0 ; la = 0d0 
       frac_sun = 0d0 ; frac_shade = 0d0 ; sun_sensible = 0d0 ; shade_sensible = 0d0 
       dead_sensible = 0d0 ; leaf_temp_intermediate = 0d0 ; leaf_temp_sun = 0d0 
       leaf_temp_shade = 0d0 ; dead_leaf_temp = 0d0 ; checkpar = 0d0

       ! calculate dew formation as a canopy whole with lai corrections
       ! function first calculates (kg.m-2.t-1)
       if (user_opts%iterative_canopy .and. (sum(lai)+sum(dead_lai)) > 0d0) then
           ! we need an initial estimate of whole canopy net radiation
           ! net radiation = shortwave + longwave radiation
           canopy_rad = sum(nirabs + longem + ((parabs_sun + parabs_shade) / ppfd_to_par))
           ! now estimate intercepted rainfall / dew formation across the bulk
           ! canopy
           call dew_or_wet_canopy_evaporation(time)
           ! now we must clear this for the canopy interative processes
           canopy_rad = 0d0
       end if ! lai condidtion

       do i = 1 , grid%canopy

          psil       = LWPstore(i) ! load layer's LWP from previous timestep
          co2amb     = coa         ! load current CO2 concentration
          psil_store = psil        ! store LWP 

          ! assign some local vars for convenience..
          frac_sun   = leaf_sun(i)
          la_sun     = frac_sun * lai(i)
          frac_shade = (dble_one-leaf_sun(i))
          la_shade   = frac_shade * lai(i)

          ! apply restriction to prevent precision error that occurs at very
          ! small lai values
          if (sum(lai) <= 1d-30) then
              la_sun = 0d0 ; la_shade = 0d0
          end if

          !-- first do sunlit leaf area --
          la = la_sun
          ! if there is active leaf area..
          if ( (la > 0d0) .and. (tot_est_evap > 0d0) ) then
               call leaf_processes(i,frac_sun,la,sunlit=.true.)
          end if

          ! now, if there is any leaf at all calculate its impact on the canopy
          ! energy balance
          ! temperature and sensible heat flux
          if (la > 0d0) then
              if ( tot_est_evap <=  0d0) gs2 = 0.00004d0
              call leaf_temperature_hfx(i,frac_sun,la,sun_sensible(i),leaf_temp_sun(i),sunlit=.true.,dead=.false.)
          else
              leaf_temp_sun(i) = 0d0
              sun_sensible(i) = 0d0
              canopy_rad = canopy_rad
          end if

          if (rad_pass == max_pass .and. lafrac(i) > 0d0) then
              sun_psil = psil   ! record LWP in sunlit leaves
          else
              ! cleared for dead foliage passes in crops
              sun_psil = 0d0
          end if

          !-- now do shaded leaf area --

          if ( lafrac(i) > 0d0 ) then
              psil = psil_store    ! reset LWP for shade leaf calculation
          else
              ! cleared for dead foliage passes in crops
              psil = 0d0
          end if

          la = la_shade
          if ( ( la > 0d0 ) .and. ( tot_est_evap > 0d0 ) ) then
               call leaf_processes(i,frac_shade,la,sunlit=.false.)
          end if

          ! now, if there is any leaf at all calculate its impact on the canopy
          ! energy balance
          ! temperature and sensible heat flux
          if (la > 0d0) then
              if ( tot_est_evap <=  0d0) gs2 = 0.00004d0
              call leaf_temperature_hfx(i,frac_shade,la,shade_sensible(i),leaf_temp_shade(i) & 
                                       ,sunlit=.false.,dead=.false.)
          else
              leaf_temp_shade(i) = 0d0
              shade_sensible(i) = 0d0
              canopy_rad = canopy_rad
          end if

          ! check for dead foliage
          if ( lafrac(i) > 0d0 ) then
              shade_psil  = psil  ! record LWP in shaded leaves
          else
              shade_psil  = 0d0
              psil        = 0d0
              LWPstore(i) = 0d0
          end if
          ! update states at end only
          if ( rad_pass == max_pass .and. lafrac(i) > 0d0) then
              psil        = frac_sun * sun_psil + frac_shade * shade_psil  ! calculate final LWP
              LWPstore(i) = psil       ! update store of LWP for next timestep
          end if

          ! now calculate, crudely, the leaf temperature contribution of dead
          ! foliage (only in crops)
          if (dead_lai(i) > 0d0 .and. user_opts%dead_lai_energy_balance) then
              gs2 = 0.00004d0
              ! determine dead leaf fractin of total leaf area
              call leaf_temperature_hfx(i,frac_sun+frac_shade,dead_lai(i),dead_sensible(i),dead_leaf_temp(i) & 
                                       ,sunlit=.false.,dead=.true.)
          else
              dead_leaf_temp(i) = 0d0
              dead_sensible(i) = 0d0
              canopy_rad = canopy_rad
          end if
 
          ! updates to leaf temperatures when in final iteration
          if (rad_pass == max_pass) then

              ! update and scale canopy temperature and sensible heat based on
              ! sun, shade and dead foliage
              leaf_temp(i) = 0d0
              ! leaf temperature (oC) is a mean value of the canopy leaf
              ! temperature between sun and shade; with convertion from
              ! potential temperatures to absolute
              leaf_temp(i)=((((leaf_temp_sun(i)*frac_sun)+(leaf_temp_shade(i)*frac_shade)) &
                        *(dble_one-dead_frac))+(dead_leaf_temp(i)*dead_frac))
              ! determine canopy averaged sensible heat exchange (W.m-2)
              can_sensible_layer(i)=sun_sensible(i)+shade_sensible(i)+dead_sensible(i)
              ! store sun/shade canopy temperature weighted by LAI for later
              ! output
              if (allocated(daily_shade_canopy_temperature)) then
                  daily_canopy_temperature(time%step) = daily_canopy_temperature(time%step) + (leaf_temp(i) * (lai(i)/sum(lai))) 
                  daily_shade_canopy_temperature(time%step) = daily_shade_canopy_temperature(time%step) + &
                                                              (leaf_temp_shade(i) * frac_shade * (lai(i)/sum(lai)))
                  daily_sun_canopy_temperature(time%step) = daily_sun_canopy_temperature(time%step) + &
                                                              (leaf_temp_sun(i) * frac_sun * (lai(i)/sum(lai)))
                  if (i == grid%canopy) then
                      daily_canopy_temperature(time%step) = daily_canopy_temperature(time%step) - temp_top
                      daily_shade_canopy_temperature(time%step) = daily_shade_canopy_temperature(time%step) - temp_top
                      daily_sun_canopy_temperature(time%step) = daily_sun_canopy_temperature(time%step) - temp_top
                  end if
                  ! also output soil conductance here to avoid further
                  ! conditional statement
                  daily_soil_conductance(time%step) = soil_conductance
              end if ! allocated shade_canopy_temperature
          else

              ! intermediate canopy layer temperature used in iteration;
              ! maintain as potential temperature for consistency with canopy
              ! energy balance
              leaf_temp_intermediate(i)=((((leaf_temp_sun(i)*frac_sun)+(leaf_temp_shade(i)*frac_shade)) &
                        *(dble_one-dead_frac))+(dead_leaf_temp(i)*dead_frac))
              ! determine canopy averaged sensible heat exchange (W.m-2)
              can_sensible_layer(i)=sun_sensible(i)+shade_sensible(i)+dead_sensible(i)
          end if ! for rad pass condition

       enddo ! end canopy loop

       ! update long wave component of radiation balance based on calculated
       ! canopy temperatures 
       if (rad_pass < max_pass) then
           call solar_update_for_iterative_canopy
       end if

    enddo ! rad pass loops

    ! calculate mean_canopy_temperature
    do i = 1, grid%canopy 
       if (abs(temp_top-leaf_temp(i)) > 10d0 .and. lafrac(i) > 0d0 .and. sum(lai) > 0.1d0) then
           write(message,*)'layer = ',i,' leaf temperature ',leaf_temp(i),'more than 10C from air ',temp_top
           call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
       endif
    enddo 

    ! check model energy balance; canopy, soil, total system
    if (abs(canopy_rad-transt(time%step)-dew_evap-sum(can_sensible_layer)) > 10d0 .and. sum(lai) > 0.1d0) then
        write(message,*)'canopy energy balance residual', (canopy_rad-transt(time%step)-dew_evap-sum(can_sensible_layer))
        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end if

    if (abs(Qn-(-Qe)-(-Qs)-(-Qh)+Qc-Qm) > 10d0) then
        write(message,*)'soil surface energy balance residual', Qn-(-Qe)-(-Qs)-(-Qh)+Qc-Qm
        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end if

    ! calculate total ecosystem sensible heat flux (W.m-2)
    sensible_heat(time%step) = ((-Qh)+sum(can_sensible_layer)) 
    ! calculate total system energy balance error (W.m-2)
    system_energy_balance(time%step) = (canopy_rad+Qn)-(-Qh+sum(can_sensible_layer))-(-Qe+transt(time%step)+dew_evap) & 
                                     -(-Qs) + Qc - Qm

    if (abs(system_energy_balance(time%step)) > 10d0 .and. sum(lai) > 0.1d0) then
        write(message,*)'system energy balance residual', &
             (canopy_rad+Qn)-(-Qh+sum(can_sensible_layer))-(-Qe+transt(time%step)+dew_evap) & 
            -(-Qs) + Qc - Qm
        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end if

    ! Determine combined surface skin temperature (K) 
    ! First calculate the canopy contribution of temperature to the surface
    ! (weighted by area) and to the conductance 
    sum_cond = dble_zero
    do i = 1, grid%canopy
       sum_cond = sum_cond+(2d0*leaf_heat_conductance(i)*lai(i))
    enddo
    ! invert total sensible and cumulative conductance to determine effective
    ! skin temperature
    ! first determine effective difference between surface and air temperature
    tsk = ((-Qh)+sum(can_sensible_layer))/(air_density_kg*cp_air*(sum_cond+soil_conductance))
    ! now add air temperature back to determine effective skin temperature (K)
    tsk = tsk + ((temp_top+273.15d0)*abs_pot_conv)
    ! convert potential skin temperature back to absolute (K)      
    tsk = tsk / abs_pot_conv

    ! its unlikely that the surface will be more than 10 degrees different
    ! from the surrounding air
    if (abs(temp_top-(tsk-273.15d0)) > 10d0 ) then
        write(message,*)'skin temperature',tsk-273.15d0, 'unusually different from air', temp_top
        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end if

    ! Calculate the total evapotranspiration during this time period..
    totet   = 0.001d0 * time%seconds_per_step * transt_kg_s(time%step)        ! transpiration (kg.m-2.s-1 -> Mg.m-2.t-1 or m t-1)
    totevap = totevap + 1000d0 * totet                                      ! (mm t-1)

    ! now call allocate carbon (dalec) if we want to simulate C pools
    if (user_opts%use_dalec) then

       ! update GPP to gC.m-2.day-1 for use in DALEC
       GPP = gppt(time%step) * time%seconds_per_day * umol_to_g_carbon 

       ! The focus here has been to make the minimum differences in the DALEC
       ! science code as possible. There are, however, large changes to the
       ! interface between the codes to make the integration efficient with SPA.
       ! Key differences between the CARDAMOM version of the src and
       ! the SPA version will be made in their respective src codes.
       if (user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3) then
          ! i.e. crops
          call CARBON_MODEL_CROP
       else

          ! default (currently) DALEC_GSI_DFOL_CWD_FR 
          call CARBON_MODEL

       endif ! GSI_DFOL_CWD_FR .or. CROP model

       ! load root C stock into SPA specific variable (DO NOT REMOVE!)
       stock_roots = POOLS(1,3)

       ! but remember to convert to gC.m-2.step-1 for SPA output
       GPP = GPP / time%steps_per_day

       ! NEE gCm-2.day-1 -> umol m-2 s-1
       neemod( time%step ) = (NEE / time%seconds_per_day) * g_to_mol_carbon * 1d6
       ! autotrophic respiration gC.m-2.day-1 -> umolC.m-2.s-1
       resp_auto = (FLUXES(time%step,3) / time%seconds_per_day) * g_to_mol_carbon * 1d6
       ! heterotrophic (litter) respiration gC.m-2.day-1 -> umolC.m-2.s-1
       resp_h_litter = (FLUXES(time%step,13) / time%seconds_per_day) * g_to_mol_carbon * 1d6
       ! heterotrophic (som) respiration gC.m-2.day-1 -> umolC.m-2.s-1
       resp_h_soilOrgMatter = (FLUXES(time%step,14) / time%seconds_per_day) * g_to_mol_carbon * 1d6

    else 

    ! update GPP to gC.m-2.step-1 for SPA output
    GPP = gppt(time%step) * time%seconds_per_step * umol_to_g_carbon

    ! NEE gCm-2.day-1 -> umol m-2 s-1
    neemod( time%step ) = -9999d0
    ! autotrophic respiration gC.m-2.day-1 -> umolC.m-2.s-1
    resp_auto = -9999d0
    ! heterotrophic (litter) respiration gC.m-2.day-1 -> umolC.m-2.s-1
    resp_h_litter = -9999d0
    ! heterotrophic (som) respiration gC.m-2.day-1 -> umolC.m-2.s-1
    resp_h_soilOrgMatter = -9999d0

    end if ! use dalec

    ! track potental evaporation
    potential_evaporation(time%step) = penman_monteith(time%seconds_per_step,(canopy_rad+Qn),temp_top,bulk_conductance,wdtop)

    ! wet leaves evaporation. convert mm t-1 to Wm-2
    wetle( time%step )  = lambda_bulk * wetev( time%step ) / time%seconds_per_step
    ! sum modelled canopy & modelled soil LE, and wet canopy evap..
    lemod( time%step ) = transt( time%step ) + ess( time%step ) + wetle( time%step )
    ! ecosystem evaporation (kg.m-2.step-1)
    mmmod( time%step ) = lemod( time%step ) * ( time%seconds_per_step / lambda_bulk )
    ! net radiation (W.m-2)
    netrad_day( time%step ) = Qn + canopy_rad
    ! sum water flux from each canopy layer..(mmol m-2 GA s-1)
    timeflux( time%step ) = timeflux( time%step ) + sum(flux(time%step, 1:grid%canopy))

  end subroutine timestep_calcs
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
  subroutine boundary

    ! determine boundary layer conductances at each canopy layer !

    use gv_clim,               only: atmos_press, gbw, temp_bot, temp_top, wind_spd, leaf_heat_conductance
    use gv_scale_declarations, only: grid, user_opts, dble_zero
    use gv_veg,                only: dimen, layer_height, tower_height, lai, dead_lai, &
                                     canopy_base, canopy_height, canopy_level, lafrac, leaf_temp
    use gv_meteo,              only: roughl, abs_virt_conv, air_density_kg, displacement
    use soil_functions,        only: exchange_coefficient, vonkarman, grav, soil_conductance
    use log_tools

    implicit none

    ! local variables..
    integer :: i
    double precision    :: alpha, mult, store, thick, u(grid%canopy), xx
    double precision    :: soil_gwv & ! soil conductance weighted by veg fraction (m.s-1)
              ,tempx    & ! Canopy layer specific air temperature, based on decay through canopy profile (oC) 
              ,Dwv        ! Diffusion coefficient of water in air (m2.s-1); air temperature and pressure dependant
    !  variables for the more advanced (but not guarenteed to be better)
    !  boundary conditions
    double precision lm        & ! mixing length for vertical turbulence momentumn within the canopy
                ,ustar_spa     & ! debug
                ,min_wind      & ! minimum wind speed allowed (m.s-1); based on ustar = 0.1
                ,min_surf_wind & ! minimum surface wind speed, set differently to canopy minimum, 
                                 ! as Qc may drive additional turbulence
                ,zh_stability  & ! value for momentum stability coefficient from Monin-Obukhov similarity theorm
                ,lc            & ! length scale for vertical momentum absorption within the canopy
                ,beta          & ! ratio of canopy top wind speed and fricition velocity above canopy
                ,beta_max, beta_min &
                ,tempvk        & ! canopy layer specific absolute virtual air temperature (K)
                ,virt_leaf_temp& ! canopy layer specific virtual leaf temperature (K)
                ,Dh            & ! diffusion coefficient of heat in air (m2.s-1); air temperature and pressure dependant
                ,Uh            & ! wind speed at top of canopy (m.s-1)
                ,min_conductance & ! minimum conductance parameter for canopy level exchange (m.s-1)
                ,scalar_mixing_length & ! scalar mixixng length at the top of the canopy
                ,thermal_conductivity & ! thermal conductivity of air (W.m-2 K.-1)
                ,volumetric_thermal_expansion & ! per degree celcius
                ,nusselt_forced & ! Nusselt value under forced convection 
                ,nusselt_free & ! Nusselt value under free convection
                ,dynamic_viscosity & ! dynamic viscosity (kg.m-2.s-1)
                ,kinematic_viscosity & ! kinematic viscosity (m2.s-1)
                ,gh_free & ! canopy conductance for heat under free convection (all m.s-1)
                ,gh_forced & ! canopy conductance for heat under forced convection
                ,gv_free & ! canopy conductance for water vapour under free convection
                ,gv_forced & ! canopy conductance for water vapour under forced convection
                ,Sh_free & ! Sherwood number under free convection 
                ,Sh_forced & ! Sherwood number under forced convection
                ,Gr & ! Grasholf number
                ,Re & ! Reynolds number
                ,Pr &  ! Prandtl number
                ,min_canopy_base ! minimum height of canopy base for wind speed calculations (m)

   if (user_opts%canopy_MDecay == 2) then

       ! use more advanced boundary conductance

       ! set parameters
       min_wind = 0.1d0 ; min_surf_wind = 0.1d0 ; min_canopy_base = 0.15d0
       thermal_conductivity = 0.0257d0 ! thermal conductivity of air (W.m-2.K-1)
       volumetric_thermal_expansion = 3.66d-3
       Pr = 0.72d0 ; beta_max = 1d0 ; beta_min = 0.2d0

       ! Reset the boundary layer conductance
       gbw = 0d0 ; leaf_heat_conductance = 0d0

       ! find level at which LAI is no longer allocated to the canopy
       do i = 2, grid%canopy
          if (lafrac(i) == 0d0 .and. lafrac(i-1) /= 0d0) then
             ! height of lowest level at which canopy exists is canopy base (m)
             canopy_base=max(min_canopy_base,layer_height(i-1))
             canopy_level=i-1
             exit ! exit loop once done
          endif
       enddo

       ! both length scale and mixing length are considered to be constant within
       ! the canopy (under dense canopy conditions)
       ! calculate length scale for momentum absorption within the canopy; Harman &
       ! Finnigan (2007)
       if ((sum(lai)+sum(dead_lai)) > 0d0) then
           lc=((4d0*canopy_height)/(sum(lai)+sum(dead_lai)))
       else
           lc=vonkarman*tower_height
       endif
       ! calculate roughness length and zero-plane displacement as defined by
       ! canopy height and lai uses assumptions from Raupach (1994)
       call z0_displacement(displacement)
       ! based on Harman & Finnigan (2008); only integrated to 1 m above canopy
       call log_law_decay(Uh,lc,displacement,min_wind,scalar_mixing_length)

       ! caclulate this steps beta coefficient; minimum value see Finnigan & Harman (2007)
       beta = min(beta_max,max(ustar/Uh,beta_min)) 

       ! calculate mixing length for vertical momentum within the canopy;
       ! Harman & Finnigan (2008)
       if ((sum(dead_lai)+sum(lai)) > 0d0) then
           lm = max(canopy_height*0.02d0,2d0*(beta**3)*lc)
       else
           lm = canopy_height*vonkarman
       endif

       ! determine wind speed decay through the canopy and exchange conductances
       do i = 1, grid%canopy

          ! Within canopy wind speed (m.s-1); Harman & Finnigan (2008)
          if (sum(dead_lai)+(sum(lai)) > 0d0) then
              if (i == 1) then
                 u(i) = max(Uh*exp((beta*((layer_height(i+1)-canopy_height)*0.5d0))/lm),min_wind)
              else
                 u(i) = max(Uh*exp((beta*(layer_height(i)-canopy_height))/lm),min_wind)
              endif
          else
              u(i) = Uh
          endif

          ! TLS: used if explicit or parameterised modelling of internal canopy
          ! temperature profiles is implimented
          ! addition held off for further development although subroutines called
          ! below are included
!          if (i == grid%canopy) then
             ! now than whole canopy momentum has been calculated, calculate the
             ! within canopy (potential) temperature profile
!             call scalar_profile(lm,temp_profile,beta,i)
             ! convert potential temperature to absolute temperature profile (oC)
!             abs_temp_profile=((temp_profile+273.15d0)/abs_pot_conv)-273.15d0
!             if ((sum(dead_lai)+sum(lai)) < 0.5d0) then
!                 abs_temp_profile=abs_temp
!                 temp_profile=temptop
!             endif ! lai condition
!          endif ! can layer condition
       enddo ! for can layer

       ! now calculate canopy conductance from modelled within canopy temperature
       ! and momentum
       do i = 1, grid%canopy
          ! Determine canopy level temperature, decays through canopy (oC)
          tempx = temp_top !abs_temp_profile(i)
          ! load and calculate absolute virtual air temperature (K)
          tempvk = (temp_top+273.15d0)*abs_virt_conv !(abs_temp_profile(i)+273.15d0)*abs_virt_conv
          ! load and calculate absolute virtual leaf temperature (K), leaf temp
          ! from previous step
          virt_leaf_temp = (leaf_temp(i)+273.15d0)*abs_virt_conv
          ! Determine diffusion coefficient (m2s-1), pressure and temperature
          ! dependant. Jones p51;
          ! note units are tempx (oC), atmos_press (Pa), 0.0000242 = conversion to
          ! make diffusion specific for water vapor (um2.s-1)
          ! N.B. conversion for heat is 0.0000215 (um2.s-1)
          Dwv = 0.0000242d0*(((tempx+273.15d0)/293.15d0)**1.75d0)*101300d0/atmos_press
          Dh = 0.0000215d0*(((tempx+273.15d0)/293.15d0)**1.75d0)*101300d0/atmos_press
          ! specify minimum canopy conductance (m.s-1)
          min_conductance = 5d-3

          ! calculate conditions (stable/unstable) specific canopy conductance
          ! (Nikolov et al., (1995))
          if (user_opts%plant_func_type == 0) then
              ! calculate free and forced conductance for cylinder leaves (i.e.
              ! needles)
              ! first calculate nusselt value under forced convection conditions
              dynamic_viscosity = (((tempx+273.15d0)**1.5d0)/((tempx+273.15d0)+120d0))*1.4963d-6
              kinematic_viscosity = dynamic_viscosity/air_density_kg
              Re = (dimen(1)*u(i))/kinematic_viscosity
              ! 2.1 coefficient represents average shelter factor proposed in Grant
              ! (1984)
              nusselt_forced = (0.69d0*(Pr**(0.33d0))*(sqrt(Re)))/2.1d0
              ! second calculate nusselt value under free convection conditions
              ! (low wind)
              Gr = (volumetric_thermal_expansion*grav*((dimen(2)**3)*abs(tempvk-virt_leaf_temp)))/(kinematic_viscosity**2)
              nusselt_free = 0.488d0*(Gr**(0.25d0))
              ! update specific Sherwood numbers
              Sh_forced = 0.962d0*nusselt_forced; Sh_free = 0.962d0*nusselt_free
              ! dimen(2) is the shoot diameter for needle trees as under free
              ! convection the individual needles act as a single unit rather than
              ! seperates
              gh_forced = ((Dh*Sh_forced)/dimen(1))*0.5d0
              gh_free = ((Dh*Sh_free)/dimen(2))*0.5d0
              gv_forced = ((Dwv*Sh_forced)/dimen(1))*0.5d0
              gv_free = ((Dwv*Sh_free)/dimen(2))*0.5d0
              ! leaf conductance is the maximum conductance between free and forced
              ! convection
              leaf_heat_conductance(i)=max(gh_forced,gh_free)
              gbw(i)=max(gv_forced,gv_free)
          else
              ! calculate forced and free conductance for non-cylinder leaves
              ! first calculate nusselt value under forced convection conditions
              dynamic_viscosity = (((tempx+273.15d0)**1.5d0)/((tempx+273.15d0)+120d0))*1.4963d-6
              kinematic_viscosity = dynamic_viscosity/air_density_kg
              Re = (dimen(1)*u(i))/kinematic_viscosity
              ! 2.1 coefficient represents average shelter factor proposed in Grant
              ! (1984)
              nusselt_forced = (1.18d0*(Pr**(0.33d0))*(sqrt(Re)))
              ! second calculate nusselt value under free convection conditions
              ! (low wind)
              Gr = (volumetric_thermal_expansion*grav*((dimen(2)**3)*abs(tempvk-virt_leaf_temp)))/(kinematic_viscosity**2)
              nusselt_free = 0.921d0*(Gr**0.25d0)
              ! update specific Sherwood numbers
              Sh_forced = 0.962d0*nusselt_forced; Sh_free = 0.962d0*nusselt_free
              ! dimen(2) is the shoot diameter for needle trees as under free
              ! convection the individual needles act as a single unit rather than
              ! seperates
              gh_forced = ((Dh*Sh_forced)/dimen(1))*0.5d0
              gh_free = ((Dh*Sh_free)/dimen(2))*0.5d0
              gv_forced = ((Dwv*Sh_forced)/dimen(1))*0.5d0
              gv_free = ((Dwv*Sh_free)/dimen(2))*0.5d0
              ! leaf conductance is the maximum conductance between free and forced
              ! convection
              leaf_heat_conductance(i) = max(gh_forced,gh_free)
              gbw(i) = max(gv_forced,gv_free)
          endif

          ! restrict canopy conductance to minimum value; no fixed minimum
          ! currently applied
          leaf_heat_conductance(i) = leaf_heat_conductance(i)
          gbw(i) = gbw(i)
       enddo ! canopy loop

!       ! if calculate new surface wind speed for surface exchanges if there is
!       ! vegetation
!       if ((sum(lai)+sum(dead_lai)) > 0d0) then
!           ! surface wind speed is selected to be 50% through canopy; selected for
!           ! stability issues when soil surface exchanges to lower heights
!           windsurf=max(Uh*exp((beta*(canopy_base-canopy_height))/lm),min_surf_wind)
!       else ! lai condition
!           ! SPA surface energy balanace equations cannot handle low wind speed
!           ! values, causes over heating
!           windsurf=max(Uh,min_surf_wind)
!       endif

       ! NOTE either use soil_surface_resistance or exchange_coefficient.
       ! exchange_coefficient is the original SPA function.
       ! Nui & Yang (2004); modelling soil surface resistance through the under
       ! canopy space and canopy.
       call soil_surface_resistance(soil_conductance,lm,displacement,beta)
       ! convert resistance to soil conductance (m.s-1)
       soil_conductance = soil_conductance**(-1d0)

    ! call for soil conductance which has been weighted with open air soil
!       call exchange_coefficient(soil_gwv) ;  soil_conductance=soil_gwv

       ! rather than bulk conductance (as this is incorrect dew to Uh calculations)
       ! calculate sum of area weighted parallel conductance to exchange across the whole
       ! system; should be equivalent to bulk (ish)
       ! calculate bulk (equiv) conductance (m.s-1); used in dew calculation
       ! WARNING neutral conductance only
       bulk_conductance = dble_zero ! clear first
       bulk_conductance = sum(lai * gbw)
!!! series cumulation of resistance
!       do i = 1,grid%canopy
!          if (lai(i) > 0d0 ) then
!              bulk_conductance = bulk_conductance+(lai(i)*gbw(i))**(-1d0)
!          end if
!       end do
!       ! add soil surface component
!       bulk_conductance = (bulk_conductance+(soil_conductance)**(-1d0))**(-1d0)

       elseif (user_opts%canopy_MDecay == 1) then

          ! use original

          ! initialise..
          alpha  = 4d0
          roughl = 0.075d0 * tower_height 
          gbw    = 0d0
          leaf_heat_conductance = 0d0
          xx     = 0d0
          store  = 1d0
          bulk_conductance = 0d0

          ! calculate..
          do i = 1 , grid%canopy
            xx    = xx + 1d0
            tempx = temp_top - 0.1111d0 * ( xx - 1d0 ) * ( temp_top - temp_bot )
            Dwv   = 0.0000242d0 * ( ( ( tempx + 273.15d0 ) / 293.15d0 )**1.75d0 ) * 101300d0 / atmos_press    !Jones p 51
            ! windspeed in layer - no decay in windspeed in this open stand
            !             u(i,t)=windsp(t)
            mult   = exp( alpha * ( layer_height( i ) / tower_height - 1d0 ) )
            u( i ) = wind_spd * mult
            thick  = 0.004d0 * ( dimen(1) / u( i ) )**0.5d0    ! boundary layer thickness
            ! conductance to water vapour (m s-1 - not mol m-2 s-1 as in Amthor p.333)
            ! i.e. remove P/RT
            gbw( i ) = Dwv / thick
            ! approximate boundary layer conductance to heat, this should be changed
            ! later
            leaf_heat_conductance(i) = gbw(i) * 0.93d0
          end do

          ! rather than bulk conductance (as this is incorrect dew to Uh calculations)
          ! calculate sum of area weighted series conductance to exchange across the whole
          ! system; should be equivalent to bulk (ish)
          ! calculate bulk (equiv) conductance (m.s-1); used in dew calculation
          ! WARNING neutral conductance only
          bulk_conductance = dble_zero ! clear first
          do i = 1,grid%canopy
             if (lai(i) > 0d0 ) then
                 bulk_conductance = bulk_conductance+(lai(i)*gbw(i))**(-1d0)
             end if
          end do
          ! get the soil exchange coefficient, at least the current estimate of
          call exchange_coefficient(soil_gwv) ; soil_conductance = soil_gwv
          ! add soil surface resistance and then convert back into conductance
          bulk_conductance = (bulk_conductance+(soil_gwv)**(-1d0))**(-1d0)

    else 

          write(message,*) "Somehow the canopy_MDecay (1 or 2 value) has become set to an invalid value"
          call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    end if ! detailed boundary or not

  end subroutine boundary
  !
  !----------------------------------------------------------------------
  !
  subroutine log_law_decay(Uh,lc,displacement,min_wind,scalar_mixing_length)

    ! subroutine uses equations from Harman & Finnigan 2008 
    ! for log law decay of above canopy wind speed with MOST
    ! stability corrections (currently set for neutral only). 
    ! tower wind speed used to calculate friction velocities 
    ! and canopy top wind speed (m.s-1)

    use gv_scale_declarations, only: dble_zero,dble_one,dble_two
    use soil_functions,   only: vonkarman, grav
    use gv_clim,          only: wind_spd, cp_air
    use gv_veg,           only: tower_height, canopy_height, can_sensible_layer
    use gv_hourscale,     only: Qh
    use gv_meteo,         only: air_density_kg, roughl

    implicit none
    integer z
    double precision Uh           & ! canopy top wind speed, decayed from reference height 
        ,min_wind     & ! minimum wind speed allowance (m.s-1)
        ,coef_1       & ! gamma coefficient 1 for similarity function calculations (Businger et al (1971))
        ,coef_2       & ! beta coeffient for similarity function calculations (Businger et al (1971))
        ,most         & ! Monin Obukhov Similarity Theorm coefficient
        ,zeta         & ! zeta is the height / obukhov length variable used in stability determination. Height dependant
        ,zref_obukhov_length & ! Obukhov length at reference (tower) height
        ,obukhov_length & ! Obukhov length used in stability correction
        ,displacement & ! zero plane displacement height (m)
        ,step_size    & ! step size of integration
        ,virtual_origin & ! virtual origin for displacement height, e.g with height set to top of canopy
        ,ref_height_above & ! reference height above canopy (m)
        ,lc           & ! lenght scale of canopy for turbulence absorbance
        ,scalar_stability & ! Monin Obukov similarity theory stability for scalars at top of canopy
        ,scalar_mixing_length & ! mxing length for scalars at the top of canopy; assumed to be constant throughout the canopy
        ,iterate      & ! number of steps to iterate above canopy
        ,U_gradiant   & ! wind speed decay gradient at height z
        ,height_above   ! height above canopy top

    ! integrate resolution (m)
    step_size = -0.01d0
    ! initialise canopy top wind speed with speed from reference height
    Uh = wind_spd
    ! load coef values for Monin - Obukov Similarity theorem
    coef_1 = 16d0 ; coef_2 = 5d0

    ! calculate friction velocity at tower height (reference height ) (m.s-1) 
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ! Monin - Obukov length calculated as equ. 2.8 Qin et al (2002)
    zref_obukhov_length=(-vonkarman*tower_height*grav*(-Qh+sum(can_sensible_layer)))/(air_density_kg*cp_air*ustar**3)

    ! determine difference between reference height and canopy top to decay
    ! through
    ! commented out as roughness sublayer calculations not being used at this
    ! time
    ref_height_above=tower_height

    ! determine step size for decay multiplication; warning stops 1 meter above
    ! canopy
    iterate=(abs((ref_height_above-(canopy_height))/step_size))
    ! integrate from reference height to canopy top if the is likely above the
    ! roughness sublayer
    ! temperary solution until roughness sublayer parameterisation is added to
    ! decay
!    if (abs(tower_height - canopy_height) > 1.5) then
        ! iterate through the space above the canopy calculating the current
        ! gradient 
        ! and subtracting that from the current wind speed
        do z = 1, int(iterate)
           ! determine the current height above canopy
           height_above = ref_height_above+(step_size*(dble(z)-dble_one))
           ! determine stability functions for momentum (and heat)
           ! first calculate the Obukhov length, as stability varies with height
           if (zref_obukhov_length == 0d0) then
               ! when no turbulent fluxes present then totally stable conditions
               ! exist
               zeta = 0d0
           else
               obukhov_length=tower_height/zref_obukhov_length
               ! calculate new zeta variable for new height in decay
               zeta=height_above/obukhov_length
           endif
           ! coefficients are only value for -2. < zeta < 0. = unstable
           ! and 0 <= zeta < 1. = stable
           if (zeta < 0d0) then
               ! unstable
               most=(dble_one-(coef_1*max(zeta,-dble_two)))**(-0.25d0)
           else
               ! stable
               most=dble_one+(coef_2*min(zeta,dble_two))
           endif
           ! as roughness sublayer correction is not being inplemented at this
           ! time a standard arrangement of this equation is being used from
           ! Jones p67
           U_gradiant = (ustar/(vonkarman*(height_above-displacement)))*most
           ! determine new wind speed for height (m.s-1); most term determines
           ! sign of term
           Uh = Uh+(U_gradiant*step_size)
        enddo
!    endif

    ! set minimum value for wind speed at canopy top (m.s-1)
    if (Uh < min_wind) then
        Uh = min_wind
    endif

    ! determine stability functions for momentum (and heat)
    ! first calculate the Obukhov length, as stability varies with height
!    obukhov_length=tower_height/zref_obukhov_length
    ! calculate new zeta variable for new height in decay
!    zeta=height_above/obukhov_length
    ! coeffients are only value for -2. < zeta < 0. = unstable
    ! and 0. <= zeta < 1. = stable
    ! 1. at beginning of each expression is actually the Pr number; which here
    ! as in WRF is assumed to be = 1.
!    if (zeta < 0d0) then
!        scalar_stability=1d0*((1d0-(coef_1*max(zeta,-2d0)))**-0.5d0)
!    else
!        scalar_stability=1d0+(coef_2*min(zeta,2d0))
!    endif

    ! calculate virtual origin
!    virtual_origin=((ustar/Uh)**2)*lc
!    ! determine the canopy top mixing lengh for scalar (heat); to allow for
!    calculation of the canopy layer conductance via Finnigan (2004)
 !   scalar_mixing_length=(vonkarman*virtual_origin)/scalar_stability

  end subroutine log_law_decay
  !
  !----------------------------------------------------------------------
  !
  subroutine z0_displacement(displacement)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    use gv_scale_declarations, only: dble_one
    use gv_veg,          only: dead_lai, lai, canopy_height
    use soil_functions,  only: vonkarman
    use gv_meteo,        only: roughl

    implicit none
    double precision :: & 
          displacement  & ! zero plane displacement
         ,cd1           & ! canopy drag parameter; fitted to data
         ,Cr            & ! Roughness element drag coefficient
         ,Cs            & ! Substrate drag coefficient
         ,ustar_Uh_max  & ! Maximum observed ratio of (friction velocity / canopy top wind speed) (m.s-1)
         ,ustar_Uh      & ! ratio of fricition velocity / canopy top wind speed (m.s-1)
         ,Cw            & ! characterises roughness sublayer depth (m)
         ,phi_h           ! roughness sublayer inflence function

    ! set parameters
    cd1 = 7.5d0 ; Cs = 0.003d0 ; Cr = 0.3d0 ; ustar_Uh_max = 0.3d0 ; Cw = 2d0

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not varied
    displacement = (dble_one-((dble_one-exp(-(cd1*max(sum(lai)+sum(dead_lai),min_lai))**0.5d0)) &
                      /(cd1*max(min_lai,sum(lai)+sum(dead_lai)))**0.5d0))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh = min((Cs+Cr*max(min_lai,sum(lai)+sum(dead_lai))*0.5d0)**0.5d0,ustar_Uh_max)
    ! calculate roughness sublayer influence function; 
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    phi_h = log(Cw)-dble_one+Cw**(-dble_one)

    ! finally calculate roughness length (m), dependant on displacement, friction
    ! velocity and lai.
    roughl = ((dble_one-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

    ! sanity check
    if (roughl /= roughl) then
        write(*,*)"ERROR roughness length calculations"
        write(*,*)"Roughness length (m)", roughl, "Displacement height (m)", displacement
        write(*,*)"Canopy height (m)", canopy_height 
        write(*,*)"lai (m2.m-2)", sum(lai),"Dead lai (m2.m-2)",sum(dead_lai)
    endif

  end subroutine z0_displacement
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_surface_resistance(soil_resistance,lm,displacement,beta)

    use soil_functions, only: vonkarman, grav, soil_roughl
    use gv_clim,        only: cp_air
    use gv_meteo,       only: air_density_kg, roughl
    use gv_veg,         only: canopy_height, lai, dead_lai
    use gv_hourscale,   only: Qh
    use gv_scale_declarations, only: dble_zero, dble_one, dble_two

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction
    ! momentum & heat are integrated through the under canopy space and canopy
    ! air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002

    implicit none
    integer iter,iter_max
    double precision   canopy_decay & ! canopy decay coefficient for soil exchange
          ,beta         & ! ratio of u*/Uh
          ,lm           & ! mixing length (m)
          ,displacement & ! zero plane displacement height (m)
          ,Kh_canht     & ! eddy diffusivity at canopy height (m2.s-1)
          ,foliage_drag & ! foliage drag coefficient
          ,coef_1       & ! gamma coefficient 1 for similarity function calculations (Businger et al (1971))
          ,coef_2       & ! beta coeffient for similarity function calculations (Businger et al (1971))
          ,zeta         & ! zeta is the height / obukhov length variable used in stability determination. Height dependant
          ,iter_step    & ! integration step size (m)
          ,height       & ! height above the soil surface (m)
          ,most_soil    & ! Monin-Obukov similarity theory stability correction coef
          ,Kh           & ! eddy diffusivity at specific height in iteration (m2.s-1)
          ,soil_resistance ! soil resistance to atmopshere for exchange (s.m-1)

    ! load coef values for Monin-Obukov similarity theory for momentum
    coef_1 = 16d0 ; coef_2 = 5d0
    ! foliage drag coefficient
    foliage_drag = 0.2
    ! integrator step size (m)
    iter_step = 0.2 !
    ! calculate distance to be integrated through
    iter_max = (int(((displacement+roughl)-soil_roughl)/iter_step))-1 ! integer
    ! assign starting point height for intergration (m)
    height = displacement+roughl 

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1) 
    Kh_canht=vonkarman*ustar*(canopy_height-displacement) !Kaimal & Finnigan 1994; for near canopy approximation
    ! reset soil resistance
    soil_resistance = dble_zero
    do iter = 1, iter_max ! 

       ! Monin - Obukov length calculated as equ. 2.8 Qin et al (2002)
       ! Note that soil sensible heat here is positive if surface is absorbing
       ! energy
       zeta = (-dble_one*vonkarman*height*grav*Qh)/(air_density_kg*cp_air*ustar**3)
       ! coefficients are only value for -2. < zeta < 0. = unstable
       ! and 0 <= zeta < 1. = stable
       if (zeta < 0d0) then
           ! unstable
           most_soil = (dble_one-(coef_1*max(zeta,-dble_two)))**(-0.25d0)
       else
           ! stable
           most_soil = dble_one+(coef_2*min(zeta,dble_one))
       endif
       ! calculate canopy decay coefficient with stability correction
       ! Note that soil sensible heat here is positive if surface is absorbing
       ! energyhis is not consistent with canopy momentum decay done by Harman &
       ! Finnigan (2008)
       canopy_decay = (((foliage_drag*canopy_height*max(min_lai,(sum(lai)+sum(dead_lai))))/lm)**0.5d0)*(most_soil**0.5d0)
       ! approximation of integral for soil resistance; maximum soil surface
       ! resistance applied at 200 s.m-1, conductance = 5e-3 m.s-1
!       soil_resistance=canopy_height/(canopy_decay*Kh_canht) & 
!                      *(min(400d0,exp(canopy_decay*(1.0-(soil_roughl/canopy_height))))- &
!                        min(200d0,exp(canopy_decay*(1.0-((roughl+displacement)/canopy_height)))))
       ! calculate eddy diffusivity decay within the canopy
       Kh = Kh_canht*exp(-canopy_decay*(dble_one-(height/canopy_height))) !
       ! integrate increase in resistance over space (s.m-1)
       soil_resistance = soil_resistance + (iter_step/Kh) !
       ! now move down the height profile
       height = height - iter_step

    end do ! iterate down canopy

  end subroutine soil_surface_resistance
  !
  !----------------------------------------------------------------------
  !
  double precision function canopy_sensible(leaf_area,fraction_temperature,gbh,temp)

    use gv_meteo, only: abs_pot_conv, air_density_kg
    use gv_clim,  only: cp_air

    implicit none

    double precision ::  leaf_area, & ! Leaf area in pass (m2) (i.e. sun or shade)
     fraction_temperature, & ! Leaf temperature (oC) for given area
      potential_leaf_temp, & ! Leaf temperature converted to potential value
                      gbh, & ! Conductance to heat (m.s-1) based on temperature at top of the canopy / lowest atmospheric layer
                     temp, & ! Canopy layer specific air temperature (oC)
     potential_temp_layer, & ! canopy layer specific potential air temperature (oC)
                      rho, & ! Determine the temperature dependant density of air (g.m-3)
          convective_heat    ! Convective heat conductance (m.s-1)

    ! density of air (kg.m-3); temperature dependant based on canopy layer
    ! temperature
    rho = air_density_kg
    ! Determine the convective heat conductance used in sensible heat
    convective_heat = 2d0*gbh
    ! convert leaf temperature to potential value used in sensible calculations
    potential_leaf_temp = ((fraction_temperature+273.15d0)*abs_pot_conv)-273.15d0
    ! convert air temperature to potential value used in sensible calculations
    potential_temp_layer = ((temp+273.15d0)*abs_pot_conv)-273.15d0
    ! Canopy sensible heat (W.m-2) determined by constants, heat conductance and
    ! temp diff between leaf and canopy. This is also weighted by the canopy lai
    ! for the layer.
    ! uses cp_air converted to (J.kg-1.K-1) to match density of air in (kg.m-3)
    ! and convective conductance (m.s-1); p232 Jones, eq. 9.3
    canopy_sensible = cp_air * rho * ( (potential_leaf_temp-potential_temp_layer)*convective_heat ) * leaf_area

  end function canopy_sensible
  !
  !----------------------------------------------------------------------
  !
  subroutine dew_or_wet_canopy_evaporation(time)

    ! Subroutine deals with the calculation of wet canopy evaporation and dew
    ! formation. These fluxes are then appropriately restricted by other canopy
    ! processes and appropriate pools are updated    

    use gv_clim,               only: wetev, wdtop, temp_top, lambda_bulk
    use gv_hourscale,          only: canopy_store
    use gv_meteo,              only: rad_pass, dew_evap
    use gv_scale_declarations, only: time_holder, dble_zero, dble_one
    use gv_veg,                only: transt, transt_kg_s, canopy_rad, dew
    use gv_soil_structure,     only: max_storage

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    double precision :: canopy_store_old, & ! previous steps canopy water store (mm)
                    dew_mass, & ! mass of water exchange from penman_monteith (kg.m-2.-t)
                 layer_store, & ! canopy layer specific contribition to canopy water store
              layer_wet_evap    ! current canopy layer wet evaporation (kg.m-2.t-1)

    ! reset outputs
    dew = 0d0 ; dew_mass = 0d0 ; layer_wet_evap = 0d0
    dew_evap = 0d0 ; wetev(time%step) = 0d0 ; canopy_store_old = canopy_store

    ! calculate potential dew or evaporation (kg.m-2.t-1)
    dew_mass = penman_monteith(time%seconds_per_step,canopy_rad,temp_top,bulk_conductance,wdtop)
    dew_potential = (dew_mass/time%seconds_per_step)*lambda_bulk

    ! evaporation or dew fall
    if (dew_mass > 0d0) then

        if (canopy_store > 0d0) then
            ! convert to W.m-2 for initial use in correction
            ! application
            dew_evap = (dew_mass/time%seconds_per_step)*lambda_bulk
            ! determine potential evap ratio; rate restricted by
            ! ratio of amount present
            dew_evap = min(dble_one,(canopy_store/max_storage))*dew_evap
            ! convert back to kg.m-2.t-1 for remaining checks
            dew_mass = (dew_evap/lambda_bulk)*time%seconds_per_step
            ! if positive then canopy evaporation has occured instead
            ! and is evaporated from the canopy water supply
            canopy_store_old = canopy_store
            if (canopy_store > dew_mass) then
                ! loss of water from canopy store on given layer
                layer_store = dew_mass
            else
                layer_store = canopy_store
            end if
            if (rad_pass == max_pass) then
                ! restrict the canopy_store to positive values only
                canopy_store = max(dble_zero,canopy_store-layer_store)
                ! evap is the difference between new and old canopy
                ! stores (kg.m-2.t-1)
                layer_wet_evap = canopy_store_old-canopy_store
                wetev(time%step) = wetev(time%step)+layer_wet_evap
            else
                layer_wet_evap = canopy_store_old-max(dble_zero,canopy_store-layer_store)
            end if
            ! now convert to W.m-2 for use in energy balance
            dew_evap = (layer_wet_evap/time%seconds_per_step)*lambda_bulk
        else ! canopy_store > 0
            dew_evap = 0d0
        end if ! canopy_store > 0

    else ! if we have dew formation / wet evap has now been switched to dew

        ! potentially dew is falling, lets see if we can balance the
        ! books / realism

        ! limit dew fall interms of possible
        ! energy input; 300 W.m-2
        if (((dew_mass/time%seconds_per_step)*lambda_bulk) < -300d0) then
            dew_mass = (-300d0/lambda_bulk)*time%seconds_per_step
        end if
        if (transt_kg_s(time%step) < 0d0) then
            ! restrict canopy dew by dew may have occured in
            ! transpiration (this is technically a bug correction)
            dew_mass = min(dble_zero,dew_mass-(transt_kg_s(time%step)*time%seconds_per_step))
        end if
        ! so dew_mass is negative, dew is falling, does it reach
        ! limit of canopy?
        if ((-dew_mass) > max(dble_zero,max_storage-canopy_store)) then
            ! more dew fall than canopy can take so limit the amount
            ! of the canopy and assume that the rest goes direct to
            ! the soil surface
            dew_mass = -(max(dble_zero,max_storage-canopy_store))
        end if
        ! if negative then dew formation has occured and this is
        ! added to the water input for next timestep
        dew = dew+(-dew_mass)
        ! if dew then pass the whole volume to canopy evap variable
        ! (kg.m-2.t-1)
        wetev(time%step) = wetev(time%step)+dew_mass
        dew_evap = (dew_mass/time%seconds_per_step)*lambda_bulk ! (W.m-2) 

    end if ! dew_evap > 0.

  end subroutine dew_or_wet_canopy_evaporation
  !
  !----------------------------------------------------------------------
  !
  subroutine leaf_processes( i , leaf_fraction , leaf_area , sunlit)

    use gv_clim,               only: gbw, ppfd_to_par, temp_bot, temp_top, wdbot, wdtop, leaf_heat_conductance
    use gv_meteo,              only: gbb, gbh, layer_LCA, nit, par, rad, temp, wdef, dew_evap, live_frac
    use gv_scale_declarations, only: time, user_opts
    use gv_veg,                only: gppt, lai, LCA,Nla, respt, transt, transt_kg_s, lafrac!, Cfol_profile
    use leaf,                  only: assimilate, set_leaf
    use light,                 only: checkpar, longem, nirabs, parabs_shade, parabs_sun

    ! subroutine prepares final inputs to the leaf subroutines and deals with
    ! outputs. All analyses are calculated assuming 1 m2/m2 leaf area and scaled
    ! at the end of assimilate(). NOTE: that any changes here will need to be
    ! mirrored in "canopy_optimisation_functions.f90"

    implicit none

    ! arguments..
    integer,intent(in) :: i  ! index
    double precision,intent(in)    :: leaf_area,     & !
                          leaf_fraction    ! fraction of leaf exposed (to sun/shade, depending on call)
    logical,intent(in) :: sunlit

    ! local variables..
    logical :: marginal = .false.
    double precision :: agr, & !
                    dew_rad, & ! 
                        gsm, & !
                      modet, & !
                 modet_kg_s, & !
                        res    !

    ! mean nitrogen in currect leaf fraction per leaf area (gN.m-2 leaf)
    nit = Nla(i)
    ! load LCA for the specific layer (gC.m-2)
    layer_LCA = LCA(i) !Cfol_profile(i) / lai(i)

    !  leaf area has all beam rad plus frac of diffuse rad..
    if ( sunlit ) then
      par = ( parabs_sun(i) + parabs_shade(i) * leaf_fraction ) / leaf_area
    else
      par = parabs_shade(i) * leaf_fraction / leaf_area
    end if

    if (user_opts%iterative_canopy) then
       ! calcalate impact on energy balance of dew to this canopy layer (W.m-2)
       dew_rad = (dew_evap*live_frac*lafrac(i)*leaf_fraction) / leaf_area
    else
       dew_rad = 0d0
    end if
    ! net radiation (kW.m-2) = shortwave + longwave radiation
    rad = 0.001d0 * ( nirabs(i) * leaf_fraction / leaf_area     &
                     + longem(i) * leaf_fraction / leaf_area  &
                       + (par / ppfd_to_par) - dew_rad )
    temp = temp_top - 0.1111d0 * ( i - 1 ) * ( temp_top - temp_bot )
    wdef = wdtop - 0.1111d0 * ( i - 1 ) * ( wdtop - wdbot )
    gbb  = gbw(i)                    ! boundary layer for water vapour (m.s-1)
    gbh  = leaf_heat_conductance(i)  ! boundary layer for heat (m.s-1)

    call set_leaf(i,leaf_fraction)    ! see leaf.f90

    call assimilate(i,time,sunlit,modet,modet_kg_s,agr,res,gsm,marginal)    ! see leaf.f90

    ! update variables
    gppt(time%step)   = gppt(time%step)   + agr     ! umol CO2 assimilation m-2 ground area s-1
    respt(time%step)  = respt(time%step)  + res     ! umol CO2 respiration m-2 ground area s-1
    transt(time%step) = transt(time%step) + modet   ! evapotranspiration in W m-2 ground area
    transt_kg_s(time%step) = transt_kg_s(time%step) + modet_kg_s   ! evapotranspiration in kg.m-2.s-1 ground area
    checkpar(i) = checkpar(i) + par * leaf_area

  end subroutine leaf_processes
  !
  !----------------------------------------------------------------------
  !
  subroutine leaf_temperature_hfx( i , leaf_fraction , leaf_area, fraction_sensible , fraction_temperature , sunlit , dead )

    use gv_clim,               only: leaf_heat_conductance, gbw, ppfd_to_par, temp_bot, temp_top, wdbot, wdtop
    use gv_meteo,              only: gbb, gbh, par, rad, temp, wdef, live_frac, dew_evap
    use gv_scale_declarations, only: user_opts
    use gv_veg,                only: lafrac, canopy_rad, lafrac_dead
    use leaf,                  only: leaf_balance, leaf_temperature, gs2, leaf_temp_precision
    use light,                 only: longem, nirabs, parabs_shade, parabs_sun
    use math_tools,            only: zbrent

    implicit none

    ! arguments..
    integer,intent(in) :: i  ! index
    double precision,intent(in)    :: leaf_area,     & !
                          leaf_fraction    ! fraction of leaf exposed (to sun/shade, depending on call)
                                           ! when dead foliage leaf_fraction is = 1
    double precision,intent(inout) :: fraction_sensible, & ! sensible heat flux (W.m-2) for current leaf layer / fraction
                          fraction_temperature ! temperature (oC) for current leaf layer / fraction
    logical,intent(in) :: sunlit, dead

    ! local variables..
    double precision    :: dew_rad, & ! energy due to dew / wet evaporation on canopy (W.m-2)
               tmprad     ! temperary holder for canopy level net radiation (kW.m-2)

    !  leaf area has all beam rad plus frac of diffuse rad..
    ! calculate PAR (umol.m-2.s-1)
    if ( sunlit ) then
      par = ( parabs_sun(i) + parabs_shade(i) * leaf_fraction ) / leaf_area
    else if ( dead ) then
      ! we assume that we can calculate the whole dead foliage temperature in
      ! one go, hence why no multiply by leaf_fraction to PAR calculation. 
      ! However it may be appropriate to change this at a later date if we get more detailed
      par = ( parabs_sun(i) / leaf_area ) + ( parabs_shade(i) / leaf_area )
    else
      par = parabs_shade(i) * leaf_fraction / leaf_area
    end if

    ! calculate contribution to energy balance of dew / wet surface evap
    if (user_opts%iterative_canopy) then
       ! special case when dealing with dead foliage
       if ( dead ) then
          dew_rad = (dew_evap*(1d0-live_frac)*lafrac_dead(i)*leaf_fraction) / leaf_area 
       else
          ! calcalate impact on energy balance of dew to this canopy layer (W.m-2)
          dew_rad = (dew_evap*live_frac*lafrac(i)*leaf_fraction) / leaf_area
       end if
    else
       dew_rad = 0d0
    end if

    ! net radiation = shortwave + longwave radiation
    rad = 0.001d0 * ( nirabs(i) * leaf_fraction / leaf_area     &
                     + longem(i) * leaf_fraction / leaf_area  &
                       + (par / ppfd_to_par) - dew_rad)
    temp = temp_top - 0.1111d0 * dble( i - 1 ) * ( temp_top - temp_bot )
    wdef = wdtop - 0.1111d0 * dble( i - 1 ) * ( wdtop - wdbot )
    gbb  = gbw(i)                    ! boundary layer for water vapour (m.s-1)
    gbh  = leaf_heat_conductance(i)  ! boundary layer for heat (m.s-1)

    if ( sunlit ) then
        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
           fraction_temperature = leaf_temperature( gs2 )
        elseif (user_opts%solve_canopy_temperature == 2) then
           fraction_temperature = zbrent('leaf_temperature_hfx_sun:leaf_balance' &
                                         ,leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if
    else if ( dead ) then
        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
           fraction_temperature = leaf_temperature( gs2 )
        elseif (user_opts%solve_canopy_temperature == 2) then
           fraction_temperature = zbrent('leaf_temperature_hfx_dead:leaf_balance' &
                                         ,leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if
    else
        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
           fraction_temperature = leaf_temperature( gs2 )
        elseif (user_opts%solve_canopy_temperature == 2) then
           fraction_temperature = zbrent('leaf_temperature_hfx_shade:leaf_balance' &
                                         ,leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if
    end if

    ! on max_pass rad has been corrected to net radiation from isothermal, Short
    ! wave, longwave and dew / wet evap corrected to (kW.m-2)
    tmprad = (rad+(dew_rad*1d-3))
    ! Assign canopy layer net radiation converted from (kW.m-2) to (W.m-2) to
    ! accumulation variable for use in checking the canopy energy balance
    canopy_rad = canopy_rad+(tmprad*1d3*leaf_area)
    ! current area sensible heat (W.m-2)
    fraction_sensible = canopy_sensible(leaf_area,fraction_temperature,gbh,temp)

  end subroutine leaf_temperature_hfx
  !
  !----------------------------------------------------------------------
  !
  double precision function penman_monteith(numsecs,canopy_radiation,canopy_temperature,layer_conductance,canopy_wdef)

    ! Penman-Monteith equation for potential evaporation (mm.t-1 or kg.m-2.t-1)
    ! # 
    ! Used for evaporation and dew of water on canopy surface; aerodynamical
    ! resistance only, no root, stomatal effects etc    #

     use gv_meteo, only: abs_pot_conv, air_density_kg
     use gv_clim,  only: cp_air, lambda_bulk

     implicit none
     double precision ::    slope & ! Rate of change of saturation vapour pressure with temperature (Pa.K-1)
                             ,rho & ! Temperature dependant density of air (kg.m-3)
                           ,psych & ! Psychrometer constant (Pa.K-1); typical value ~66.1 Pa/K
                                    ! coefficient between air temperature and wet bulb temp, therefore determining potential latent energy exchange (?)
                               ,s & ! Straight line approximation of the true slope; used in determining relationship slope
                     ,canopy_wdef & ! canopy level water deficit (g.m-3)
                      ,canopy_vpd & ! canopy level vapour pressure deficit (Pa)
               ,layer_conductance & ! canopy layer conductance to exchange (m.s-1)
                ,canopy_radiation & ! canopy layer net radiation (W.m-2)
              ,canopy_temperature & ! Temperature of air (oC); potential
     ,canopy_temperature_absolute & ! Temperature of air (oC); absolute
                         ,numsecs   ! Number of seconds in timestep

     ! Determine density of air (kg.m-3); temperature dependant
     rho = air_density_kg
     ! convert canopy temperature to absolulte for use in water linked equations
     canopy_temperature_absolute = ((canopy_temperature+273.15d0)/abs_pot_conv)-273.15d0
     ! calculate canopy level VPD (Pa)
     canopy_vpd = canopy_wdef*(canopy_temperature_absolute+273.15d0)/2.17d0
     ! Straight line approximation of the true slope; used in determining
     ! relationship slope
     s = 6.1078d0*17.269d0*237.3d0*exp(17.269d0*canopy_temperature_absolute/(237.3d0+canopy_temperature_absolute))
     ! Slope of saturation vapour pressure curve; temperature depentent (Pa.K-1)
     slope = 100d0*(s/(237.3d0+canopy_temperature_absolute)**2)
     ! Psychrometer constant (Pa.K-1); coefficient between air temperature and
     ! wet bulb temp, therefore determining potential latent energy exchange 
     psych = (64.6d0*exp(0.00097d0*canopy_temperature_absolute))

     ! If evaporation takes place equation is positive, dew formation is
     ! negative
     penman_monteith = ((slope*canopy_radiation)+(rho*cp_air*canopy_vpd*layer_conductance))/(lambda_bulk*(slope+psych))
     ! Penman-Monteith equation, kgH2O m-2 s-1
     ! convert from kgH2O.m-2.s-1 to kgH2O.m-2.t-1
     penman_monteith = (penman_monteith*numsecs)

  end function penman_monteith
  !
  !----------------------------------------------------------------------
  !
end module canopy
!
!------------------------------------------------------------------------
!
 
