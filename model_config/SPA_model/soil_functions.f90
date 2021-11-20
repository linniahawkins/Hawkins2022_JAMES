! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module soil_functions

  !! SPA_DALEC version, incorporates calculation of field capacity   !!
  !! Crank Nicholson model to solve temperature flux in soil profile !!
  !! new latent energy flux model (based on SWP of surface layer)    !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: soil_processes,grav, vonkarman, soil_conductance, exchange_coefficient, lhs, soil_roughl

  ! physical constants used in module..
  ! (if want to duplicate them elsewhere then perhaps declare higher up)
  double precision,parameter :: grav      = 9.8067,   & ! acceleration due to gravity, m s-2
                              soil_roughl = 0.05,     & ! soil roughness length (m) values ~0.01-0.05
                                vonkarman = 0.41e0,   & ! von Karman's constant
                                Vw        = 18.05e-6, & ! partial molal volume of water, m3 mol-1 at 20C
                                hcap_ice  = 2100.0,   & ! specific heat capacity of ice, J K-1 kg-1
                                hcap_wat  = 4180.0,   & ! specific heat capacity of water, J K-1 kg-1
                                lhf       = 334000.0, & ! latent heat of fusion of ice, J kg-1
                                lhs       = 2.835e6,  & ! latent heat of sublimation of ice, J kg-1
                                rho_ice   = 917.0,    & ! density of ice, kg m-3
                                rho_wat   = 1000.0      ! density of water, kg m-3

  double precision :: soil_conductance & ! soil surface conductance for heat / water vapour (m.s-1)
                    ,snow_energy_stock & ! energy equivalent in the snow stock
                    ,infi             &
                    ,lambda_soil        ! latent heat of vapourisation (J kg-1) for the soil surface

  logical :: snow_scheme = .false.
  save

contains
  !
  !-------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_processes( time , totestevap )

    ! > subroutine summary? < !

    use gv_clim,                only: atmos_press, snowfall, ppt, temp_bot, temp_top, wdbot, wetev, wind_spd
    use gv_hourscale,           only: freeze, evap_store, hour, hourppt, hourpress, hourrad,  & 
                                      hourrnet, hourtemp, hourtime,hourts, hourvpd, hourwind, overflow,& 
                                      Qc, Qe, Qh, Qn, Qm, Qs, surface_watermm, totet, underflow, canopy_store
    use gv_irradiance_sunshade, only: soilnet
    use gv_scale_declarations,  only: grid, time_holder, user_opts
    use gv_snow_info,           only: snow_watermm, snowheight
    use gv_soil_structure,      only: conduc, pptgain, soil_temp, soil_temp_nplus1, thickness, waterfrac, watergain, watericemm, &
                                      waterloss, soilR, soilR1, soilR2, rooted_layers, root_length, root_mass
    use gv_veg,                 only: ess, lai, modrnet, soiletmm, dew
    use gv_daily_averages,      only: daily_ground_heat_flux
    use math_tools,             only: zbrent
    use soil_air,               only: slayer, soil_conductivity, soil_porosity, soil_resistance, soil_water_potential, &
                                      water_uptake_layer
    use spa_io,                 only: handle_output
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    double precision, intent(out)   :: totestevap

    ! local variables..
    integer :: i 
    double precision    :: checkdiff, checkerror, cw, delw, fluxsum, pw, t1, &
               t2, ts, tb, waterchange, xacc, snow_watermm_previous

    !--- main calculations begin below ---

    pw         = sum(watericemm)
    xacc       = 0.0001d0
    Qh         = 0d0
    Qe         = 0d0
    Qn         = 0d0
    Qc         = 0d0
    Qm         = 0d0
    Qs         = 0d0
    hour       = time%step
    i          = time%step
    hourtemp   = temp_bot + freeze      !convert from oC to K
    hourvpd    = 0.1d0 * wdbot * hourtemp / 217.0d0
    hourwind   = wind_spd
    hourrad    = soilnet
    hourpress  = atmos_press
    hourtime   = time%daytime
    hourppt    = ppt
    hourrnet   = modrnet
    snow_watermm_previous = snow_watermm
    snow_energy_stock = (snow_watermm / time%seconds_per_step)*lhs ! estimate potential Qs rate
    surface_watermm = 0d0      ! assume that water infiltrates within a timestep
    waterloss  = 0d0
    watergain  = 0d0
    pptgain    = 0d0

    ! Re-calculate soil water..
    do slayer = 1 , grid%core   ! loop over layers..
      conduc(slayer) = soil_conductivity( waterfrac(slayer) )
    enddo

    ! determine soil porosity and conductivity 
    call soil_porosity
    ! determine soil water potential based on type and soil moisture content 
    call soil_water_potential
    ! determine soil-root hydraulic resistance
    call soil_resistance(soilR,soilR1,soilR2,root_length,root_mass,rooted_layers)

    ! current assumption is that even if there is liquid rain (<0oC air
    ! temperature) that dew formation is infact to sublimation and thus frozen
    if ( temp_top <= 0d0) then
       if (snow_scheme) then
           ! if the snow scheme is active then add dew to snowfall for addition
           ! to that scheme (below)
           snowfall = snowfall + (dew / time%seconds_per_step)
       else
           ! if snow scheme not active then just accumulate snow in simple pool
           snow_watermm = snow_watermm + dew + (snowfall * time%seconds_per_step)
       endif
       dew = 0d0
    endif
  
    ! Assuming dew is infact liquid add direct to the canopy water pools.
    canopy_store = canopy_store + dew ; dew = 0d0

    if (.not.snow_scheme) then
       ! remove from the snow pack assuming constant decay, this is actually
       ! rubbish and should at least have a time step dependent component.
       if ( (snow_watermm .gt. 0d0 ) .and. ( temp_top .gt. 1d0 ) ) then ! melt snow
          ! snow decay rate
          snow_watermm    = snow_watermm * 0.9d0    ! decay rate
          surface_watermm = surface_watermm + snow_watermm * 0.1d0 ! how much water reaches the soil surface from melting snow?
       end if
       ! remove small amount of snow to avoid numeric solver problems
       if ( snow_watermm .lt. 0.01d0 ) snow_watermm = 0d0
    end if ! snow scheme?

    if ( ppt > 0d0 .or. canopy_store > 0d0) then ! there appears to be liquid rain irrespective of temperature
        call canopy_balance           ! interception of rain and canopy evap
        ! only do this if we are using the old, uncoupled wet surface evaporation
        if ( .not. user_opts%iterative_canopy) then
            wetev(i) = evap_store     ! mm t-1
        end if
    end if

    ! parameters for solving soil surface energy balance
    t1 = hourtemp - 50d0
    t2 = hourtemp + 50d0   ! set temp bounds based on air temp
    ts = zbrent( 'soil_processes:energy' , energy , t1 , t2 , xacc )    ! calculate surface temp based on energy balance
    if (allocated(daily_ground_heat_flux)) daily_ground_heat_flux(time%step) = Qc
    if (snow_scheme .and. snow_watermm > 0d0) then
        lambda_soil = lhs
    else if (ts <= freeze) then
        ! sublimation occurs, we assume latent heat of sublimation is (J kg-1)
        lambda_soil = lhs !2.835e+6
    else
        ! latent heat of vapourisation (J kg-1)
        lambda_soil = 1000d0 * ( 2501.0d0 - 2.364d0 * ( ts - freeze ) )
    endif
    if (snow_scheme) then
       ! if we have snow cover see if some of it will melt...
       if (snow_watermm > 0d0) call melt(time, ts)
       ! calculate snow temperature and mass profiles anyway
       call snow(time, snowfall, ts, tb)
    endif

    ! if we have leaves estimate where in the soil profile transpiration is coming from
    if ( sum( lai ) .gt. 0d0 ) then 
      call water_uptake_layer( totestevap )
    end if
    ! load soil surface temperature to memory
    hourts = ts

    ! 
    if (user_opts%prescribed_swc) then
        call prescribed_water_fluxes( time )
    else
        ! losses and gain terms for each soil layer and profile as a whole
        call water_fluxes( time ) 
        ! update temperature profile for water movement within bounds of the soil
        ! surface and core temperatures
        call water_thermal
    endif  
    ! update variables for soil temperature profile update
    soil_temp_nplus1       = 0d0
    soil_temp(1)           = ts
    soil_temp_nplus1(1)    = ts
    soil_temp_nplus1(grid%core) = soil_temp(grid%core)
    ! update soil temperature profile
    call crank_nicholson   
    ! determine thaw depth and icefractions for each soil layer
    call thaw_depth(ts) 

    ! load soil surface evaporation and snow sublimation into memory (W.m-2)
    ! NOTE: sign convention changed such that + if evaporating. 
    ess(time%step)        = (-Qe) + (-Qs)
    soiletmm( time%step ) = ess( time%step ) / lambda_soil * time%seconds_per_step  ! convert to mm t-1

    ! net change in profile soil water content
    waterchange = 1d3 * ( sum(pptgain) + sum(watergain) - watergain(grid%core) ) &
                 - 1d3 * ( sum(waterloss) - waterloss(grid%core) ) + (snow_watermm_previous-snow_watermm)
    ! sum flux into or out of the system based on flux terms
    fluxsum    = soiletmm(time%step) + 1d3 * ( overflow + totet + underflow ) - surface_watermm
    ! snow pack correction
    fluxsum    = fluxsum - ( (Qm / lambda_soil) * time%seconds_per_step)
    checkerror = waterchange + fluxsum
    ! surrent water content
    cw         = sum( watericemm )
    ! difference between current and previous water contents
    delw       = pw - cw
    ! the difference between current and previous time step should equal the
    ! calculated flux terms, if not then we have a leak
    checkdiff  = delw - fluxsum 

    ! report errors in state water balance
!    if (abs(checkdiff) > 0.01 .and. time%steps_count > 1) then
!        write(message,*) 'Soil water fluxes and state error; imbalance ',checkdiff
!        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
!    end if

    ! write output..
    if (user_opts%soils_csv_output) then
        ! output per model step
        call handle_output( 6 , output_data=(/ ts, waterchange, &
                               fluxsum, checkerror, delw, checkdiff /) )
        ! output per day
        if ( time%step .eq. time%steps_per_day )  &
            call handle_output( 7 , output_data=(/ ts, waterchange, &
                                 fluxsum, checkerror, delw, checkdiff /) )
    endif ! 

  end subroutine soil_processes
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine canopy_balance

    ! > subroutine summary? < !

    use gv_hourscale,      only: canopy_store, surface_watermm, hourppt
    use gv_soil_structure, only: max_storage, through_fall
    use math_tools,        only: dxsav, kmax, ode_int
    use log_tools

    implicit none

    ! local variables..
    integer,parameter :: nvar = 2
    integer           :: nbad, nok
    double precision  :: eps, h1, hmin, x1, x2, ystart(nvar)

    eps  = 1.0d-3 ! accuracy of integrator
    h1   = 0.001d0  ! first guess step size 
    hmin = 0d0    ! minimum step size allowed in integration
    kmax = 100    ! maximum number of iterations?
    x1   = 1d0     ! x1 and x2 define the integrator range for the process to occur; i.e. 2-1 = 1, 1 step for integrator
    x2   = 2d0
    dxsav = ( x2 - x1 ) * 0.05d0
    ! initial conditions
    if ( canopy_store .lt. 1d-5 * max_storage ) canopy_store = 0d0   ! empty store if it's tiny
    ystart( 1 ) = canopy_store  
    ystart( 2 ) = surface_watermm

    call ode_int( 'canopy_balance:canopy_water_store' , ystart , nvar , x1 , &
                   x2 , eps , h1 , hmin , nok , nbad , canopy_water_store )

    ! error condition for if water is disappearing in the process
    if ( abs(((ystart(1) + ystart(2)) - (canopy_store+surface_watermm)) - hourppt) > 0.1d0) then
        write(message,*) 'Input water to surface / canopy does not balance within tolerance storage and losses (0.1 mm)'
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*) 'Precip input (mm.t-1)', hourppt
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'water balance (mm)', ( (ystart(1) + ystart(2)) - (canopy_store+surface_watermm) ) - hourppt
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'canopy store (mm) start', canopy_store
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'surface water (mm) start', surface_watermm
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'canopy store (mm) after', ystart(1)
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'surface water (mm) after', ystart(2)
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'current canopy store', canopy_store
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
        write(message,*)'through_fall', through_fall
        call write_log( message  , msg_warning , __FILE__ , __LINE__ )
    end if

    canopy_store    = ystart( 1 )
    surface_watermm = ystart( 2 )

  end subroutine canopy_balance
  !
  !----------------------------------------------------------------------
  !
  subroutine canopy_water_store( time_dummy, y , dydt )

    ! determines canopy water storage and evaporation, !
    ! and water reaching soil surface.                 !

    use gv_hourscale,          only: evap_store, hourppt
    use gv_scale_declarations, only: max_nos_iterations, time, user_opts
    use gv_soil_structure,     only: max_storage, through_fall
    use log_tools

    implicit none

    ! arguments..
    double precision,intent(in)  :: y(max_nos_iterations)
    double precision,intent(in)  :: time_dummy ! dummy argument, provided for ode_int
    double precision,intent(out) :: dydt(max_nos_iterations)

    ! local variables..
    double precision :: a, add_ground, add_store, b, drain_store, potential_evap, ratio, drainage_coefficient

    add_store  = (( 1d0 - through_fall ) * hourppt) ! rate of input of water to canopy storage per min
    add_ground = through_fall * hourppt             ! rate of input of water to ground

    ! if not iterating canopy we must use the old uncoupled (incorrect) function
    ! to get wet canopy evap. Plus there is no point in linking drainage
    ! coefficients to LAI, via the max_storage if not iterating either
    if (user_opts%iterative_canopy) then 
       ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
       ! 0.002 is canopy specific coefficient modified by
       ! 0.002*(max_storage/1.05)
       ! where max_storage is the canopies maximum capacity (mm) (LAI based) and
       ! 1.05 is the original canopy capacitance
       drainage_coefficient = 0.002d0 * ( max_storage / 1.05d0 )
       ! zero the evap store when iterative canopy in use as the dew / wet
       ! surface evap is carried out directly in the canopy.f90
       evap_store = 0d0
    else ! use uncoupled
       ! function gives potential evaporation rate..
       potential_evap = wetted_surface_evap()
       ! rate of evaporation from storage NB store cannot exceed max storage..
       ratio      = min( 1d0 , y( 1 ) / max_storage )
       evap_store = potential_evap * ratio
       drainage_coefficient = 0.002d0
    end if

    b = 3.7d0
    a = log( drainage_coefficient ) - b * max_storage
    if ( y( 1 ) .gt. max_storage ) then
      ! rate of drainage from store (mm per timestep (in units of minutes))
      drain_store = exp( a + b * y(1) ) * ( time%seconds_per_step / 60d0 )
    else
      drain_store = 0d0
    end if

    if ( ( y(2) .gt. 0d0 ) .and. ( y(2) .lt. 7d-37 ) ) then
      write(message,*) 'problem in Canopy_water_store: y(2) (',y(2),') is too small.'
      call write_log( message , msg_warning , __FILE__ , __LINE__ )
    end if

    dydt(1) = add_store - drain_store - evap_store  ! change in canopy storage
    dydt(2) = add_ground + drain_store              ! addition to soilwater
    if ( dydt(1) .gt. 100d0 ) then
      write(message,*) 'Problem in Canopy_water_store: rate of change of y (',dydt(1), &
                       ') is too large to be physically plausible.'
      call write_log( message  , msg_warning , __FILE__ , __LINE__ )
    end if

  end subroutine canopy_water_store
  !
  !----------------------------------------------------------------------
  !
  subroutine crank_nicholson()
 
    ! Finite difference Partial Differential Equation solver for soil temperature
    ! profile. Crank Nicholson equations, determines temperature profile using surface
    ! temperature and constant core temperature to interpolate between.
    ! Process iterates: runs until the difference between adjacent soil layers throughout the profile
    ! changes less than the max_error value      
 
    use gv_scale_declarations, only: grid, time
    use gv_soil_structure,     only: soil_temp, soil_temp_nplus1, thermal, thickness

    implicit none

    ! local variables..
    integer :: i
    double precision    :: beta, D, error, max_error, old_value, tdiffuse

    ! --calculations begin below--
    max_error = 0.0000005d0
    beta      = 0.5d0

    ! Continue looping until error is sufficiently reduced
    ! (exit condition at end of loop)
    do
      error = 0d0  ! reset error

      ! (Loop over all x-dimension nodes, except first and last)
      do i = 2 , grid%core - 1 
        thermal = thermal_conductivity( i )
        ! Walbroeck thermal conductivity, W m-1 K-1 is converted to J m-1 K-1 timestep-1...
        tdiffuse = time%seconds_per_step * thermal / heat_capacity( i )
        D = tdiffuse / ( thickness( i ) * thickness( i ) )
        old_value = soil_temp_nplus1( i )          ! Store value of previous iteration
        ! Calculate the temperature at the new time step using an implicit method..
        soil_temp_nplus1( i ) = ( D / ( 1d0 + 2d0 * beta * D ) ) &
                             * ( beta * ( soil_temp_nplus1( i + 1 ) + soil_temp_nplus1( i - 1 ) ) &
                             +  ( 1 - beta ) * ( soil_temp( i + 1 ) - 2d0 * soil_temp( i ) + soil_temp( i - 1 ) ) ) &
                             +  soil_temp( i ) / ( 1d0 + 2d0 * beta * D )

        error = error + abs( old_value - soil_temp_nplus1( i ) )  ! Calculate the error

      enddo

      ! exit if total error sufficiently small..
      if ( error .le. max_error ) exit

    enddo

    ! Set the values at time n equal to the values at time n+1 for the next time step..

    ! (Loop over all x-dimension nodes, except first and last)
    do i = 2 , grid%core - 1  
      soil_temp( i ) = soil_temp_nplus1( i )
    enddo

  end subroutine crank_nicholson
  !
  !----------------------------------------------------------------------
  !
  double precision function energy( ts )

    ! Determines the surface energy balance. !
    ! (not canopy)                           !
    ! (see Hinzmann et al. 1998)             !

    use gv_hourscale,          only: freeze, hourrad, hourtemp, Qc, Qe, Qh, Qn, Qs
    use gv_scale_declarations, only: boltz
    use gv_soil_structure,     only: soil_temp, thickness, waterfrac
    use gv_clim,               only: cp_air
    use gv_meteo,              only: air_density_kg
    use gv_veg,                only: emiss
    use gv_snow_info,          only: Dsnow, Tsnow, snow_watermm, snowheight, fsnow

    implicit none

    ! arguments..
    double precision,intent(in) :: ts ! incoming soil surface temperature (K)

    ! local variables..
    double precision ::  downwelling_rad, gah, rho, upwelling_rad, &
             Dtop, Ktop, Ttop
    infi = 0d0

    ! Sensible heat flux..
    rho = air_density_kg                  ! density of air kg m-3 (t-dependent)
    gah = soil_conductance                ! load soil conductance to heat from memory (m.s-1)
!    call exchange_coefficient( gah )      ! conductance to heat, m s-1
 
    ! sensible heat flux (W.m-2); postive if surface absorbing heat
    Qh  = cp_air * rho * gah * ( hourtemp - ts )

    if (snow_scheme .and. snow_watermm > 0d0) then
        lambda_soil = lhs
    else if (ts <= freeze) then
        ! sublimation occurs, we assume latent heat of sublimation is (J kg-1)
        lambda_soil = lhs 
    else
        ! latent heat of vapourisation (J kg-1)
        lambda_soil = 1000d0 * ( 2501.0d0 - 2.364d0 * ( ts - freeze ) )
    endif

    ! calculate soil surface evaporation or snow sublimation
    if (snow_scheme) then
       if (snow_watermm > 0d0) then
          ! sublimation occurs, we assume latent heat of sublimation is (J kg-1)
          Qe = qe_flux( ts )*(1d0-fsnow) !  = 0.0
          Qs = qs_flux( ts ) ! qs_flux( ts )
          Qs = max(-snow_energy_stock,Qs)*fsnow
       else
          Qs = 0d0
          Qe = qe_flux( ts )
       endif
    else 
       ! Latent energy flux (W.m-2); negative when evaporating
       Qe = qe_flux( ts )
       Qs = 0d0
    endif

    ! Net radiation (emitted LW varies with surface temp)..
    ! The brackets stop the enormous ts value from swamping
    ! the much smaller emissivity and boltzman values.
    upwelling_rad = ( 1d0 * emiss * boltz )*ts**4
    ! PAR+NIR+LW incident on soil from light.f90
    downwelling_rad = hourrad
    ! net radiation of the soil surface (W.m-2); positive when gaining energy
    Qn = downwelling_rad - upwelling_rad

    if (snow_scheme) then
       ! Thermal conductivity of top snow or soil layer..
       Dtop = max(thickness(1), Dsnow(1))
       Ktop = thickness(1) / ( 2d0*Dsnow(1)/snow_conductivity(1) +  &
                              (thickness(1) - 2d0*Dsnow(1))/thermal_conductivity(1) )
       Ttop = soil_temp(2) + (Tsnow(1) - soil_temp(2))*Dsnow(1)/thickness(1)
       if ( Dsnow(1) .gt. 0.5d0*thickness(1) ) Ktop = snow_conductivity(1)
       if ( Dsnow(1) .gt. thickness(1) ) Ttop = Tsnow(1)
       Qc = -Ktop * ( ts - Ttop ) / ( 0.5d0 * Dtop )
    else
       ! (Hillel) Thermal conductivity of top layer, to calculate the ground heat
       ! flux between 1st and 2nd soil layers (strictly speaking this should be
       ! between skin and 1st soil); positive when heat moving up the soil profile
       Qc = -thermal_conductivity(1) * ( ts - soil_temp(2) ) / ( 0.5d0 * thickness(1) )
    endif

    ! Energy balance..
    energy = Qh + Qe + Qs + Qn + Qc
if (energy /= energy .or. Qs /= Qs .or. Qh /= Qh .or. Qe /= Qe .or. Qn /= Qn .or. Qc /= Qc ) then
print*,"energy",energy,"Qh",Qh,"Qe",Qe,"Qn",Qn,"Qc",Qc,"Qs",Qs
print*,"ts",ts,"rho",rho,"gah",gah

stop
end if
  end function energy
  !
  !----------------------------------------------------------------------
  !
  subroutine melt( time, ts )

    ! Determines snowmelt from the surface energy balance. !

    use gv_hourscale,          only: freeze, hourrad, hourtemp, Qc, Qe, Qh, Qn, &
                                     Qm, Qs
    use gv_scale_declarations, only: boltz, time_holder
    use gv_snow_info,          only: Dsnow, Tsnow, snow_watermm, snowheight, fsnow
    use gv_soil_structure,     only: soil_temp, thickness, waterfrac
    use gv_clim,               only: cp_air
    use gv_meteo,              only: air_density_kg
    use gv_veg,                only: emiss

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    double precision,intent(inout) :: ts

    ! local variables..
    double precision ::  downwelling_rad, gah, rho, upwelling_rad
    double precision ::  dSWE, dt, Dtop, Ktop, Ttop

    if (ts .gt. freeze) then

      ! ts = freeze ! orginal

      ! estimate exposed soil surface latent flux  
      Qe = qe_flux( ts ) * (1d0-fsnow)  ! added
      ! estimate snow covered latent flux
      Qs = qs_flux( freeze )                ! added
      Qs = max(-snow_energy_stock,Qs)*fsnow ! added

      ! update surface temperature based on snow fraction
      ts = (freeze*fsnow) + ( ts*(1d0-fsnow) )
      ! Sensible and latent heat fluxes..
      rho = air_density_kg                  ! density of air kg m-3 (t-dependent)
      gah = soil_conductance                ! load soil conductance to heat from memory (m.s-1) 
      Qh  = cp_air * rho * gah * ( hourtemp - ts )
!      Qs  = qs_flux( ts )

      ! Net radiation
      upwelling_rad = ( 1d0 * emiss * boltz ) * ts**4
      downwelling_rad =  hourrad
      Qn = downwelling_rad - upwelling_rad

      ! Thermal conductivity of top snow layer..
      Dtop = max(thickness(1), Dsnow(1))
      Ktop = thickness(1) / ( 2d0*Dsnow(1)/snow_conductivity(1) +  &
                              (thickness(1) - 2d0*Dsnow(1))/thermal_conductivity(1) )
      Ttop = soil_temp(2) + (Tsnow(1) - soil_temp(2))*Dsnow(1)/thickness(1)
      if ( Dsnow(1) .gt. 0.5d0*thickness(1) ) Ktop = snow_conductivity(1)
      if ( Dsnow(1) .gt. thickness(1) ) Ttop = Tsnow(1)
      Qc = -Ktop * ( ts - Ttop ) / ( 0.5d0 * Dtop )

      ! Melt from energy balance...positive when melting occurs
      Qm = Qh + Qs + Qn + Qc + Qe ! TLS added +Qe
      dt = time%seconds_per_step ! sign flipped for Qs here as Qs is negative when water is being removed
      dSWE = (Qm/lhf - Qs/lhs) * dt
      ! if snow melt is greater than stock then limit to total
      if (dSWE > snow_watermm) then
         ! to preseve the energy balance place the imbalance here to sensible heat
         Qh = Qh - ( Qm - (lhf * (snow_watermm/dt + Qs/lhs)) )
         ! then recalculate melting energy
         Qm = lhf * (snow_watermm/dt + Qs/lhs)
      endif

    end if

  end subroutine melt
  !
  !----------------------------------------------------------------------
  !
  subroutine exchange_coefficient( exco_out )

    ! heat or vapour exchange coefficient/boundary layer !
    ! cond, m s-1.   See Mat William's ref 1028          !

    use gv_clim,      only: wind_spd
    use gv_veg,       only: canopy_height, tower_height

    implicit none

    ! arguments..
    double precision, intent(out) :: exco_out

    ! local variables..
    double precision :: answer, log_frac, numerator, bulk_roughl

    ! Boundary layer conductance at ground level for NEUTRAL conditions.
    ! Substitute lower altitude than canopy_height
    ! Heat exchange coefficient from Hinzmann
    ! 0.13 * canopy_height gives roughness length
    ! NOTE: 0.13 is coefficient for bulk surface conductance
    ! with moderately dense canopy. Soil roughl itself should
    ! be ~ 0.01->0.05 m

    bulk_roughl = 0.01d0
    numerator   = 1d0 * wind_spd * ( 1d0 * vonkarman * vonkarman)
    log_frac    = log( tower_height / ( bulk_roughl * canopy_height ) )
    answer      = ( 1d0 * numerator ) / ( 1d0 * log_frac * log_frac)
    exco_out    = 1d0 * answer

  end subroutine exchange_coefficient
  !
  !----------------------------------------------------------------------
  !
  subroutine snow( time, snow_fall, ts, tb )

    ! finite difference PDE solver for snow temperature profile !

    use gv_scale_declarations,  only: grid, time_holder
    use gv_hourscale,           only: freeze, hourtemp, Qc, Qm, Qs, surface_watermm
    use gv_snow_info,           only: Dsfix, Dsnow, Nsnow, Sice, Sliq, snowheight,  &
                                      snowalb_nir, snowalb_par, snow_watermm, &
                                      snowweight,Tsnow,new_snownirref,new_snowparref
    use gv_soil_structure,      only: soil_temp, thickness
    use math_tools,             only: tridiag

    implicit none

    ! arguments..
      type(time_holder),intent(in) :: time
      double precision, intent(in)  :: snow_fall, ts
      double precision, intent(out) :: tb

    ! local variables
    double precision :: a(grid%snow),   & ! Below-diagonal matrix elements
            b(grid%snow),   & ! Diagonal matrix elements
            c(grid%snow),   & ! Above-diagonal matrix elements
        csnow(grid%snow),   & ! Areal heat capacity of snow (J/K/m^2)
          dTs(grid%snow),   & ! Temperature increments (k)
         D0(0:grid%snow),   & ! Layer thickness before adjustment (m)
          E(0:grid%snow),   & ! Energy contents before adjustment (J/m^2)
           Gs(grid%snow),   & ! Thermal conductivity between layers (W/m^2/k)
        ksnow(grid%snow),   & ! Thermal conductivity of snow (W/m/K)
        Dsnew(grid%snow),   & ! Available thickness in new layer (m)
          rhs(grid%snow),   & ! Matrix equation rhs
          S(0:grid%snow),   & ! Ice contents before adjustment (kg/m^2)
            U(grid%snow),   & ! Layer internal energy contents (J/m^2)
          W(0:grid%snow)      ! Liquid contents before adjustment (kg/m^2)

    double precision :: alim,           & ! Limiting albedo
        coldcont,           & ! Layer cold content (J/m^2)
              dl,           & ! Local snow depth (m)
              dt,           & ! Timestep (s)
           dSice,           & ! Change in layer ice content (kg/m^2)
           Dsold,           & ! Remaining thickness in old layer (m)
             phi,           & ! Porosity
            rhos,           & ! Density of snow layer (kg/m^3)
              rt,           & ! Reciprocal timescale for albedo adjustment (1/s)
           Sice0,           & ! Ice content of fresh snow (kg/m^2)
         SliqMax,           & ! Maximum liquid content for layer (kg/m^2)
             tau,           & ! Timescale (s)
          Tsnow0,           & ! Temperature of fresh snow (K)
             Win,           & ! Liquid water entering a layer (kg/m^2)
              wt              ! Layer weighting
    integer :: k,kold,knew,kstart, & ! Level pointers
               Nold                  ! Previous number of snow layers

    ! initial values
    dt = time%seconds_per_step
    tb = soil_temp(1)

    if (Nsnow > 0) then   ! Existing snowpack

        ! Snow albedo; note that 0.73 is nir reflectanc of new snow and 0.95 is
        ! par reflectance of new snow
        tau = 3600d0*100d0
        rt = 1d0/tau + snow_fall/10d0
        alim = (0.56d0/tau + snow_fall*new_snownirref/10d0)/rt
        if (ts >= freeze) alim = (0.38d0/tau + snow_fall*new_snownirref/10d0)/rt
        snowalb_nir = alim + (snowalb_par - alim)*exp(-rt*dt)
        alim = (0.84d0/tau + snow_fall*new_snowparref/10d0)/rt
        if (ts >= freeze) alim = (0.62d0/tau + snow_fall*new_snowparref/10d0)/rt
        snowalb_par = alim + (snowalb_par - alim)*exp(-rt*dt)

        ! Heat capacity and thermal conductivity
        do k = 1, Nsnow
           csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
           ksnow(k) = snow_conductivity(k)
        end do

        ! Heat conduction
        if (Nsnow == 1) then
            Gs(1) = 2d0 / (Dsnow(1)/ksnow(1) + thickness(1)/thermal_conductivity(1))
            dTs(1) = (Gs(1)*(soil_temp(2) - Tsnow(1)) - Qc)*dt /  &
                     (csnow(1) + Gs(1)*dt)
        else
            do k = 1, Nsnow - 1
               Gs(k) = 2d0 / (Dsnow(k)/ksnow(k) + Dsnow(k+1)/ksnow(k+1))
            end do
            a(1) = 0d0
            b(1) = csnow(1) + Gs(1)*dt
            c(1) = - Gs(1)*dt
            rhs(1) = - (Qc + Gs(1)*(Tsnow(1) - Tsnow(2)))*dt
            do k = 2, Nsnow - 1
               a(k) = c(k-1)
               b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
               c(k) = - Gs(k)*dt
               rhs(k) = Gs(k-1)*(Tsnow(k-1) - Tsnow(k))*dt  &
                        + Gs(k)*(Tsnow(k+1) - Tsnow(k))*dt
            end do
            k = Nsnow
            Gs(k) = 2d0 / (Dsnow(k)/ksnow(k) + thickness(1)/thermal_conductivity(1))
            a(k) = c(k-1)
            b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
            c(k) = 0d0
            rhs(k) = Gs(k-1)*(Tsnow(k-1) - Tsnow(k))*dt  &
                     + Gs(k)*(soil_temp(2) - Tsnow(k))*dt
            call tridiag( Nsnow, grid%snow, a, b, c, rhs, dTs)
        end if
        do k = 1, Nsnow
           Tsnow(k) = Tsnow(k) + dTs(k)
        end do
        k = Nsnow
        tb = (ksnow(k)*thickness(1)*Tsnow(k) + thermal_conductivity(1)*Dsnow(k)*soil_temp(2))  &
             / (ksnow(k)*thickness(1) + thermal_conductivity(1)*Dsnow(k))

        ! Convert melting ice to liquid water
        dSice = (Qm / lhf) * dt
        do k = 1, Nsnow
           coldcont = csnow(k)*(freeze - Tsnow(k))
           if (coldcont < 0d0) then
               dSice = dSice - coldcont / lhf
               Tsnow(k) = freeze
           end if
           if (dSice > 0d0) then
               if (dSice > Sice(k)) then  ! Layer melts completely
                   dSice = dSice - Sice(k)
                   Dsnow(k) = 0d0
                   Sliq(k) = Sliq(k) + Sice(k)
                   Sice(k) = 0d0
               else                       ! Layer melts partially
                   Dsnow(k) = (1d0 - dSice/Sice(k))*Dsnow(k)
                   Sice(k) = Sice(k) - dSice
                   Sliq(k) = Sliq(k) + dSice
                   dSice = 0d0              ! Melt exhausted
             end if
           end if
        end do

        ! Remove snow by sublimation 
!        dSice = max(Qs/lhs, 0.)*dt
        dSice = max(-Qs/lhs,0d0)*dt ! TLS: sign change
        if (dSice > 0d0) then
            do k = 1, Nsnow
               if (dSice > Sice(k)) then  ! Layer sublimates completely
                   dSice = dSice - Sice(k)
                   Dsnow(k) = 0d0
                   Sice(k) = 0d0
               else                       ! Layer sublimates partially
                   Dsnow(k) = (1d0 - dSice/Sice(k))*Dsnow(k)
                   Sice(k) = Sice(k) - dSice
                   dSice = 0d0                ! Sublimation exhausted
               end if
            end do
        end if

        ! Snow hydraulics
        Win = surface_watermm
        do k = 1, Nsnow
           phi = 0d0
           if (Dsnow(k) > epsilon(Dsnow)) phi = 1d0 - Sice(k)/(rho_ice*Dsnow(k))
           SliqMax = rho_wat*Dsnow(k)*phi*0.03d0
           Sliq(k) = Sliq(k) + Win
           Win = 0d0
           if (Sliq(k) > SliqMax) then  ! Liquid capacity exceeded
               Win = Sliq(k) - SliqMax    ! so drainage to next layer
               Sliq(k) = SliqMax
           end if
           coldcont = csnow(k)*(freeze - Tsnow(k))
           if (coldcont > 0d0) then       ! Liquid can freeze
               dSice = min(Sliq(k), coldcont/lhf)
               Sliq(k) = Sliq(k) - dSice
               Sice(k) = Sice(k) + dSice
               Tsnow(k) = Tsnow(k) + lhf*dSice/csnow(k)
           end if
        end do
        surface_watermm = Win

        ! Snow compaction
        tau = 3600d0*200d0
        do k = 1, Nsnow
           rhos = 100d0
           if (Dsnow(k) > epsilon(Dsnow)) then
               rhos = (Sice(k) + Sliq(k)) / Dsnow(k)
               if (Tsnow(k) >= freeze) then
                   if (rhos < 500d0) rhos = 500d0 + (rhos - 500d0)*exp(-dt/tau)
               else
                   if (rhos < 300d0) rhos = 300d0 + (rhos - 300d0)*exp(-dt/tau)
               end if
               Dsnow(k) = (Sice(k) + Sliq(k)) / rhos
           end if
        end do

    end if  ! Existing snowpack

    ! Add snowfall and frost as layer 0
    ! TLS: sign change for Qs
    Sice0 = max( snow_fall*dt - min(-Qs/lhs, 0d0)*dt, 0d0 )
    Tsnow0 = min(hourtemp, freeze)
    D0(0) = Sice0 / 100d0
    E(0) = Sice0*hcap_ice*(Tsnow0 - freeze)
    S(0) = Sice0
    W(0) = 0d0

    ! Calculate new snow depth
    snowheight = D0(0)
    do k = 1, Nsnow
       snowheight = snowheight + Dsnow(k)
    end do

    ! Store state of old layers
    do k = 1, Nsnow
       csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
       D0(k) = Dsnow(k)
       E(k) = csnow(k)*(Tsnow(k) - freeze)
       S(k) = Sice(k)
       W(k) = Sliq(k)
    end do
    Nold = Nsnow

    ! Initialise new layers
    Dsnow(:) = 0d0
    Sice(:) = 0d0
    Sliq(:) = 0d0
    Tsnow(:) = freeze
    U(:) = 0d0
    Nsnow = 0

    if (snowheight > 0d0) then  ! Existing or new snowpack

        ! Re-assign and count snow layers
        dl = snowheight
        Dsnow(1) = dl
        k = 1
        if (Dsnow(1) > Dsfix(1)) then
            do k = 1, grid%snow
               Dsnow(k) = Dsfix(k)
               dl = dl - Dsfix(k)
               if (dl <= Dsfix(k) .or. k == grid%snow) then
                   Dsnow(k) = Dsnow(k) + dl
                   exit
               end if
            end do
        end if
        Nsnow = k
        Dsnew(:) = Dsnow(:)

        ! Fill new layers from the top downwards
        knew = 1
        do kold = 0, Nold                     ! Loop over old layers
           Dsold = D0(kold)
           kstart = knew
           do k = kstart, Nsnow                ! Loop over new layers with remaining space
              if (Dsold > Dsnew(k)) then  ! New layer filled
                  Dsold = Dsold - Dsnew(k)
                  if (D0(kold) > epsilon(D0)) then
                      wt = Dsnew(k) / D0(kold)
                      Sice(k) = Sice(k) + S(kold)*wt
                      Sliq(k) = Sliq(k) + W(kold)*wt
                      U(k) = U(k) + E(kold)*wt
                  end if
                  knew = k + 1                    ! Move pointer to next new layer
              else                              ! Old layer will be exhausted by this increment
                  Dsnew(k) = Dsnew(k) - Dsold
                  wt = 1
                  if (D0(kold) > epsilon(D0)) wt = Dsold / D0(kold)
                  Sice(k) = Sice(k) + S(kold)*wt
                  Sliq(k) = Sliq(k) + W(kold)*wt
                  U(k) = U(k) + E(kold)*wt
                  exit  ! Proceed to next old layer by exiting new layer loop
               end if   
           end do       ! New layers
        end do          ! Old layers

    end if  ! Existing or new snowpack

    ! Diagnose snow layer temperatures and bulk SWE
    snow_watermm = 0d0
    do k = 1, Nsnow
       csnow(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
       if (csnow(k) > epsilon(csnow)) Tsnow(k) = freeze + U(k) / csnow(k)
       snow_watermm  = snow_watermm  + Sice(k) + Sliq(k)
    end do
    snowweight = snow_watermm

  end subroutine snow
  !
  !----------------------------------------------------------------------
  !
  double precision function heat_capacity( i )

    ! Walbroeck: impacts of freeze/thaw on energy fluxes !
    ! heat capacity, J m-3 K-1                           !

    use gv_hourscale, only: freeze
    use gv_soil_structure, only: iceprop, mineralfrac, organicfrac, soil_temp, thickness, waterfrac

    implicit none

    ! arguments..
    integer,intent(in) :: i

    ! local variables..
    double precision    :: delt, lw, volhc

    delt = 1d0         ! temperature range over which freezing occurs

    volhc = volumetric_heat_capacity(mineralfrac(i),organicfrac(i) &
                                    ,iceprop(i),waterfrac(i))

    ! soil_temp is in kelvin, so convert to celcius and test if it lies in freezing range (-1<x<0)..
    if ( ( (soil_temp(i)-freeze) .le. 0d0 ) .and. ( (soil_temp(i)-freeze) .gt. -delt ) ) then
       lw      = 1000d0 * waterfrac(i) * ( thickness(i) * 10d0 )       ! liquid water content (kg m-3 soil)
       heat_capacity = volhc + lhf * lw / delt
    else
       heat_capacity = volhc
    end if

  end function heat_capacity
  !
  !----------------------------------------------------------------------
  !
  subroutine infiltrate

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    use gv_hourscale,          only: runoff, surface_watermm, overflow
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: porosity, pptgain, thickness, waterfrac, watergain, waterloss

    implicit none

    integer :: i
    double precision    :: add   & ! surface water available for infiltration (m)
              ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert surface water from mm -> m
    add = surface_watermm * 0.001d0
    pptgain = 0d0
    do i = 1 , grid%soil
       ! determine the available pore space in current soil layer
       wdiff = max(0d0,(porosity(i)-waterfrac(i))*thickness(i)-watergain(i)+waterloss(i))

       if(add .gt. wdiff)then
          ! if so fill and subtract from input and move on to the next layer
          pptgain(i) = wdiff
          add = add-wdiff
       else
          ! otherwise infiltate all in the current layer
          pptgain(i) = add
          add = 0d0
       end if
       ! if we have added all available water we are done
       if(add .le. 0d0) exit
    enddo

    ! if after all of this we have some water left assume it is runoff
    if(add .gt. 0d0)then
       overflow = add
    else
       overflow = 0d0
    end if
    runoff = runoff+overflow

  end subroutine infiltrate
  !
  !----------------------------------------------------------------------
  !
  double precision function qe_flux( ts )

    ! latent energy loss from soil surface !

    use gv_hourscale,      only: freeze, gaw, gws, hourpress, hourtemp, hourvpd
    use gv_meteo,          only: Rcon, air_density_kg
    use gv_soil_structure, only: drythick, porosity, SWP, thickness

    implicit none

    ! arguments..
    double precision, intent(in) :: ts

    ! local variables..
    double precision :: diff, ea, esat_air,esat_soil, esurf, por, rho, tort

    tort = 1d0       ! tortuosity
    por = 0.9d0
    !por  = porosity(1) ! porosity of surface layer

!    call exchange_coefficient( gaw )                      ! determine boundary layer conductance (m s-1)
    gaw = soil_conductance                                 ! load soil conductance from memory (m.s-1)
    rho = air_density_kg                                   ! density of air (kg m-3) (t-dependent)

    diff = 24.2d-6 * ( ts / 293.2d0 )**1.75d0                  ! (m2 s-1) diffusion coefficient for water

    ! calculate saturated vapour pressure (function of temperature) air above and soil
    ! air space (kPa)
    ! note: both 0.622 and 0.611213 or saturation vapour pressure at 0c (MPa)
    if (hourtemp <= freeze) then
        ! saturation vapour pressure of air (kPa)...
        esat_air  = 0.611213d0 * exp( 22.4422d0 * (hourtemp-freeze) / ( 272.186d0 + (hourtemp-freeze) ) )
        ! saturation vapour pressure (kPa) at surface - assume saturation...
        esat_soil = 0.611213d0 * exp( 22.4422d0 * (ts-freeze) / ( 272.186d0 + (ts-freeze) ) )
    else 
        esat_air  = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * hourtemp - 4717.306081d0 ) / ( hourtemp - 35.86d0 ) )
        esat_soil = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * ts - 4717.306081d0 ) / ( ts - 35.86d0 ) )
    endif
    ! vapour pressure of air (kPa)
    ea   = esat_air - hourvpd  
    ! vapour pressure in soil airspace (kPa), dependent on soil water potential - jones p.110. vw=partial molal volume of water...
    esurf = esat_soil * exp( 1d6 * (sum(swp(1:2) * thickness(1:2)) / sum(thickness(1:2))) * vw / ( rcon * ts ) )
    ! soil conductance to water vapour diffusion (m s-1)...
    gws = por * diff / ( tort * drythick )
    ! generate final output
    qe_flux = lambda_soil * rho * 0.622d0 / ( 1d-3 * hourpress ) * ( ea - esurf ) / ( 1d0 / gaw + 1d0 / gws )

  end function qe_flux
  !
  !---------------------------------------------------------------------- 
  !
  double precision function qs_flux( ts )

    ! latent energy loss from snow surface !
    ! negative is water changing from ice to water

    use gv_hourscale,      only: freeze, gaw, hourtemp, hourvpd
    use gv_meteo,          only: air_density_kg, Rcon

    implicit none

    ! arguments..
    double precision, intent(in) :: ts

    ! local variables..
    double precision :: ea, esat, esurf, rho, tc

    gaw = soil_conductance                                 ! load soil conductance from memory (m.s-1)
    rho = air_density_kg                                   ! density of air (kg m-3) (t-dependent)

    ! saturation vapour pressure of air (kPa)...
    esat  = 0.611213d0 * exp( 22.4422d0 * (hourtemp-freeze) / ( 272.186d0 + (hourtemp-freeze) ) )
    ea   = esat - hourvpd         ! vapour pressure of air
    ! we only want to try this if things are frozen
    tc   = ts - freeze
    if (tc > 0d0) then
        qs_flux = 0d0
        return
    endif
    ! saturation vapour pressure (kPA) with respect to ice at snow surface...
    esurf = 0.611213d0 * exp( 22.4422d0 * (ts-freeze) / ( 272.186d0 + (ts-freeze) ) )
    qs_flux = gaw * ( (lhs * 0.611213d0) / (Rcon * ts) ) * (ea-esurf)

    return

  end function qs_flux
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_balance( soil_layer )

    ! integrator for soil gravitational drainage !

    use gv_soil_structure, only: field_capacity, iceprop, porosity, thickness, waterfrac, watergain, waterloss
    use log_tools
    use math_tools,        only: dxsav, kmax, ode_int
    use soil_air,          only: drainlayer, liquid, slayer, unsat

    implicit none

    ! arguments..
    integer,intent(in) :: soil_layer

    ! local variables..
    integer,parameter :: nvar = 1
    integer           :: nbad, nok
    double precision  :: change, eps, h1, hmin, newwf, x1, x2, ystart(nvar)

    ! --calculations begin below--
    eps   = 1.0d-4  ! accuracy
    h1    = 0.001d0 ! first guess at step size
    hmin  = 0d0     ! minimum step size
    kmax  = 100     ! maximum number of iterations (?)
    x1    = 1d0     ! x1 and x2 define the integrator range for process to
                    ! occur. i.e. 2-1 = 1, 1 step for integrator
    x2    = 2d0
    dxsav = ( x2 - x1 ) * 0.05d0

    ! liquid content of the soil layer, i.e. fraction avaiable for drainage
    liquid     = waterfrac( soil_layer ) * ( 1d0 - iceprop( soil_layer ) )     ! liquid fraction
    ! soil water capacity of the current layer
    drainlayer = field_capacity( soil_layer )
    
    ! unsaturated volume of layer below (m3 m-2)..
    unsat      = max( 0d0 , ( porosity( soil_layer+1 ) - waterfrac( soil_layer+1 ) ) &
                          * thickness( soil_layer+1 ) / thickness( soil_layer )     )
    ! soil layer passed in common block for integrator
    slayer     = soil_layer

    ! initial conditions; i.e. is there liquid water and more water than layer
    ! can hold
    if ( ( liquid .gt. 0d0 ) .and. ( waterfrac( soil_layer ) .gt. drainlayer ) ) then
      ! there is liquid water..
      ystart(1) = waterfrac( soil_layer )    ! total layer   
      ! call integrator
      call ode_int( 'soil_balance:soil_water_store' , ystart , nvar , x1 , &
                     x2 , eps , h1 , hmin , nok , nbad , soil_water_store )
      newwf  = ystart(1)
      ! convert from waterfraction to absolute amount (m)
      change = ( waterfrac( soil_layer ) - newwf ) * thickness( soil_layer )
      ! update soil layer below with drained liquid
      watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
      waterloss(    soil_layer  ) = waterloss(   soil_layer   ) + change
    end if

    if ( waterloss( soil_layer ) .lt. 0d0 ) then
        call write_log( 'waterloss probem in soil_balance' , msg_error , __FILE__ , __LINE__ )
    end if

  end subroutine soil_balance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_store( time_dummy , y , dydt )

    ! determines gravitational water drainage !

    use gv_scale_declarations, only: max_nos_iterations, time
    use soil_air,              only: drainlayer, liquid, soil_conductivity, unsat

    implicit none

    ! arguments..
    double precision,intent(in)  :: y(max_nos_iterations)
    double precision,intent(in)  :: time_dummy ! dummy argument, provided for ode_int
    double precision,intent(out) :: dydt(max_nos_iterations)

    ! local variables..
    double precision    :: drainage

    drainage = soil_conductivity( y(1) ) * time%seconds_per_step

    if ( y(1) .le. drainlayer ) then  ! gravitational drainage above field_capacity 
      drainage = 0d0
    end if
    if ( drainage .gt. liquid ) then  ! ice does not drain
      drainage = liquid
    end if
    if ( drainage .gt. unsat ) then   ! layer below cannot accept more water than unsat
      drainage = unsat
    end if

    dydt(1) = -drainage               ! waterloss from this layer

  end subroutine soil_water_store
  !
  !----------------------------------------------------------------------
  !
  subroutine thaw_depth(ts)

    ! determines layer ice fraction, and depth of thaw !

    use gv_hourscale,          only: freeze
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: iceprop, soil_temp, thickness

    implicit none

    ! arguments..
    double precision,intent(in) :: ts

    ! local variables..
    integer :: k, liquid(0:grid%soil), numthaw
    double precision    :: depthtotop(grid%soil), middepth(grid%soil), &
               root, split, thaw(grid%soil), topt, topthick

    thaw       = -9999d0
    numthaw    = 1          ! There may be two or more thaw points, so thaw is an array
    depthtotop = 0d0        ! Records depth to start of each layer
    iceprop    = 0d0

    ! has my surface fronzen or not
    if ( ts .gt. freeze ) then
        liquid(0) = 1
    else
        liquid(0) = 0
    end if

    ! determine which soil layers are above freezing and therefore should be
    ! liquid
    do k = 1 , grid%soil   ! Check if layer centres are frozen
       if ( soil_temp(k) .gt. freeze ) then
           liquid(k) = 1
       else
           liquid(k) = 0
       end if
    enddo

    ! therefore determine the thaw depth
    do k = 1 , grid%soil    ! Locate thaw depth
       if ( k .eq. 1 ) then  ! For layer 1, upper boundary is surface temperature
           topthick      = 0d0
           middepth(k)   = 0.5d0 * thickness(k)
           depthtotop(k) = 0d0
       else
           topthick        = thickness( k-1 )
           depthtotop( k ) = sum( depthtotop) + thickness( k-1 ) ! Increment
           middepth( k )   = depthtotop( k ) + 0.5d0 * thickness( k )
       end if

       ! detewrmine of the soil layer above is a different state from current
       if ( liquid( k-1 ) .ne. liquid( k ) ) then
           ! Change of state - locate thaw
           if ( k .eq. 1 ) then  !for layer 1, upper boundary is surface temperature
             topt = ts
           else
             topt = soil_temp( k-1 )
           end if
           ! Location of thaw relative to layer midpoint (m)
           root = ( freeze - soil_temp(k) ) * ( 0.5d0 * thickness(k) + 0.5d0 * topthick ) &
                        / ( topt -soil_temp(k) )
           thaw( numthaw ) = middepth( k ) - root  ! Determine thaw depth

           ! Now ice fractions of the thaw layer
           if ( thaw( numthaw ) .gt. depthtotop( k ) ) then 
               ! Thaw is in layer k
               ! fraction of top half of layer thawed/unthawed..
               split = ( thaw( numthaw ) - depthtotop( k ) ) / ( 0.5d0 * thickness( k ) )
               if ( soil_temp(k) .lt. freeze ) then
                   ! soil layer is frozen, but how much
                   iceprop( k ) = iceprop( k ) + 0.5d0 * ( 1d0 - split )
               else
                   ! soil layer with be mostly melted so...
                   iceprop( k ) = iceprop( k ) + 0.5d0 * split
                   ! if soil layer is not the surface and some of the layer is frozen,
                   ! therefore this is a phase change and must be frozen at the bottom
                   if ( k .gt. 1 ) iceprop( k-1 ) = iceprop( k-1 ) + 0.5d0  ! bottom half of k-1 is frozen
               end if
           else
               ! Thaw is in layer k-1
               ! fraction of bottom half of layer-1 thawed/unthawed..
               split = ( depthtotop( k ) - thaw( numthaw ) ) / ( 0.5d0 * thickness( k - 1 ) )
               if ( soil_temp( k-1 ) .lt. freeze ) then
                   iceprop( k-1 ) = iceprop( k-1 ) + 0.5d0 * ( 1d0 - split )
               else
                   iceprop( k-1 ) = iceprop( k-1 ) + 0.5d0 * split
                   iceprop( k   ) = iceprop(  k  ) + 0.5d0  ! top half of layer k is frozen
               end if
           end if
           numthaw = numthaw + 1       ! next thaw has separate storage location
       else
           ! No change of state
           if ( liquid( k-1 ) + liquid( k ) .eq. 2d0 ) then
               ! Both water..
               iceprop( k ) = iceprop( k )
               if ( k .gt. 1 ) iceprop( k-1 ) = iceprop( k-1 )
           else
               ! Both ice..
               iceprop( k ) = iceprop( k ) + 0.5d0
               if ( k .gt. 1 ) iceprop( k-1 ) = iceprop( k-1 ) + 0.5d0
           end if
       end if
    enddo

  end subroutine thaw_depth
  !
  !----------------------------------------------------------------------
  !
  double precision function snow_conductivity( i )

    ! thermal conductivity of snow W m-1 K-1  !

    use gv_snow_info,           only : Dsnow, Sice
    use gv_soil_structure,      only : k_ice

    implicit none

    ! arguments..
    integer,intent(in) :: i

    double precision :: rho

    rho = 100d0
    if ( Dsnow(i) .gt. 0d0 ) rho = Sice(i) / Dsnow(i)
    snow_conductivity = k_ice * (rho / rho_ice)**2

  end function snow_conductivity
  !
  !----------------------------------------------------------------------
  !
  double precision function thermal_conductivity( i )

    ! thermal conductivity (W m-1 K-1)  !
    ! Hillel p.295                      !

    use gv_soil_structure, only: iceprop

    implicit none

    ! arguments..
    integer,intent(in) :: i

    if ( i .lt. 4 ) then
      thermal_conductivity = 0.34d0 + ( 0.58d0 - 0.34d0 ) * iceprop(i)
    else if ( ( i .ge. 4 ) .and. ( i .lt. 10 ) ) then
      thermal_conductivity = 1.5d0  + ( 2.05d0 -  1.5d0 ) * iceprop(i)
    else
      thermal_conductivity = 1.69d0 + ( 3.63d0 - 1.69d0 ) * iceprop(i)
    end if

  end function thermal_conductivity
  !
  !----------------------------------------------------------------------
  !
  subroutine water_fluxes( time )

    ! waterfrac is m3 m-3, soilwp is MPa !

    use gv_hourscale,          only: freeze, hourts, Qe, Qm, surface_watermm, totet, unintercepted
    use gv_scale_declarations, only: grid, time_holder
    use gv_soil_structure,     only: abovebelow, drythick, fraction_uptake &
                                    ,rooted_layers, soilR, soilR1, soilR2  &
                                    ,root_mass, root_length, thickness     &
                                    ,snow, watergain, waterloss
    use log_tools
    use soil_air,              only: soil_resistance, soil_water_potential

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    integer            :: i, rr

    ! snow not currently in use but is available for wettinglayer calculation,
    ! for possible inclusion
    snow = 0d0

    ! determine wettinglayers and the dryzone
    call wetting_layers( time )

    ! From which layer is evap water withdrawn?..
    if ( drythick .lt. thickness( 1 ) ) then
        rr = 1     ! The dry zone does not extend beneath the top layer
    else
        rr = 2     ! The dry zone does extend beneath the top layer
    end if

    ! determine water loss in upper layers due to evaporation.
    ! NOTE: that Qm and Qs are not added here due to their being added to the
    ! soil layers directly in "snow"
    if ( Qe .lt. 0d0 ) then
        ! Evaporation (t m-2 t-1, m t-1)..
        waterloss( rr ) = waterloss( rr ) + 0.001d0 * ( -Qe / lambda_soil ) * time%seconds_per_step
    else
        ! Dew formation (t m-2 t-1, m t-1)..
        watergain( 1 )  = watergain( 1 )  + 0.001d0 * (  Qe / lambda_soil ) * time%seconds_per_step
    end if

    ! ensure evapotransporation cannot be less than zero, i.e. roots do not add
    ! water to the soil
    totet = max( 0d0 , totet )
    do i = 1 , rooted_layers            !water loss from each layer 
       waterloss( i ) = waterloss( i ) + totet * fraction_uptake( i ) * abovebelow
    enddo

    do i = 1 , grid%soil
       ! determines water movement between soil layers to due drainage down the
       ! profile
       call soil_balance( i )
       if ( waterloss( i ) .lt. 0d0 ) then
           write(message,*)'trouble with water loss = ',waterloss(i),' in layer ',i
           call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
       end if
    enddo

    ! how much surface water infiltrantes the first soil layer in the current
    ! step. Water which does not infiltrate in a single step is considered run
    ! off
    call infiltrate

    ! Find SWP & soil resistance without updating waterfrac yet (do that in waterthermal)
    call soil_water_potential
    call soil_resistance(soilR,soilR1,soilR2,root_length,root_mass,rooted_layers)

    ! count water which is unintercepted of output to file
    if (time%step == 1) then
        unintercepted = surface_watermm
    else
        ! determine total water reaching ground for day summation
        unintercepted = unintercepted + surface_watermm
    end if

  end subroutine water_fluxes
  !
  !----------------------------------------------------------------------
  !
  subroutine water_thermal

    ! redistribute heat according to !
    ! water movement in the soil     !

    use gv_hourscale,          only: discharge, hourtemp, hourts, Qe, underflow,freeze
    use gv_scale_declarations, only: grid, time
    use gv_soil_structure,     only: iceprop, mineralfrac, organicfrac, pptgain, soil_temp, thickness, &
                                     waterfrac, watergain, watericemm, waterloss
    use log_tools

    implicit none

    ! local variables..
    integer :: i
    double precision    :: hold_var    & ! 
              ,evap        & ! total surface evaporation in timestep (tonnes.m-2.t-1) 
              ,heat        & ! Thermal energy contained within soil layer (J.m-2)
              ,heatgain    & ! Thermal energy gained to soil layer due to water movement in (J.m-2)
              ,heatloss    & ! Thermal energy lost from layer due to water movement out (J.m-2)
              ,icecontent  & ! absolute ice content (m)
              ,newheat     & ! Net thermal energy change due to water movement (J.m-2)
              ,volhc       & ! Volumetric heat capacity (J.m-3.K-1)
              ,watercontent  ! absolute water content of soil layer (m or tonnes.m-2)

    do i = 1 , grid%soil    
       ! Determine volumetric heat capacity (J.m-3.K-1) of soil layer; based on
       ! water, ice, organic and mineral fractions of soil
       volhc = volumetric_heat_capacity(mineralfrac(i),organicfrac(i) &
                                       ,iceprop(i),waterfrac(i))
       ! Current thermal energy in soil layer (J.m-2)
       heat           = soil_temp(i)*volhc*thickness(i) 
       ! water loss in m * heat capacity of water * temp => J m-2
       heatloss       = waterloss( i ) * 4.2d6 * soil_temp( i )
       ! water gain in m * heat capacity of water * temp => J m-2
       heatgain       = watergain( i ) * 4.2d6 * soil_temp( i ) + pptgain( i ) * 4.2d6 * hourtemp
       ! Liquid water content of layer (m or tonnes.m-2)
       watercontent   = ( waterfrac(i) * ( 1d0 - iceprop(i) ) ) * thickness(i) 
       ! Ice water content of layer (m or tonnes.m-2)
       icecontent     = ( waterfrac(i) *     iceprop(i)      ) * thickness(i) 
       ! check what the difference is for the water content. To ensure mass
       ! balance in cases where sublimation occurs.
       ! special condition for when the soil profile is entirely frozen but
       ! sublimation occurs at the surface. We remove the water from the
       ! icecontent instead. This is only done for the top 2 soil layers to be
       ! consistent with the calculation of Qe
       ! ACTUALLY SHOULD BE DONE WITH DEPTH ACTUALLY REACHED BASED ON DRYTHICK
       hold_var = watercontent + watergain(i) + pptgain( i ) - waterloss( i )
       if (hold_var < 0d0 .and. icecontent > 0d0 .and. i < 3) then
           icecontent = max(0d0, icecontent-max(0d0,waterloss(i)-(waterloss(i)+hold_var)))
       endif
       ! Net change in water content (m or tonnes.m-2); max condition to ensure
       ! no negative content values
       watercontent   = max( 0d0 , watercontent + watergain(i) + pptgain( i ) - waterloss( i ) ) 
       ! Determine new total water content of layer (m or tonnes.m-2)
       waterfrac( i ) = ( watercontent + icecontent ) / thickness( i )
       ! Determine new ice proportion of water in layer; max condition added
       ! 22/03/2011 to prevent NaN values
       if ( ( watercontent + icecontent ) == 0d0 ) then
           iceprop(i) = 0d0
           call write_log("water_thermal: iceprop is zero!" , msg_warning , __FILE__ , __LINE__ )
       else
           iceprop(i) = icecontent / ( watercontent + icecontent )
       end if
 
       ! If thermal exchange has a net value then calculate new soil layer
       ! temperature
       if ( heatgain + heatloss .ne. 0d0 )then     
           ! Net thermal change in soil layer (J.m-2)
           newheat     = heat - heatloss + heatgain
           ! Determine NEW volumetric heat capacity (J.m-3.K-1) of soil layer;
           ! based on water, ice, organic and mineral fractions of soil
           volhc = volumetric_heat_capacity(mineralfrac(i),organicfrac(i) &
                                           ,iceprop(i),waterfrac(i))
           soil_temp(i) = newheat / ( volhc * thickness(i) )
       end if
       ! determine mm or kg.m-2 of water in each soil layer
       watericemm(i) = 1d3 * waterfrac(i) * thickness(i) 
    enddo
    ! determine water which has drained out the bottom of the soil column
    underflow = watergain( grid%core )
    discharge = discharge + watergain( grid%core ) * 1d3
 
  end subroutine water_thermal
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine prescribed_water_fluxes( time )
  
    ! soil water and thermal fluxes based on prescribed soil water content in
    ! meteorology forcing file. 
   
    use gv_hourscale,          only: freeze, hourts, discharge, hourtemp, Qe &
                                     ,underflow
    use gv_scale_declarations, only: grid, met, time_holder
    use gv_soil_structure,     only: iceprop, mineralfrac, organicfrac, pptgain &
                                     ,soil_temp, thickness, waterfrac &
                                     ,watergain, watericemm, waterloss
    use gv_clim,               only: swc10, swc20, swc30, swc40, swc50, swc60 &
                                     ,swc70, swc80, swc90, swc100, swc110 &
                                     ,swc120, swc130, swc140, swc150

    use log_tools
    
    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    
    ! local variables
    integer          :: i
    double precision :: hold_var  & !
              ,evap      & ! total surface evaporation in timestep (tonnes m-2 t-1)
              ,heat        & ! Thermal energy contained within soil layer (J.m-2)
              ,heatgain    & ! Thermal energy gained to soil layer due to water movement in (J.m-2)
              ,heatloss    & ! Thermal energy lost from layer due to water movement out (J.m-2)
              ,icecontent  & ! absolute ice content (m)
              ,newheat     & ! Net thermal energy change due to water movement (J.m-2)
              ,volhc       & ! Volumetric heat capacity (J.m-3.K-1)
              ,watercontent  ! absolute water content of soil layer (m or tonnes.m-2)

    ! calculate evaporation from soil layers, infiltration, drainage, soil water
    ! potential, soil resistance, runoff
    call water_fluxes( time )

    ! assign waterfrac (m3/m3) per layer from met forcing file

    do i = 1 , grid%soil   

       ! assign water content from met forcing file
       if ( i == 1) then
           waterfrac( i ) = swc10
       endif
       if (i == 2) then
           waterfrac( i ) = swc20
       endif
       if (i == 3) then
           waterfrac( i ) = swc30
       endif
       if (i == 4) then
           waterfrac( i ) = swc40
       endif
       if (i == 5) then
           waterfrac( i ) = swc50
       endif
       if (i == 6) then
           waterfrac( i ) = swc60
       endif
       if (i == 7) then
           waterfrac( i ) = swc70
       endif
       if (i == 8) then
           waterfrac( i ) = swc80
       endif
       if (i == 9) then
           waterfrac( i ) = swc90
       endif
       if (i == 10) then
           waterfrac( i ) = swc100
       endif
       if (i == 11) then
           waterfrac( i ) = swc110
       endif
       if (i == 12) then
           waterfrac( i ) = swc120
       endif
       if (i == 13) then
           waterfrac( i ) = swc130
       endif
       if (i == 14) then
           waterfrac( i ) = swc140
       endif
       if (i == 15) then
           waterfrac( i ) = swc150
       endif

       ! waterloss(i)
       ! watergain(i)
       call soil_balance( i )
       if ( waterloss( i ) .lt. 0d0 ) then
           write(message,*)'trouble with water loss = ',waterloss(i),' in layer',i
           call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
       end if 


       ! Determine volumetric heat capacity (J.m-3.K-1) of soil layer; based on
       ! water, ice, organic and mineral fractions of soil
       volhc = volumetric_heat_capacity(mineralfrac(i),organicfrac(i) &
                                       ,iceprop(i),waterfrac(i))
       ! Current thermal energy in soil layer (J.m-2)
       heat           = soil_temp(i)*volhc*thickness(i) 
       ! water loss in m * heat capacity of water * temp => J m-2
       heatloss       = waterloss( i ) * 4.2d6 * soil_temp( i )
       ! water gain in m * heat capacity of water * temp => J m-2
       heatgain       = watergain( i ) * 4.2d6 * soil_temp( i ) + pptgain( i ) * 4.2d6 * hourtemp
       ! Liquid water content of layer (m or tonnes.m-2)
       watercontent   = ( waterfrac(i) * ( 1d0 - iceprop(i) ) ) * thickness(i) 
       ! Ice water content of layer (m or tonnes.m-2)
       icecontent     = ( waterfrac(i) *     iceprop(i)      ) * thickness(i) 
       ! check what the difference is for the water content. To ensure mass
       ! balance in cases where sublimation occurs.
       ! special condition for when the soil profile is entirely frozen but
       ! sublimation occurs at the surface. We remove the water from the
       ! icecontent instead. This is only done for the top 2 soil layers to be
       ! consistent with the calculation of Qe
       ! ACTUALLY SHOULD BE DONE WITH DEPTH ACTUALLY REACHED BASED ON DRYTHICK
       hold_var = watercontent + watergain(i) + pptgain( i ) - waterloss( i )
       if (hold_var < 0d0 .and. icecontent > 0d0 .and. i < 3) then
           icecontent = max(0d0, icecontent-max(0d0,waterloss(i)-(waterloss(i)+hold_var)))
       endif
       ! Net change in water content (m or tonnes.m-2); max condition to ensure
       ! no negative content values
       !watercontent   = max( 0d0 , watercontent + watergain(i) + pptgain( i ) - waterloss( i ) ) 
       ! Determine new total water content of layer (m or tonnes.m-2)
       !waterfrac( i ) = ( watercontent + icecontent ) / thickness( i )

       ! Determine new ice proportion of water in layer; max condition added
       ! 22/03/2011 to prevent NaN values
       if ( ( watercontent + icecontent ) == 0d0 ) then
           iceprop(i) = 0d0
           call write_log("water_thermal: iceprop is zero!" , msg_warning , __FILE__ , __LINE__ )
       else
           iceprop(i) = icecontent / ( watercontent + icecontent )
       end if
 
       ! If thermal exchange has a net value then calculate new soil layer
       ! temperature
       !if ( heatgain + heatloss .ne. 0d0 )then     
           ! Net thermal change in soil layer (J.m-2)
       !    newheat     = heat - heatloss + heatgain
           ! Determine NEW volumetric heat capacity (J.m-3.K-1) of soil layer;
           ! based on water, ice, organic and mineral fractions of soil
       !    volhc = volumetric_heat_capacity(mineralfrac(i),organicfrac(i) &
       !                                    ,iceprop(i),waterfrac(i))
       !    soil_temp(i) = newheat / ( volhc * thickness(i) )
       !end if
       ! determine mm or kg.m-2 of water in each soil layer
       watericemm(i) = 1d3 * waterfrac(i) * thickness(i) 
    enddo
    ! determine water which has drained out the bottom of the soil column
    underflow = watergain( grid%core )
    discharge = discharge + watergain( grid%core ) * 1d3

  end subroutine
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  double precision function wetted_surface_evap()

    ! evaporation from wetted surfaces (mm t-1) !

    use gv_hourscale,          only: freeze, hourrnet, hourtemp, hourwind, hourvpd
    use gv_scale_declarations, only: time
    use gv_veg,                only: canopy_height
    use gv_clim,               only: cp_air
    use log_tools

    implicit none

    ! local variables..
    double precision :: d, htemp, psych, ra, rho, s, slope, z, z0, lambda

    ! boundary layer constants..
    d  = 0.75d0 * canopy_height
    z0 = 0.1d0 * canopy_height
    z  = canopy_height + 2d0
    htemp  = hourtemp - freeze ! convert Kelvin to Celcius

    rho    = 353.0d0 / hourtemp                          ! density of air (g m-3) (t-dependent)
    ra     = ( 1d0 / vonkarman**2 * hourwind ) * ( log( z - d ) / z0 )**2 ! aerodynamic resistance
    s      = 6.1078d0 * 17.269d0 * 237.3d0 * exp( 17.269d0 * htemp / ( 237.3d0 + htemp ) )
    slope  = 100d0 * ( s / ( 237.3d0 + htemp )**2 )       ! slope of saturation vapour pressure curve (t-dependent)
    ! latent heat of vapouriation / sublimation (J.kg-1)
    if (htemp > 0d0) then
       lambda=1000d0*(2501.0d0-2.364d0*htemp)
    else
       ! sublimation must be occuring
       lambda = 2.83d6
    endif
    psych  = 100d0 * ( 0.646d0 * exp( 0.00097d0 * htemp ) ) ! psychrometer constant
    if ( slope * hourrnet .gt. rho * cp_air * hourvpd / ra ) then
        ! Penman-Monteith equation (kg m-2 s-1)
        wetted_surface_evap = ( slope * hourrnet + rho * cp_air * hourvpd / ra ) &
                                 / ( lambda * ( slope + psych ) )
    else
        ! Actually dewfall occurs
        wetted_surface_evap = 0d0
    end if

    ! convert from kg m-2 s-1 to mm t-1..
    wetted_surface_evap = wetted_surface_evap * time%seconds_per_step

    if ( wetted_surface_evap .lt. -1000d0 ) then
      call write_log( 'Problem in wetted_surface_evap' , msg_warning , __FILE__ , __LINE__ )
    end if

  end function wetted_surface_evap
  !
  !----------------------------------------------------------------------
  !
  subroutine wetting_layers( time )

    ! surface wetting and drying determines !
    ! thickness of dry layer and thus Qe    !

    use gv_hourscale,          only: freeze, hourtemp, Qe, surface_watermm
    use gv_scale_declarations, only: time_holder, grid
    use gv_soil_structure,     only: drythick, porosity, snow, thickness, wettingbot, wettingtop, &
                                     waterfrac
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
!    integer :: ar1(1), ar2, ar1b, 
    integer :: soil_layer, top_20cm
    double precision    :: dmin!airspace, diff, netc

    ! minimum thickness (m) of dry layer in soil
    dmin = 0.01d0

    ! something simplier
    ! first calculate the number of soil layers to 20 cm depth
    do soil_layer = 1, grid%soil
      if (sum(thickness(1:soil_layer)) >= 0.2d0) then
          top_20cm = soil_layer
          exit
      endif
    enddo
    ! now based on ratio of soil water content over its field capacity (or porosity) determine
    ! the drythick
    drythick = 0.2d0 * min(1d0,1d0 - ( sum(waterfrac(1:top_20cm)) / sum(porosity(1:top_20cm))))
    drythick = max(dmin,drythick)
    if (drythick >= 0.2d0) then
       ! There is no water left in top 20 cm for soil evaportation 
       write(message,*) "There is no water left in top 20 cm for direct soil evap!"
       call write_log( message , msg_warning , __FILE__ , __LINE__ )
       print*,"waterfrac",waterfrac
    endif 
!    ! TLS: MAT'S original code commented on 3rd July 2016
!    airspace = porosity( 1 )
!    ! Soil LE should be withdrawn from the wetting layer with the smallest depth..
!    ar1  = minloc( wettingbot , MASK = wettingbot .gt. 0. )
!    ar1b = sum( ar1 )   ! convert ar1 to scalar, and make sure always greater than or equal to 1
!    if ( ar1b .eq. 0. ) then
!      ! There is no water left in soils!!
!      write(message,*) "There is no water left in any soil layer!"
!      call write_log( message , msg_warning , __FILE__ , __LINE__ )
!    end if
!
!    ! Calulate the net change in wetting in the top zone..
!    netc = ( 0.001 * Qe / lambda_soil * time%seconds_per_step ) / airspace + &
!              ( surface_watermm * 0.001 + snow ) / airspace   ! m
!    if ( netc .gt. 0. ) then      ! Wetting
!      ! resaturate the layer if top is dry and recharge is greater than drythick
!      if ( ( netc .gt. wettingtop( ar1b ) ) .and. ( wettingtop( ar1b ) .gt. 0. ) ) then
!        diff = netc - wettingtop( ar1b )  ! extra water to deepen wetting layer
!        wettingtop( ar1b ) = 0.
!        if ( ar1b .gt. 1 ) then
!          ! Not in primary layer (primary layer can't extend deeper)..
!          wettingbot( ar1b ) = wettingbot( ar1b ) + diff
!        end if
!        drythick = dmin
!      else
!        if ( wettingtop( ar1b ) .eq. 0. ) then
!          ! surface is already wet, so extend depth of this wet zone..
!          if ( ar1b .gt. 1 ) then
!            ! not in primary layer (primary layer can't extend deeper)..
!            wettingbot( ar1b ) = wettingbot( ar1b ) + netc          
!            if ( wettingbot(ar1b) .ge. wettingtop(ar1b-1) ) then
!              ! Layers are conterminous..
!              wettingtop( ar1b - 1 ) = wettingtop(ar1b)
!              wettingtop(   ar1b   ) = 0.     ! remove layer
!              wettingbot(   ar1b   ) = 0.
!            end if
!          end if
!        else
!          ! Create a new wetting zone..
!          wettingtop( ar1b + 1 ) = 0.
!          wettingbot( ar1b + 1 ) = netc
!        end if
!        drythick = dmin
!      end if
!    else    ! drying
!      wettingtop( ar1b ) = wettingtop( ar1b ) - netc         ! Drying increases the wettingtop depth
!      if ( wettingtop( ar1b ) .ge. wettingbot( ar1b ) ) then ! Wetting layer is dried out.
!        diff = wettingtop( ar1b ) - wettingbot( ar1b )       ! How much more drying is there?                 
!        wettingtop( ar1b ) = 0.
!        wettingbot( ar1b ) = 0.
!        ar2 = ar1b - 1
!        if ( ar2 .gt. 0 ) then  ! Move to deeper wetting layer
!          wettingtop( ar2 ) = wettingtop( ar2 ) + diff    ! dry out deeper layer
!          drythick = max( dmin , wettingtop( ar2 ) )
!        else    ! no deeper layer
!          drythick = thickness( 1 )   ! layer 1 is dry
!        end if
!      else
!        drythick = max( dmin , wettingtop( ar1b ) )
!      end if
!    end if
!    if ( drythick .eq. 0. ) then
!      call write_log( 'Problem in drythick' , msg_warning , __FILE__ , __LINE__ )
!    end if

  end subroutine wetting_layers
  !
  !----------------------------------------------------------------------
  !  
  double precision function volumetric_heat_capacity(mineral_fraction,organic_fraction &
                                        ,ice_proportion,water_fraction)

    ! calculates soil layer specific volumetric heat capacity based on mineral,
    ! organic ice and water fractions (Hillel 1980 p. 294)

    implicit none

    ! declare inputs
    double precision, intent(in) :: mineral_fraction, &
                        organic_fraction, &
                          water_fraction, &
                          ice_proportion 

    ! volumetric heat capacity of soil layer (J m-3 K-1)
    volumetric_heat_capacity = 2d6 * mineral_fraction + 2.5d6 * organic_fraction  &
                             + 4.2d6 * ( water_fraction * ( 1d0 - ice_proportion ) ) &
                             + 1.9d6 * water_fraction * ice_proportion
!print*,"vhc",volumetric_heat_capacity 
  end function volumetric_heat_capacity
  !
  !----------------------------------------------------------------------
  !
end module soil_functions
!
!------------------------------------------------------------------------
!
