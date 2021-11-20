! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module main_spa

implicit none

private

public :: spa

contains

!
!---------------------------------------------------------------------------------------------------
!
  !
  !-------------------------------------------------------------------------------------------------
  !
  subroutine spa(met,pheno,pars,deltat,nosteps,out_steps,nodays,lat,long,year_in &
                ,lai,GPP,FLUXES,POOLS,nopheno,nopars,nomet,nopools,nofluxes &
                ,soilwater,evapotrans,sensible,daily_weighted_SWP           &
                ,daily_weighted_soilR,soilevap,wetcanopyevap,potevap,netrad &
                ,porosity_out,sun_canopy_temperature_out,shade_canopy_temperature_out&
                ,canopy_temperature_out,ground_heat_out,soil_conductance_out)

  !! Soil Plant Atmosphere model coupled to DALEC !!

  use soil_air,              only: roots
  use gv_soil_structure,     only: rooted_layers, root_mass, root_length
  use gv_veg,                only: stock_roots
  use gv_hourscale,          only: canopy_store
  use canopy,                only: timestep_calcs
  use gv_scale_declarations, only: time, user_opts
  use gv_soil_structure,     only: waterfrac,field_capacity,porosity 
  use gv_daily_averages,     only: daily_sun_canopy_temperature, &
                                   daily_shade_canopy_temperature, &
                                   daily_canopy_temperature, &
                                   daily_ground_heat_flux, &
                                   daily_soil_conductance
  use log_tools
  use spa_io,                only: handle_output, start_spa, & 
                                   update_met_drivers, update_phenology_drivers

  implicit none

  ! arguments
  ! declare input variables
  integer, intent(in) :: nopars    & ! number of paremeters in vector
                        ,year_in   & ! start year of analysis
                        ,out_steps & 
                        ,nomet     & ! number of meteorological fields
                        ,nopheno   & ! number of phenology fields
                        ,nofluxes  & ! number of model fluxes
                        ,nopools   & ! number of model pools
                        ,nodays    & ! number of days in simulation
                        ,nosteps    ! number of steps in simulation

  double precision, intent(in) :: met(nomet,nosteps) & ! met drivers
                       ,pheno(nodays,nopheno) & ! phenology drivers
                       ,deltat            & ! time step in decimal days
                       ,pars(nopars)      & ! number of parameters
                       ,long              & ! site longitude
                       ,lat                 ! site latitude (degrees)

  double precision, dimension(out_steps), intent(inout) :: lai & ! leaf area index
                                             ,daily_weighted_SWP   & ! MPa
                                             ,daily_weighted_soilR & !
                                             ,sun_canopy_temperature_out &
                                             ,shade_canopy_temperature_out &
                                             ,canopy_temperature_out &
                                             ,ground_heat_out &
                                             ,soil_conductance_out &
                                             ,sensible   & ! sensible heat
                                             ,evapotrans & ! evapotranspiration (mm.day-1 .or. W.m-2)
                                             ,soilevap   & ! soil evaporation (mm.day-1 .or. W.m-2)
                                             ,wetcanopyevap & ! wet canopy evaporation (mm.day-1 .or. W.m-2)
                                             ,potevap & ! potential evaporation (mm.day-1 or W.m-2)
                                             ,netrad & ! net radiation (W.m-2)
                                             ,porosity_out &
                                             ,soilwater  & ! soil water content (mm)
                                             ,GPP  ! Gross primary productivity

  double precision, dimension(nosteps,nopools), intent(inout) :: POOLS ! vector of ecosystem pools
  double precision, dimension(out_steps,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

  ! local variables
  integer            :: yr, day, step, counter_daily, counter_hourly, start_day
  double precision               :: iwater

  ! read user config, open files & initialise (spa_io.f90)
  call initialise_spa_coupling(nosteps,nodays,nomet,met,pheno,nopheno,nopars,pars,lat,long,year_in)

  ! increment day component of model timing information
  call spa_couple_time_update(time,int(met(1,1)))

  ! initial conditions
  counter_daily = 0 ; counter_hourly = 0

  ! starting point and adjust the time%day back one day. This is cause time% day
  ! will be incremented to the correct day again in "increment_day" 
  start_day = time%day ; time%day = time%day - 1 
  ! for each year..
  do yr = 1 , time%nos_of_years
!     if (yr /= 1)call increment_time( 'year' , time )
     ! in first year don't do this to avoid resetting user inputted start day
     if (yr == 1) then
        call increment_time( 'first_year' , time , start_day)
     else
        call increment_time( 'year' , time , start_day)
     end if

     ! for each day...
     do day = start_day, time%days_in_year(yr)

        ! for runs where we want to remove moisture stress update the soil to
        ! field capacity each day
!        waterfrac = porosity !(porosity+field_capacity)*0.5
!        waterfrac = field_capacity
!        canopy_store = 0d0

        ! update output counter
        counter_daily = counter_daily + 1

        ! increment day component of model timing information
        call increment_time( 'day' , time, start_day)

        ! Reset some varables for the new day
        daily_sun_canopy_temperature = 0d0
        daily_shade_canopy_temperature = 0d0
        daily_canopy_temperature = 0d0
        daily_ground_heat_flux = 0d0
        daily_soil_conductance = 0d0

        if (.not. user_opts%use_dalec) then
           ! update phenology for lai and root C
           call update_phenology_drivers
        end if ! use_dalec

        ! Update root structure..
        call roots(stock_roots,rooted_layers,root_length,root_mass)

        ! for each sub-daily slice...
        do step = 1 , time%steps_per_day
           ! for runs where we want to remove moisture stress update the soil to
           ! field capacity each step
!           waterfrac = field_capacity !* 0.75

           ! update output counter
           counter_hourly = counter_hourly + 1

           call increment_time( 'step' , time , start_day)

           ! get met driver for step (spa_io.f90)
           call update_met_drivers( time )

           ! transfer carbon and water (canopy.f90)
           call timestep_calcs( time )

           ! write output if needed (spa_io.f90)
           call handle_output( 1 , time )

           if (out_steps == nosteps) then
               ! handle output through interface files
               call spa_couple_output_hourly(counter_hourly,step,GPP,lai,FLUXES,POOLS,nopools &
                                            ,nofluxes,nosteps,soilwater,evapotrans,sensible   &
                                            ,daily_weighted_SWP,daily_weighted_soilR,soilevap &
                                            ,wetcanopyevap,potevap,netrad,porosity_out        &
                                            ,sun_canopy_temperature_out,shade_canopy_temperature_out &
                                            ,canopy_temperature_out,ground_heat_out,soil_conductance_out)
           endif

        enddo

        if (out_steps == nodays) then
            ! handle output through interface files
            call spa_couple_output_daily(counter_daily,step,GPP,lai,FLUXES,POOLS,nopools &
                                        ,nofluxes,nosteps,soilwater,evapotrans,sensible  &
                                        ,daily_weighted_SWP,daily_weighted_soilR,soilevap&
                                        ,wetcanopyevap,potevap,netrad,porosity_out       &
                                        ,sun_canopy_temperature_out,shade_canopy_temperature_out&
                                        ,canopy_temperature_out,ground_heat_out,soil_conductance_out)
        endif

     enddo

#ifdef USE_NETCDF
     ! check if it is time to write a restart file...
     call check_restart
#endif

  enddo

!$$ SHOULD THIS BE HERE, OR SOMEWHERE ELSE?
  ! write the crops output..
  if (user_opts%plant_func_type .eq. 2 .or. user_opts%plant_func_type == 3)  call handle_output( 8 , time )
!$$

#ifdef USE_NETCDF
  ! Check that, if restarts were wanted, we wrote at least one!..
  call check_restart( end_of_run = .true. )
#endif

  end subroutine spa
  !
  !-------------------------------------------
  !
  subroutine increment_time( type , time , start_day)

    ! Increment parts of the time-holder variable. !
  
    use gv_scale_declarations, only: time_holder
    use log_tools,             only: write_log, message, msg_info, msg_warning

    implicit none

    ! arguments..
    character(*),intent(in)         :: type
    type(time_holder),intent(inout) :: time
    integer, intent(inout) :: start_day

    select case (type)
    case ('first_year')
      time%year = time%year + 1
      time%step = 0
!      write(message,*)'year is: ',time%year
!      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('year')
      ! next year, reset day & step..
      time%year = time%year + 1
      time%day  = 0
      time%step = 0
!      write(message,*)'year is: ',time%year
!      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('day')
      ! next day, reset step...
      time%day  = time%day + 1
      time%run_day = time%run_day + 1
      time%step = 0
!      write(message,*)'day is: ',time%day
!      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('step')
      ! next step.
      ! Also update the count of total-steps, and the daytime..
      time%step = time%step + 1
!      write(message,*)'step-of-day is: ',time%step
!      call write_log( trim(message), msg_info  , __FILE__ , __LINE__ )
      time%steps_count = time%steps_count + 1
!      write(message,*)'number of steps so far is: ',time%steps_count
!      call write_log( trim(message), msg_info , __FILE__ , __LINE__ )

      ! Re-calculate the current time in units of days..
      ! Consider step 1 to be at 00:00 on the day, and the
      ! last step of the day to be at or before 23.59...
      time%daytime = ( time%year -1 ) * time%days_in_year(time%year)    &
                      + time%day                                  &
                       + ( dble(time%step) - 1 ) / dble(time%steps_per_day)
!      write(message,*)'daytime is: ',time%daytime
!      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case default
!      write(message,*)'Do not recognise type of update: ',type
!      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    endselect

  end subroutine increment_time
  !
  !-------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine check_restart( end_of_run )

    ! Check if it is time to write a restart file !

    use gv_scale_declarations, only: time, user_opts
    use log_tools
    use spa_restart

    implicit none

    ! arguments..
    logical, intent(in), optional :: end_of_run

    ! local variables..
    logical :: final_call

    ! If user didn't want dumps, then exit the s/r now..
    if ( .not. user_opts%make_restart_file ) return

    if (present(end_of_run)) then
      final_call = end_of_run
    else
      final_call = .false.
    end if

    if ( final_call ) then
      ! check that we wrote at least one restart dump.
      if ( prev_restart_write .eq. 0. ) then
        ! Not written yet. Do write now..
        call write_restart_file( time , user_opts%plant_func_type )
        ! And tell user..
        write(message,*) "The run completed before any restart files" &
              //" were written. A restart file has been written now."
        call write_log( message , msg_warning , __FILE__ , __LINE__ )
      end if      
    else
      ! Check if it is time to write a restart file yet..
      if ( time%year .eq. prev_restart_write + user_opts%restart_write_frequency ) then
         call write_restart_file( time , user_opts%plant_func_type )
         prev_restart_write = time%year
      end if
    end if

  end subroutine check_restart
#endif
  !
  !-------------------------------------------
  !
  subroutine initialise_spa_coupling(nos_steps,nos_days,nos_met,met_in,pheno_in &
                                    ,nos_pheno,nopars_in,pars_in,lat_in,long_in,year_in)

    use spa_config
    use log_tools, only: set_msg_level
    use spa_initialise
    use gv_clim
    use gv_hourscale
    use gv_Hydrol
    use gv_scale_declarations
    use gv_soil_structure
    use gv_Snow_Info
    use gv_veg
    use gv_meteo
    use soil_air, only: saxton_parameters, water_retention, soil_porosity
    use carbon_model_mod
    use carbon_model_crop_mod
    use gv_carbon_model
    use leaf
    use gv_daily_averages
    use spa_io, only: met_slice, pheno_slice, calculate_simulation_days

    implicit none

    ! declare input variables
    integer, intent(in) :: nos_steps,nos_met,nopars_in,nos_days,nos_pheno,year_in
    double precision :: dummy,dl
    double precision :: lat_in,long_in
    double precision, dimension(nos_met,nos_steps), intent(in) :: met_in
    double precision, dimension(nos_days,nos_pheno), intent(in) :: pheno_in
    double precision, dimension(nopars_in), intent(in) :: pars_in

    ! declare local variables
    integer :: i, month, day

    ! Get the General options...
    user_opts%soil_hydrology_opts = 1 ! 1 = Saxton; 2 = Oregon soils (Sandy)
    user_opts%canopy_MDecay       = 2
    user_opts%use_dalec = .false. 
    user_opts%dead_lai_energy_balance = .false. 
    user_opts%iWUE = 1  !# see Bonan et al., 2014 for details
    user_opts%iterative_canopy = .true. ! iterate canopy to close energy balance
    user_opts%solve_canopy_temperature = 2 ! analytic(1) or numerical solution(2)
    user_opts%load_restart = .false.
    user_opts%plant_func_type = 1 ! evergreen (0), decid (1), cereal (2), tuber (3)
    grid%canopy = 10 ; grid%soil = 20 ; grid%core = grid%soil + 1
    grid%latitude = lat_in ; grid%longitude = long_in
    ! convert latitude to module variable in radians
    latitude_radians = grid%latitude * (pi / 180d0)
    time%nos_of_years = nint(nos_days/365.25) !  1
    time%steps_per_day = nos_steps / nos_days ! 24

    ! allocate default dimensions to carbon model 
    if (allocated(pars)) deallocate (pars)
    allocate(pars(nopars_default))

    ! work out what month we are in from the day of year
    if (met_in(1,1) < days_in_month(1)) then
       ! we must be in january
       month = 1 ; day = int(met_in(1,1))
    else
       ! we are not in january
       do i = 2, 12
          if ((sum(dble(days_in_month(1:(i-1))))) < met_in(1,1) .and. sum(dble(days_in_month(1:i))) >= met_in(1,1)) then
             ! ok so we have now found the month we need
             month = i ; day = met_in(1,1) - sum(dble(days_in_month(1:(i-1))))
          endif
       enddo
    endif ! month

    ! combine these into the datestime object
    write(time%start_date,"(I4,A1,I2,A1,I2,A9)")year_in,"-",month,"-",day,"_00:00:00"

    user_opts%include_leap_years = .true.
    if (allocated(time%days_in_year)) deallocate(time%days_in_year)
    allocate(time%days_in_year(time%nos_of_years))
    if (.not.user_opts%include_leap_years) then
       ! also used to control how long the simulation is if we are doing less than one year
       time%days_in_year = nos_days ! int(met_in(1,1)+real(nos_days-1d0))
    endif

    user_opts%loop_met_driver = .false.
    user_opts%loop_pheno_driver = .false.
    user_opts%print_msg_at_severity = 1
    call set_msg_level( user_opts%print_msg_at_severity )

    ! Get the initialisation input files...
    ! Standard initialisation.
    user_opts%restart_filename="hardcoded"
    user_opts%soils_filename = "hardcoded" ; user_opts%veg_filename = "hardcoded"
    user_opts%carbon_filename = "hardcoded"
    user_opts%phenology_filename = "hardcoded" ; user_opts%crops_filename="hardcoded"

    ! Meteorology driver
    user_opts%met_filename = "via_interface"
    user_opts%temp_in_kelvin = .false. ; user_opts%met_file_has_lw_rad = .false. 
    user_opts%met_file_is_nc = .false. ; user_opts%use_co2_from_met_file = .true.
    user_opts%use_ppfd_from_met_file = .true. ; user_opts%met_file_has_sfc_press = .true.
    user_opts%precip_is_rate = .false. ; user_opts%vpd_or_specific_humidity = "vpd"

   ! silence all SPA outputs
    user_opts%std_csv_output       = .false. 
    user_opts%Ecofluxes_csv_output = .true.  ! Set if you want to do calcs for deciduous veg
    user_opts%canopy_csv_output    = .true.  ! Set if you want to do calcs for deciduous veg
    user_opts%Cstock_csv_output    = .true.  ! Set if you want to do calcs for deciduous veg
    user_opts%soils_csv_output     = .true.  ! Set if you want to do calcs for deciduous veg
    user_opts%make_restart_file    = .false.
    user_opts%restart_write_frequency = 1

    ! set up arrays
    call allocate_arrays

    ! set some variables as zero at initialisation
    met_slice = 0   ; pheno_slice = 0 ; totet = 0d0 
    leaf_temp = 0d0 ; lai = 0d0  ; dead_lai = 0d0 ; lafrac = 0d0 ; lafrac_dead = 0d0
    stock_surface_litter = 0d0   ; drythick = 0d0 
    nlink = 1   ; max_raso = 0d0 ; raso = 0d0
    gppt = 0d0  ; respt = 0d0    ; transt = 0d0 ; wetev = 0d0 
    LWPstore = 0d0 ; weighted_SWP = 0d0 ; rad_pass = 0 ; temp_top = 0d0 
    psil = 0d0 ; psil_pd = 0d0 ; dew = 0d0  
    SWP = 0d0 ; gs2 = 0d0 ; rad = 0d0
    ! root resistivity MPa s g mmol −1 H2O
    rootresist = 25.0d0  ! Bonan et al 2014

    ! now load up met drivers
    if (allocated(met%ambient_co2)) then
        deallocate(met%ambient_co2 , &
                         met%lw_rad, &
                           met%ppfd, &
                         met%precip, &
                         met%sw_rad, &
                   met%sfc_pressure, &
                           met%temp, &
                            met%vpd, &
                       met%wind_spd, &
             met%clearance_fraction, &
           met%clearance_management, &
                  met%fire_fraction )
    endif
    if (.not.allocated(met%ambient_co2)) then
       allocate( met%ambient_co2( nos_steps ) , &
                      met%lw_rad( nos_steps ) , &
                        met%ppfd( nos_steps ) , &
                      met%precip( nos_steps ) , &
                      met%sw_rad( nos_steps ) , &
                met%sfc_pressure( nos_steps ) , &
                        met%temp( nos_steps ) , &
                         met%vpd( nos_steps ) , &
                    met%wind_spd( nos_steps ) , &
          met%clearance_fraction( nos_steps ) , &
        met%clearance_management( nos_steps ) , &
               met%fire_fraction( nos_steps ) )
    endif
    ! pair incoming variables with SPA local
    met%ambient_co2=dble(met_in(3,1:nos_steps)) ; met%ppfd=dble(met_in(7,1:nos_steps)) ; met%precip=dble(met_in(8,1:nos_steps))
    met%sw_rad=dble(met_in(5,1:nos_steps)) ; met%sfc_pressure=dble(met_in(9,1:nos_steps)) ; met%temp=dble(met_in(2,1:nos_steps))
    met%vpd=dble(met_in(6,1:nos_steps)) ; met%wind_spd=dble(met_in(4,1:nos_steps))

    ! Is precip in volume-per-timestep (mm), or rate (mm/sec)?
    if ( user_opts%precip_is_rate ) then
      met%precip = met%precip * time%seconds_per_step
    end if
    ! Adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2d0 ) met%wind_spd = 0.2d0

    ! calculate which years of meteorology are leap years
    call calculate_simulation_days( met )

    ! year cannot be less than 365 days if we want more than 1 year of
    ! simulation
    if (time%nos_of_years > 1 .and. time%days_in_year(1) < 365) then
        print*,"cannot have less than 365 days days_per_year if you want more than 1 year of simulation"
        stop
    endif

    ! now load up met drivers
    if (.not.allocated(pheno%lai)) then
       allocate( pheno%lai( int(nos_days) ) , &
                 pheno%rootC( int(nos_days) )   )
    endif
    ! pair incoming variables with SPA local
    pheno%lai = dble(pheno_in(1:nos_days,1)) ; pheno%rootC = dble(pheno_in(1:nos_days,2))

    ! read in the soil parameters
    thickness = 0.1d0 
    do i = 1, grid%core
       layer_depth(i) = sum(thickness(1:i))
    enddo
    soil_temp = dble(pars_in(2))
    iceprop = 0d0 ; if (dble(pars_in(2)) < 273.15) iceprop = 1d0
    soil_frac_sand = dble(pars_in(3:23))  !40.0
    soil_frac_clay = dble(pars_in(24:44)) !15.0
    where (soil_frac_sand == 0) 
       soil_frac_sand = soil_frac_sand(3)
    end where
    where (soil_frac_clay == 0) 
       soil_frac_clay = soil_frac_clay(3)
    end where
    resp_rate_temp_coeff = 0.0693d0
    snowweight = 0d0 ; snowheight = 0d0

    ! initialize multi-layer snow model
    Dsfix(:) = thickness(1)
    Dsnow(:) = 0d0
    Sice(:) = 0d0
    Sliq(:) = 0d0
    Tsnow(:) = 270d0 ! Temperature (K)
    snowalb_nir = 0.73d0 ! Near infer-red reflectance
    snowalb_par = 0.95d0 ! PAR reflectance
    Nsnow = 0
    if (snowheight > 0) then
      dl = snowheight
      Dsnow(1) = dl
      i = 1
      if (Dsnow(1) > Dsfix(1)) then
        do i = 1, grid%snow
          Dsnow(i) = Dsfix(i)
          dl = dl - Dsfix(i)
          if (dl <= Dsfix(i) .or. i == grid%snow) then
            Dsnow(i) = Dsnow(i) + dl
            exit
          end if
        end do
      end if
      Nsnow = i
      do i = 1, Nsnow
        Sice(i) = snowweight * Dsnow(i) / snowheight
      end do
    end if

    ! we will assume that each day the soil starts at field capacity, so we will
    ! work out what that is
    call saxton_parameters
    call water_retention
    call soil_porosity

    ! sanity check - for when saxton used out of range
    where (porosity <= field_capacity) 
        porosity = field_capacity + 0.05d0
    !    print*,"WARNING: Saxton equations out of bounds: porosity < field capacity"
    end where

    ! now assume initial water conditions are same as field capacity
    waterfrac = field_capacity
    organicfrac = dble(pars_in(45:65)) !0.1
    mineralfrac = 1d0-porosity-organicfrac

    ! read the veg parameters
    lafrac = 0d0 ; lafrac(1:4) = 0.25d0 
    lafrac_dead = 0d0
    nfrac = 0d0
    nfrac(1) = 0.4d0 ; nfrac(2) = 0.25d0
    nfrac(3) = 0.2d0 ; nfrac(4) = 0.15d0
    layer_height = 0d0 ; layer_height(1) = 9.0d0
    layer_height(2) = 7.0d0 ; layer_height(3) = 5.0d0
    layer_height(4) = 2.0d0 ; layer_height(5) = 1.0d0
    layer_height(6) = 0.9d0 ; layer_height(7) = 0.8d0
    layer_height(8) = 0.7d0 ; layer_height(9) = 0.6d0
    layer_height(10) = 0.5d0
    c3 = .true. ; pars(17) = 80d0 !50.0
    avN = dble(pars_in(1)) ; gplant = 5d0 ; minlwp = -2d0 ; iWUE = 0.007d0 ; WUE = 900d0
    capac = 2500d0 ; dimen = 0.02d0  ; tower_height = 11d0
    conductivity = 0 ; kappac = 33.6d0 ; kappaj = kappac*1.6d0
    max_depth = 2d0 ; root_k = 100d0; 
    cond_slope=9d0 ; cond_slope_b =1d0 ; cond_slope_c = 0.1d0; 
    swp_params=0d0; swp_params(1)= -4d0; swp_params(2)=1.6d0 ;
    swp_params(3)=0.096d0 ; swp_params(4) = 0.0184d0
    through_fall = 0.7d0 ; max_storage = 1d0
    PHUem = 0d0    ! growing degree sum at emergence
    DR_pre = 0d0   ! development rate before flowering
    DR_post = 0d0  ! development rate after flowering
    tmin = 0d0     ! minimum temperature for development
    topt = 0d0     ! optimal temperature for development
    tmax = 0d0     ! maximum temperature for development
    tmin_v = 0d0   ! cardinal temperature for vernalization: minimum temperature
    topt_v = 0d0   ! cardinal temperature for vernalization: optimum temperature
    tmax_v = 0d0   ! cardinal temperature for vernalization: maximum temperature
    VDh = 0d0      ! effective vernalization days when plants are 50% vernalized
    PHCR = 0d0     ! critical photoperiod below which no development occurs
    PHSC = 0d0

    ! call some default initialisation subroutines
    call initialise_soils
    call initialise_veg    

    ! allocate final variables
    if (.not.allocated(daily_weighted_SWP)) then
        allocate(daily_weighted_SWP(time%steps_per_day),daily_weighted_soilR(time%steps_per_day), &
                 daily_shade_canopy_temperature(time%steps_per_day),&
                 daily_sun_canopy_temperature(time%steps_per_day),&
                 daily_canopy_temperature(time%steps_per_day),&
                 daily_ground_heat_flux(time%steps_per_day),&
                 daily_soil_conductance(time%steps_per_day))
    endif

  end subroutine initialise_spa_coupling
  !
  !-------------------------------------------
  !
  subroutine spa_couple_output_daily(counter,step,GPP_out,lai_out,FLUXES,POOLS,nopools,nofluxes,nosteps   &
                              ,soilwater_out,evapotrans_out,sensible_out,daily_weighted_SWP_out           &
                              ,daily_weighted_soilR_out,soilevap_out,wetcanopyevap_out,potevap_out        &
                              ,netrad_out,porosity_out,sun_canopy_temperature_out,shade_canopy_temperature_out&
                              ,canopy_temperature_out,ground_heat_out,soil_conductance_out)

    use gv_veg, only: lai,gppt,system_energy_balance,soiletmm,transt_kg_s,minlwp, &
                      sensible_heat,potential_evaporation, netrad_day
    use gv_irradiance_sunshade, only: daylength
    use gv_clim, only: wetev
    use gv_scale_declarations, only: time, mol_to_g_carbon
    use gv_soil_structure, only: soil_temp, waterfrac, watericemm, &
                                 rooted_layers, porosity
    use gv_daily_averages

    implicit none

    ! declare input variables
    ! declare input variables
    integer, intent(in) :: counter  & ! simulation step
                          ,step     & ! model time in day
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nosteps    ! number of steps in simulation

    double precision, dimension(nosteps), intent(inout) :: lai_out & ! leaf area index
                                               ,daily_weighted_SWP_out   & 
                                               ,daily_weighted_soilR_out &
                                               ,sun_canopy_temperature_out & 
                                               ,shade_canopy_temperature_out &
                                               ,canopy_temperature_out &
                                               ,ground_heat_out &
                                               ,soil_conductance_out &
                                               ,evapotrans_out & ! evapotranspiration (mm.day-1)
                                               ,soilevap_out &   ! soil evaporation (mm.day-1) 
                                               ,wetcanopyevap_out & ! wet canopy evaporation (mm.day-1)
                                               ,potevap_out & ! potential evaporation (kg.m-2.day-1)
                                               ,netrad_out  & ! net radiation (W.m-2)
                                               ,porosity_out & ! mean porosity in rooting zone
                                               ,sensible_out & ! sensible heat
                                               ,soilwater_out & ! soil water content (mm)
                                               ,GPP_out   ! Gross primary productivity

    double precision, dimension(nosteps,nopools), intent(inout) :: POOLS ! vector of pools
    double precision, dimension(nosteps,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare local variables
    integer :: i, dusk, dawn
    double precision :: steps_per_day_real

    steps_per_day_real = dble(time%steps_per_day)
    dawn = nint((steps_per_day_real * 0.5d0)-(daylength*0.5d0))    
    dusk = nint((steps_per_day_real * 0.5d0)+(daylength*0.5d0))

    ! prepare the output
    lai_out(counter) = dble(sum(lai))
    ! umolC.m-2.s-1 -> gC.m-2.t-1
    GPP_out(counter) = dble(sum(gppt*time%seconds_per_step * 1d-6 * mol_to_g_carbon))
    FLUXES(counter,1) = dble(sum(system_energy_balance))
    soilwater_out(counter) = dble((watericemm(1)))
    evapotrans_out(counter) = dble(sum(transt_kg_s*time%seconds_per_step)) ! add transpiration
    evapotrans_out(counter) = evapotrans_out(counter)+dble(sum(wetev+soiletmm)) ! add wet surface
    soilevap_out(counter) = dble(sum(soiletmm)) !explicit soil evaporation out
    wetcanopyevap_out(counter) = dble(sum(wetev)) ! explicit wet canopy surface out
    potevap_out(counter) = dble(sum(potential_evaporation))
    netrad_out(counter) = dble(sum(netrad_day))/steps_per_day_real ! mean net radiation (W.m-2)
    sensible_out(counter) = dble(sum(sensible_heat))*time%seconds_per_step*1d-6 ! W.m-2 -> MJ.m-2.day-1
    ground_heat_out(counter) = dble(sum(daily_ground_heat_flux*time%seconds_per_step)*1d-6)
    daily_weighted_SWP_out(counter) = dble(sum(daily_weighted_SWP)/steps_per_day_real)
    daily_weighted_soilR_out(counter) = dble(sum(daily_weighted_soilR)/steps_per_day_real)
    sun_canopy_temperature_out(counter) = dble(sum(daily_sun_canopy_temperature(dawn:dusk))/dble(nint(daylength)))
    shade_canopy_temperature_out(counter) = dble(sum(daily_shade_canopy_temperature(dawn:dusk))/dble(nint(daylength)))
    canopy_temperature_out(counter) = dble(sum(daily_canopy_temperature(dawn:dusk))/dble(nint(daylength)))
    porosity_out(counter) = dble(porosity(1)) 
    soil_conductance_out(counter) = dble(sum(daily_soil_conductance)/steps_per_day_real)

  end subroutine spa_couple_output_daily
  !
  !-------------------------------------------
  !
  subroutine spa_couple_output_hourly(counter,step,GPP_out,lai_out,FLUXES,POOLS,nopools,nofluxes,nosteps &
                              ,soilwater_out,evapotrans_out,sensible_out,daily_weighted_SWP_out          &
                              ,daily_weighted_soilR_out,soilevap_out,wetcanopyevap_out,potevap_out       &
                              ,netrad_out,porosity_out,sun_canopy_temperature_out,shade_canopy_temperature_out &
                              ,canopy_temperature_out,ground_heat_out,soil_conductance_out)

    use gv_veg, only: lai, gppt, system_energy_balance, soiletmm, transt, minlwp, &
                      sensible_heat, potential_evaporation, netrad_day
    use gv_clim, only: wetev, lambda_bulk
    use gv_scale_declarations, only: time, mol_to_g_carbon
    use gv_soil_structure, only: soil_temp, waterfrac, watericemm, &
                                 rooted_layers, porosity
    use gv_daily_averages

    implicit none

    ! declare input variables
    ! declare input variables
    integer, intent(in) :: counter  & ! simulation step
                          ,step     & ! model time in day
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nosteps    ! number of steps in simulation

    double precision, dimension(nosteps), intent(inout) :: lai_out & ! leaf area index
                                               ,daily_weighted_SWP_out   &
                                               ,daily_weighted_soilR_out &
                                               ,sun_canopy_temperature_out &
                                               ,shade_canopy_temperature_out &
                                               ,canopy_temperature_out &
                                               ,ground_heat_out &
                                               ,soil_conductance_out &
                                               ,evapotrans_out & ! evapotranspiration (W.m-2)
                                               ,soilevap_out & ! soil evaporation (W.m-2)
                                               ,wetcanopyevap_out & ! wet canopy evaporation (W.m-2)
                                               ,potevap_out & ! potential evaporation (W.m-2)
                                               ,netrad_out & ! net radiation (W.m-2)
                                               ,porosity_out & ! mean porosity in rooting zone
                                               ,sensible_out   & !  sensible heat
                                               ,soilwater_out & ! soil water content (mm)
                                               ,GPP_out   ! Gross primary productivity

    double precision, dimension(nosteps,nopools), intent(inout) :: POOLS ! vector of pools
    double precision, dimension(nosteps,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare local variables
    integer :: i

    ! prepare the output
    lai_out(counter) = dble(sum(lai))
    ! umolC.m-2.s-1
    GPP_out(counter) = dble(gppt(step))
    FLUXES(counter,1) = dble(system_energy_balance(step))
    soilwater_out(counter)=dble(waterfrac(1))!dble(sum(watericemm(1:rooted_layers)))
    evapotrans_out(counter)=dble(transt(step)) ! add transpiration
    evapotrans_out(counter)=evapotrans_out(counter)+dble(((wetev(step)+soiletmm(step))/time%seconds_per_step)*lambda_bulk) ! add wet surf & soil
    soilevap_out(counter)=dble((soiletmm(step)/time%seconds_per_step)*lambda_bulk) 
    wetcanopyevap_out(counter)=dble((wetev(step)/time%seconds_per_step)*lambda_bulk)
    potevap_out(counter)=dble((potential_evaporation(step)/time%seconds_per_step)*lambda_bulk)
    netrad_out(counter)=dble(netrad_day(step))
!    evapotrans_out(counter)=dble((wetev(step)/time%seconds_per_step)*lambda_bulk)
    sensible_out(counter)=dble(sensible_heat(step))
    ground_heat_out(counter) = dble(daily_ground_heat_flux(step))
    daily_weighted_SWP_out(counter)=dble(daily_weighted_SWP(step))
    daily_weighted_soilR_out(counter)=dble(daily_weighted_soilR(step))
    sun_canopy_temperature_out(counter) = dble(daily_sun_canopy_temperature(step))
    shade_canopy_temperature_out(counter) = dble(daily_shade_canopy_temperature(step))
    canopy_temperature_out(counter) = dble(daily_canopy_temperature(step))
    porosity_out(counter) = dble(sum(porosity(1:rooted_layers))) / dble(rooted_layers)
    soil_conductance_out(counter) = dble(daily_soil_conductance(step))

  end subroutine spa_couple_output_hourly
  !
  !-------------------------------------------
  !
  subroutine spa_couple_time_update(time,doy_start)

    ! Updates parts of the time-holder variable. !
    ! to account for the coupling situation
    
    use gv_scale_declarations, only: time_holder
    use spa_io, only: prev_step_of_day

    implicit none

    ! arguments..
    type(time_holder),intent(inout) :: time ! datetime object
    integer, intent(in) :: doy_start

    ! set the beginning of the day variable to the DOY beginning
    time%day  = int(doy_start)
    time%run_day = 0
    time%year = 0
    time%step = 0
    time%steps_count = 0
    time%daytime = 0d0
    prev_step_of_day = -9999

  end subroutine spa_couple_time_update
  !
  !-------------------------------------------
  !
end module main_spa

