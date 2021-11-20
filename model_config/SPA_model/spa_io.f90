! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_io

  !! coordinates calling of all io routines, including    !!
  !! reading user-config and opening/reading input files. !!

  use gv_scale_declarations, only: fname_length
#ifdef USE_NETCDF
  use spa_io_netcdf,         only: nc_met_file
#endif

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: handle_output, start_spa, update_met_drivers, update_phenology_drivers &
           ,met_slice, pheno_slice, calculate_simulation_days, prev_step_of_day

#ifdef USE_NETCDF
  ! Not sure this is the right place to put this..
  type(nc_met_file),pointer :: met_nc
#endif

  ! variables.. 
  character(fname_length) :: spa_config_filename = "SPA.config" ! default value
  integer                 :: pheno_slice = 0,  &
                             met_slice = 0,    &
                             prev_step_of_day

  save
  
contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine handle_output( flag , time , output_data )

    ! wrapper to the various possible writes. !

    use gv_scale_declarations, only: time_holder, user_opts
    use log_tools
    use spa_io_csv
#ifdef USE_NETCDF
     use spa_io_netcdf
#endif

    implicit none

    integer,intent(in)      :: flag
    ! flag takes the following values...
    ! 0   = open output files
    ! 1   = Ecosystem fluxes (model time step)
    ! 2   = Daily time step 
    ! 3-7 = specific writes for subroutines within SPA
    ! 8   = crops
    ! 9   = close output files
    type(time_holder),optional,intent(in) :: time
    double precision,dimension(:),optional            :: output_data

    if (user_opts%std_csv_output) then

       if ( ( flag .ge. 3 ) .and. ( flag .le. 7) &
            .and. ( .not. present(output_data) ) ) then
         write(message,*)"When handle_output is called with flag=",flag,&
                           " you MUST supply output_data."
         call write_log( trim(message), msg_fatal, __FILE__ , __LINE__ )
       end if

       select case (flag)
       case (0) !open output files and write headers..
           call open_output_csv( user_opts%output_directory , user_opts%plant_func_type )
       case (1) ! write the standard output..
           call write_ecosystem_fluxes_output_csv( time )
       case (2)
           call write_daily_stocks_output_csv( time )
           call write_daily_fluxes_output_csv( time , user_opts%plant_func_type )
           call write_daily_water_output_csv( time )
           call write_daily_ecosystem_fluxes_output_csv( time )
       case (3)
           call write_assimilate_output_csv( time , output_data , .true. )
       case (4)
           call write_assimilate_output_csv( time , output_data , .false. )
       case (5)
           call write_solar_output_csv( output_data )
       case (6)
           call write_soil_output_csv( .false. , output_data )
       case (7)
           call write_soil_output_csv( .true. , output_data )
       case (8)
           call write_crops_output_csv
       case (9)
           call close_output_csv
       case default
           write(message,*)"flag supplied to handle_output:",flag," was not recognised"
           call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
       end select

    endif ! full output

  end subroutine handle_output
  !
  !----------------------------------------------------------------------
  !
  subroutine start_spa

    ! Read the user config, then call the relevant routines !
    ! to open files and do any initialisation required.     !

    use config_tools,          only: ConfigSection
    use gv_scale_declarations, only: met, pheno, user_opts
    use log_tools
    use spa_cmd_line,          only: parse_spa_cmd_line
    use spa_config,            only: read_user_config, update_parameters_from_user_config
    use spa_initialise
    use spa_io_csv
#ifdef USE_NETCDF
    use spa_io_netcdf
    use spa_restart
#endif

    implicit none

    ! local variables..
    type(ConfigSection), pointer :: section_names

    ! Find out the config filename..
    call parse_spa_cmd_line( spa_config_filename )

    ! Read user's config file..
#ifdef USE_NETCDF
!    allocate( met_nc, met_nc%header )
    allocate( met_nc )
    call read_user_config( spa_config_filename , section_names , met_nc )
#else
    call read_user_config( spa_config_filename , section_names )
#endif

    ! Open the meteorology file..(but don't load any data from it yet)
    if ( user_opts%met_file_is_nc ) then
#ifdef USE_NETCDF
       ! load basic data associated with file...
       call write_log( "opening the NC meteorology input file" , msg_info , __FILE__ , __LINE__ )
       call open_met_nc( met_nc%header , met_nc%time , met_nc%lat , met_nc%lon )

       ! ! Check latitude and longitude desired by user are inside bounds of this file..
       ! call write_log( "checking it contains the desired lat-lon" , msg_info , __FILE__ , __LINE__ )
       ! call check_LatLon( met_nc%lat%values, met_nc%lon%values, met_nc%grid )

       ! Load all the meteorological data into the pointer..
       call load_met_nc_data( met_nc )
       call write_log_div
#else
       write(message,*) "SPA was NOT compiled with the Netcdf library, so cannot open the NC " &
                      //"meteorology input file.  See the Makefile for the switch to change this." 
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
#endif
    else

       ! Load all the meteorological data into the pointer..
       call write_log( "opening the CSV meteorology input file" , msg_info , __FILE__ , __LINE__ )
       call read_met_csv( user_opts%met_filename , met )
       call write_log_div

       ! also now read phenology file if needed
       if (.not. user_opts%use_dalec) then
          call read_phenology_csv( user_opts%phenology_filename , pheno )
       end if

    end if

    ! calculate which years of meteorology are leap years
    call calculate_simulation_days( met )

    ! calculate the mean air temperature of the warmest quarter, averaged across
    ! years for use in the maintenance respiration model (Reich et al., 2008;
    ! Thomas et al., in prep)
!    call calculate_air_temperature_of_warmest_quarter( met )

    if ( .not. user_opts%loop_pheno_driver .and. .not. user_opts%use_dalec ) then
      ! check the phenology-driver file has enough timeslices to perform the run..
      write(message,*) "checking phenolgy file contains enough time-slices for run length"
      call write_log( message , msg_info , __FILE__ , __LINE__ )
      call check_enough_pheno_timeslices( pheno )
    end if

    if ( .not. user_opts%loop_met_driver ) then
      ! check the met-driver file has enough timeslices to perform the run..
      write(message,*) "checking met file contains enough time-slices for run length"
      call write_log( message , msg_info , __FILE__ , __LINE__ )
      call check_enough_timeslices( met )
    end if

#ifdef USE_NETCDF
    ! Setup lists of SPA variables/dimensions
    ! (useful in various input/output files)..
    call setup_spa_lists
#endif

    ! How are we starting model?..
    if ( user_opts%load_restart ) then

#ifdef USE_NETCDF
      ! Restarting from a previous dump..
      call write_log( "Starting SPA from restart file." , msg_info , __FILE__ , __LINE__ )
      call restart_model( user_opts%restart_filename )
#else
       write(message,*) "SPA was NOT compiled with the Netcdf library, so cannot create or read" &
                      //" from restart files.  See the Makefile for the switch to change this." 
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
#endif

    else

      ! Standard initialisation..

      ! Read the soils file..    
      call write_log( "opening the soils file" , msg_info , __FILE__ , __LINE__ )
      call read_soil_csv( user_opts%soils_filename, user_opts%plant_func_type )

      ! Read the veg file..
      call write_log( "opening the veg file" , msg_info , __FILE__ , __LINE__ )
      call read_veg_csv( user_opts%veg_filename )

      ! Read the carbon file..
      if (user_opts%use_dalec) then
          call write_log( "opening the carbon file" , msg_info , __FILE__ , __LINE__ )
          call read_carbon_csv( user_opts%carbon_filename )
      endif

      ! Initialise the soils..
      call write_log( "Initialising the soils" , msg_info , __FILE__ , __LINE__ )
      call initialise_soils

      ! Initialise the veg.. (On day 1 set leaf WP in each layer)
      call write_log( "Initialising the veg" , msg_info , __FILE__ , __LINE__ )
      call initialise_veg

      ! If this is a crops run..
      if ( user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3 ) then

        ! Read the crops file..
        call read_crops_csv( user_opts%crops_filename )

        ! initialise the crops..
        call initialise_crops

      end if

    end if

    ! Make sure parameters read in from the user-config overwrite any in the files/default values..
    call write_log( "Reading user-defined params from user-config" , msg_info , __FILE__ , __LINE__ )
    call update_parameters_from_user_config( section_names )

    ! Open output files..
    call write_log( "Opening output files" , msg_info , __FILE__ , __LINE__ )
    call handle_output( 0 )

  end subroutine start_spa
  !
  !----------------------------------------------------------------------
  !
  subroutine update_met_drivers( time )

    ! wrapper to call appropriate nc/csv routine !
    ! for reading next slice of meteorology data !

    use gv_clim,               only: atmos_press, coa, daypar, dayppt, par_top, ppt, sw_rad, lw_rad, sw_diffuse, &
                                     temp_bot, temp_top, vpd_bot, vpd_top, wdbot, wdtop, wetev, wind_spd, &
                                     snowfall, gas_constant_d, cp_air, swc10, swc20, swc30, swc40, swc50, &
                                     swc60, swc70, swc80, swc90, swc100, swc110, swc120, swc130, swc140, swc150
    use gv_meteo,              only: air_density_kg, abs_pot_conv, abs_virt_conv
    use gv_hourscale,          only: discharge, runoff, freeze
    use gv_scale_declarations, only: met, time_holder, user_opts, boltz
    use gv_veg,                only: flux, gppt, respt, totass, totevap, totres, transt, transt_kg_s, &
                                     sensible_heat,latitude_radians
    use gv_carbon_model,       only: clearance_fraction, clearance_management, fire_fraction, &
                                     avg_min_airt_store, avg_dayl_store, avg_vpd_Pa_store, &
                                     avg_min_airt,avg_dayl,avg_vpd_Pa
    use math_tools,            only: calculate_daylength_hours
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time 

    ! local variables
    integer :: i
    double precision :: vpsat, & ! Saturation vapour pressure (Pa)
            vpair, & ! Vapour pressure of the air (Pa)
           vpd_pa, & ! vapour pressure deficit (Pa)
         sp_moist    ! specific humidity (ratio; kg/kg)

    ! check we haven't been called already this step...
    if ( time%step .ne. prev_step_of_day ) then

      ! iterate forward..
      met_slice = met_slice + 1

      if ( user_opts%loop_met_driver ) then
          ! if exceeding length of data then reset to start of data..
          if ( met_slice .gt. size(met%temp,1) ) met_slice = 1
      end if

      ! On first hour of day, set counters to zero..
      if ( time%step .eq. 1 ) then
          daypar  = 0d0 ; dayppt  = 0d0 ; discharge = 0d0 ; flux   = 0d0
          gppt    = 0d0 ; respt   = 0d0 ; runoff    = 0d0 ; totass = 0d0
          totevap = 0d0 ; totres  = 0d0 ; transt    = 0d0 ; wetev  = 0d0
          transt_kg_s = 0d0 ; sensible_heat = 0d0 
      end if

      ! load this time-step's driver into each var..
      atmos_press = met%sfc_pressure( met_slice )
      coa         =  met%ambient_co2( met_slice ) 
      par_top     =         met%ppfd( met_slice )
      ppt         =       met%precip( met_slice )
      sw_rad      =       met%sw_rad( met_slice )
      temp_bot    =         met%temp( met_slice ) 
      temp_top    = temp_bot
      vpd_bot     =          met%vpd( met_slice )
      vpd_top     = vpd_bot
      wind_spd    =     met%wind_spd( met_slice )
      swc10       =  met%swc10( met_slice )
      swc20       =  met%swc20( met_slice )
      swc30       =  met%swc30( met_slice )
      swc40       =  met%swc40( met_slice ) 
      swc50       =  met%swc50( met_slice )
      swc60       =  met%swc60( met_slice )
      swc70       =  met%swc70( met_slice )
      swc80       =  met%swc80( met_slice )
      swc90       =  met%swc90( met_slice )
      swc100      =  met%swc100( met_slice )
      swc110      =  met%swc110( met_slice )
      swc120      =  met%swc120( met_slice )
      swc130      =  met%swc130( met_slice )
      swc140      =  met%swc140( met_slice )
      swc150      =  met%swc150( met_slice )

      ! calculate 21 day averaged GSI inputs if needed
      if (user_opts%use_dalec .and. &
          user_opts%plant_func_type == 0 .or. user_opts%plant_func_type == 1) then
          ! first allocate the number of days needed
          if (.not.allocated(avg_min_airt)) then
              allocate(avg_min_airt(1),avg_dayl(1),avg_vpd_Pa(1) &
                      ,avg_min_airt_store(21),avg_dayl_store(21),avg_vpd_Pa_store(21))
              avg_min_airt_store = temp_top ; avg_vpd_Pa = vpd_top*1d3 
              avg_dayl_store = calculate_daylength_hours(dble(time%day),latitude_radians)
              ! convert to seconds
              avg_dayl_store = avg_dayl_store * 3600d0
          else 
              ! at the beginning of the day we need to shift averaging values
              ! along and reset some variables
              if (time%step == 1) then
                  ! if we are at the beginning of the day shift the averaging
                  ! vectors down one
                  do i = 20, 1, -1
                     avg_min_airt_store(i+1) = avg_min_airt_store(i)
                     avg_dayl_store(i+1) = avg_dayl_store(i)
                     avg_vpd_Pa_store(i+1) = avg_vpd_Pa_store(i)
                  enddo
                  ! as daylength doesn't change each day we do the calculation
                  ! one per day here
                  avg_dayl_store(1) = calculate_daylength_hours(dble(time%day),latitude_radians)
                  avg_dayl_store(1) = avg_dayl_store(1)*3600d0
                  ! not forgetting to reset the temperature and vpd components
                  ! for correct min and averaging respectivily
                  avg_min_airt_store(1) = temp_top
                  avg_vpd_Pa_store(1) = (vpd_top*1d3) / time%steps_per_day
              else
                  ! now we need to begin generating the 21 day averages for
                  ! avg_min_airt and avg_vpd_Pa
                  ! is the current 
                  avg_min_airt_store(1) = min(temp_top,avg_min_airt_store(1))
                  avg_vpd_Pa_store(1) = avg_vpd_Pa_store(1) &
                                      + ((vpd_top*1d3) / time%steps_per_day)
              endif ! beginning or not
              ! at the end of the day generate the new average values
              ! This update MUST occur at the same frequency as the GSI call.
              if (time%step == time%steps_per_day) then
                  ! if we are at the end of the day calculate the new 21 day averages
                  ! for use is the carbon_model.f90 (0.04761905 = 1/21)
                  avg_min_airt = sum(avg_min_airt_store)*0.04761905d0
                  avg_dayl = sum(avg_dayl_store)*0.04761905d0
                  avg_vpd_Pa = sum(avg_vpd_Pa_store)*0.04761905d0
              endif ! end of the day
          endif ! first time or not
      endif ! use dalec

      ! if snowfall information has been provided
      if (user_opts%met_file_has_snowfall) then
          snowfall = met%snowfall(met_slice)
      else
          if (temp_top <= 0d0) then
             ! converting rain to rate
             snowfall = ppt / time%seconds_per_step ; ppt = 0d0
          else
             snowfall = 0d0
          endif
      endif

      ! if disturbance information is available      
      if (user_opts%met_file_has_disturbance) then
          clearance_fraction = met%clearance_fraction( met_slice )
          clearance_management = met%clearance_management( met_slice )
          fire_fraction = met%fire_fraction( met_slice )
      else 
          clearance_fraction = 0
          clearance_management = 0
          fire_fraction = 0
      endif

      ! if diffuse short wave was in the file calculate the fraction for use in
      ! solar*
      if (user_opts%met_file_has_sw_diffuse) then
          sw_diffuse = met%sw_diffuse( met_slice ) / sw_rad
          if (sw_diffuse > 1) then 
              print*,"diffuse short wave radiation cannot be greater than total short wave radiation"
              stop
          endif
      endif

      ! If incoming longwave radiation wasn't in the met-file
      !  then calculate it according to Tair..
      if ( user_opts%met_file_has_lw_rad ) then
          lw_rad = met%lw_rad( met_slice )
      else
         ! WARNING! This formulation doesn't really work but this is the best we can do.
         ! Suggest always using measured lw_rad if at all possible.
          lw_rad = boltz * ( temp_top + freeze - 20d0 )**4
      endif

      ! calculate new air density (Kg m-3)
      air_density_kg = 353d0 / (temp_top + 273.15d0 )
      ! calculate ratio for conversion between absolute and potential
      ! temperatures (K)
      abs_pot_conv = (100000d0/atmos_press)**(gas_constant_d/cp_air)

      ! Saturation vapour pressure (Pa) calculation from Jones p110; uses
      ! absolute air temperature (oC)
      vpsat = (613.75d0*exp((17.502d0*temp_top)/(240.97d0+temp_top)))

      ! negative VPDs are not possible (in reality, although the instruments can
      ! measure this due to condensation) therefore set these to zero
      if (vpd_top < 0d0) then
          vpd_top = 0d0 ; vpd_bot = 0d0
          call write_log( "Negative VPD value found an replaced with zero" , msg_info , __FILE__ ,__LINE__ )    
      endif

      if (user_opts%vpd_or_specific_humidity == "specific_humidity") then
         ! Therefore vpd_top is actually specific humidity (kg/kg)
         sp_moist = vpd_top
         ! Determine vapour pressure (Pa); based on specific humidity, air
         ! pressure
         ! (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
         ! a physical outline)
         vpair = ((sp_moist*(atmos_press*1.0d-2))/0.62197d0)*1.0d2
         ! Difference between the vapour pressure of saturation and air, i.e. the
         ! VPD (Pa)
         vpd_pa = vpsat-vpair
         ! Vapour pressure (Pa) of the air determined from saturation VP and humidity
         ! VPD is SPA is used as kPa, Pa -> kPa  
         vpd_top = vpd_pa*1d-3 ; vpd_bot = vpd_top
      else if (user_opts%vpd_or_specific_humidity == "vpd") then
         ! Calculate vapour pressure of air (Pa) from VPD; vpd_top (kPa -> Pa), 
         vpair = vpsat - (vpd_top * 1d3)
         ! calculate specific humidity; 0.62197 is ratio of gas constants (p95
         ! McIlveen 1986, Basic Meteorology - a physical outline)
         sp_moist = ((vpair / 1d2)*0.62197d0)/(atmos_press*1d-2)
      end if ! specific_humidity or vpd

      ! calculate virtual (water vapour corrected) temperatures for air and leaf
      ! (K)
      abs_virt_conv = (1d0 + 0.61d0 * sp_moist)

      ! Calculate the absolute water deficit (g m-3)
      wdtop = vpd_top * 217.0d0 / ( 0.1d0 * ( temp_top + 273.4d0 ) )
      wdbot = vpd_bot * 217.0d0 / ( 0.1d0 * ( temp_bot + 273.4d0 ) )

      ! Expectation is that saturated mixing ratio and therefore subsequently
      ! calculated vpsat may be more accurate than Jones equation.
      ! However according to NOAH scheme the QGH variable from WRF breaks down
      ! if surface is covered with snow. 

      ! Saturation vapour pressure (Pa) calculation from Jones p110; uses
      ! absolute air temperature (oC)
!      vpsat=(0.061375*exp((17.502*abs_temp)/(240.97+abs_temp)))*1.0e4

      ! Need the saturation atmospheric mixing ratio (kg.kg-1). This equation
      ! requires pressure values as (hPa); so *1.e-2.
      ! This equation then provides (g.kg-1); we need (kg.kg-1) so divide by
      ! 1000 
      ! vpsat and atmospheric pressure values are in (Pa)
      ! Replaced 21/01/2011 with saturation mixing ratio from WRF
      ! sat_mix=(621.97*((vpsat*1.0e-2)/((pressure*1.0e-2)-(vpsat*1.0e-2))))*1.0e-3

      ! Calculate relative humidity as a proportion, not %, as this is what is
      ! required for the VPD calculation
      !rel_humidity=moist_mix/sat_mix

      ! if SNOW conditions
      ! Calculate saturation specific humidity from saturation mixing ratio
      ! (kg/kg)
!       sat_sp_moist=water_saturation_ratio/(1.0+water_saturation_ratio))
      ! Determine saturation vapour pressure (Pa); based on saturation specific
      ! humidity, air pressure (hPa) and ratio of gas constants (p95 McIlveen
      ! 1986, Basic Meteorology - a physical outline)
!      vpsat=((sat_sp_moist*(pressure*1.0e-2))/0.62197)*1.0e2
      !---------
      ! Calculate specific humidity
 !     sp_moist=water_mixing_kg/(1.0+water_mixing_kg))
      ! Determine vapour pressure (Pa); based on specific humidity, air pressure
      ! (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
      ! a physical outline)
 !     vpair=((sp_moist*(pressure*1.0e-2))/0.62197)*1.0e2
      ! Diagnostic tool
 !     rel_humidity=vpair/vpsat
      ! Difference between the vapour pressure of saturation and air, i.e. the
      ! VPD (Pa)
 !     vpd_pa=vpsat-vpair
 !     ! Vapour pressure (Pa) of the air determined from saturation VP and
 !     humidity
      ! VPD is SPA is used as kPa, Pa -> kPa  
 !     vpdtop=vpd_pa*1e-3
 !     ! Vapour pressure deficit at canopy top (kPa)
 !     vpdbot=vpd_pa*1e-3
      ! Diagnostic tool; VPD (kPa)
!      wrf_spa_vpd(X,Y)=vpdtop
      ! Convert VPD (kPa) into absolute water deficit (g.m-3)
!      wdtop=vpdtop*217./(0.1*(abs_temp+273.15))
!      wdbot=vpdbot*217./(0.1*(abs_temp+273.15))

      !--- The following are only used for output files.. ---

      ! keep a (daily) running total of the precip..
      dayppt = dayppt + ppt

      ! sum the daily energy (PAR multiplied by nos of seconds per timestep)..
      daypar = daypar + par_top * time%seconds_per_step

      prev_step_of_day = time%step

!      call write_log( "Loaded met data for step" , msg_info , __FILE__ , __LINE__ )

    end if ! are we at the beginning of a new day?

  end subroutine update_met_drivers
  !
  !--------------------------------------------------------------------
  !
  subroutine update_phenology_drivers

    ! wrapper to call appropriate nc/csv routine !
    ! for reading next slice of meteorology data !

    use gv_scale_declarations, only: pheno, time, user_opts, dble_one, dble_zero
    use gv_veg,                only: lai, stock_roots, lafrac, &
                                     LCA, avN, nfrac, Nla
    use gv_carbon_model,       only: pars,leaf_growth, leaf_death
    use log_tools

    implicit none

    ! local variables
    double precision :: leaf_change

    ! iterate forward..
    pheno_slice = pheno_slice + 1
  
    if ( user_opts%loop_pheno_driver ) then
        ! if exceeding length of data then reset to start of data..
        if ( pheno_slice .gt. size(pheno%lai,1) ) pheno_slice = 1
    end if
!  pheno%lai(pheno_slice) = 1.0 ! pheno%lai(pheno_slice) * 3.0
    if (time%steps_count == 0) then
        ! allocate initialisation (for consistency across all runs) leaf area across layers 
!        lai = lafrac*dble_one
!        ! estimate initial LCA based on maintaining C:N ratio through canopy
!        LCA = 0.0 ; Nla = (nfrac * avN * sum( lai )) !; nfrac_top = nfrac(1)
!        where (lai > 0.0) Nla = Nla / lai ! avN per canopy layer (gN.m-2)
!        where (lai > 0.0) LCA = pars(17) * (Nla / avN)
!        ! allocate actual leaf area across layers
        lai = pheno%lai( pheno_slice )*lafrac 
        leaf_growth = dble_zero ; leaf_death = dble_zero
    else 
        leaf_change = pheno%lai( pheno_slice )-sum(lai)
        if (leaf_change > 0d0) leaf_growth = leaf_change*pars(17)
        if (leaf_change < 0d0) leaf_death = -leaf_change*pars(17)
    endif

    ! load root carbon
    stock_roots = pheno%rootC( pheno_slice )

!    call write_log( "Loaded phenology data for day" , msg_info , __FILE__ ,__LINE__ )

  end subroutine update_phenology_drivers
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
    subroutine calculate_simulation_days( met )

    use gv_scale_declarations, only: met_drivers, time, user_opts, days_in_month

    implicit none

    ! arguments
    type(met_drivers),intent(in) :: met

    ! local variables..
    integer :: year_start,year, month, day, i, nos_of_timeslices, approx_nos_years
    double precision :: dummy

    ! calculate how much data we have been given
    nos_of_timeslices = size( met%temp )
    ! calculate the approximate number of years in the met file.
    ! this should be good enough to work out the sequence of normal to leap
    ! years
    dummy = nos_of_timeslices / (time%steps_per_day * 365.25)
    approx_nos_years = nint(dummy)

    ! read in the start datetime and extract the date components
    read(time%start_date(1:4),  fmt='(I4)') year
    read(time%start_date(6:7),  fmt='(I2)') month
    read(time%start_date(9:10), fmt='(I2)') day

    if (user_opts%include_leap_years) then
       ! Never forget where you began
       year_start=year
       ! loop through years to determine which are leap years
       do i = 1, time%nos_of_years
          ! current year a leap or not
          time%days_in_year(i) = 365
          if (mod(year,4).eq.0) then
             time%days_in_year(i) = 366
             if (mod(year,100).eq.0) then
                time%days_in_year(i) = 365
                if (mod(year,400).eq.0) then
                   time%days_in_year(i) = 366
                endif
              endif
           endif
          ! update year to the next one for checking
          year=year+1
          if ((year-year_start) == approx_nos_years) year=year_start
       end do

    endif ! end if for leap year condition

    ! read in the start datetime and extract the date components
    read(time%start_date(1:4),  fmt='(I4)') year
    read(time%start_date(6:7),  fmt='(I2)') month
    read(time%start_date(9:10), fmt='(I2)') day

    ! adjust for leap year of one has been specified
    if (time%days_in_year(1) == 366 .and. month > 2) day = day + 1 
    ! use the start date information to update the starting day
    time%day = (sum(days_in_month(1:month)) - days_in_month(month)) + day

    ! now if this is a run of less than one year and we are not using leap
    ! years...
    if (time%nos_of_years == 1 .and. time%days_in_year(1) < 365 &
       .and. .not.user_opts%include_leap_years) then
       time%days_in_year=time%day+(time%days_in_year-1)
    endif

  end subroutine calculate_simulation_days
  !
  !----------------------------------------------------------------------
  !
!  subroutine calculate_air_temperature_of_warmest_quarter(met)
!
!    ! subroutine calculates the mean air temperature (oC) for the warmest
!    ! quarter of the year. This TWQ is then averaged across all years in the
!    ! simulation period. The TWQ is used in the adjustment of scaler baseling
!    ! parameter linking foliar N to maintenance respiration (Reich et al 2008; Quinn et al 2017). 
!    ! Historical values could be used if it is believed that plant metabolism is sensitive on shorter time scales.
!
!    use gv_carbon_model, only: twq
!    use gv_scale_declarations, only: met_drivers, time
!
!    ! arguments
!    type(met_drivers),intent(in) :: met
!
!    ! local variables
!    integer :: d, yr, &
!               nos_of_timeslices, &
!               nos_slices_in_quarter, &
!               est_nos_years_in_drivers, &
!               start,finish,start_offset
!    double precision :: tmp_store,tmp
!    double precision, dimension(:), allocatable :: twq_tmp
!
!    ! calculate how much data we have been given
!    nos_of_timeslices = size( met%temp )
!
!    ! estimate the number of years in the actual driver file; here is is 
!    est_nos_years_in_drivers = floor((dble(nos_of_timeslices)/dble(time%steps_per_day))/365.25d0)
!    ! estimate how many time steps per quarter
!    nos_slices_in_quarter = floor(90d0*time%steps_per_day)
!    allocate(twq_tmp(est_nos_years_in_drivers))
!    start_offset = 0 ; tmp_store = -9999d0
!    ! loop through each year
!    do yr = 1, est_nos_years_in_drivers
!       ! loop through a rolling window of the air temperature one quarter in
!       ! length
!       do d = 1,time%days_in_year(yr)
!          start=(start_offset+d)*time%steps_per_day ; finish = start + nos_slices_in_quarter
!          tmp = sum(met%temp(start:finish))/dble(nos_slices_in_quarter)
!          if (tmp > tmp_store) tmp_store = tmp
!          ! if finish has hit the end of the available data then escape the loop
!          if (finish == nos_of_timeslices ) exit
!       end do !  loop rolling quarter
!       ! update offset of rolling loop
!       start_offset = start_offset + time%days_in_year(yr)
!       ! store this years warmest quarter value
!       twq_tmp(yr)=tmp_store
!       ! reset tmp_store
!       tmp_store = -9999d0
!    end do ! loop years
!
!    ! output
!    twq = sum(twq_tmp) / dble(est_nos_years_in_drivers)
!
!    ! tidy up
!    deallocate(twq_tmp)
!
!  end subroutine calculate_air_temperature_of_warmest_quarter
  !
  !----------------------------------------------------------------------
  !
  subroutine check_enough_timeslices( met )

    ! Check that there are enough timeslices in the driver  !
    ! file to run the model for the number of days the user !
    ! wants to run for.                                     !

    use gv_scale_declarations, only: met_drivers, time
    use log_tools

    implicit none

    ! arguments..
    type(met_drivers),intent(in) :: met

    ! local variables..
    integer :: nos_of_timeslices

    ! nos of available timeslices = length of time dimension.
    ! nos of required timeslices  = user_run_length in days * nos of steps per day (we read at every step)

    nos_of_timeslices = size( met%temp )

    if ( nos_of_timeslices .lt. (sum(time%days_in_year) * time%steps_per_day) ) then
      write(message,*)"There are not enough timeslices in the input file to run the model "//&
                      "for the number of days declared in the config file (default is 365)" 
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

  end subroutine check_enough_timeslices
  !
  !----------------------------------------------------------------------
  !
  subroutine check_enough_pheno_timeslices( pheno )

    ! Check that there are enough timeslices in the driver  !
    ! file to run the model for the number of days the user !
    ! wants to run for.                                     !

    use gv_scale_declarations, only: phenology_drivers, time
    use log_tools

    implicit none

    ! arguments..
    type(phenology_drivers),intent(in) :: pheno

    ! local variables..
    integer :: nos_of_timeslices

    ! nos of available timeslices = length of time dimension.
    ! nos of required timeslices  = user_run_length in days 

    nos_of_timeslices = size( pheno%lai )

    if ( nos_of_timeslices .lt. sum(time%days_in_year) ) then
      write(message,*)"There are not enough timeslices in the input file to run the model "//&
                      "for the number of days declared in the config file (default is 365)"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

  end subroutine check_enough_pheno_timeslices
  !
  !----------------------------------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine check_LatLon( file_lats , file_lons , file_grid )

    ! Check the lat/lon defined in the user-config falls within !
    ! that of the input file. If they are, retrieve the indices !
    ! of the points that encompass the user's desired location. !

    use gv_scale_declarations, only: grid
    use log_tools
    use netcdf_tools
    use spa_io_netcdf,         only: gridcalc

    implicit none

    ! arguments..
    double precision,dimension(:),intent(in) :: file_lats, file_lons
    type(gridcalc),pointer       :: file_grid

    ! local variables..

    logical :: status 
    integer :: tmp(1), grid_index
    double precision    :: spacing, value

    status = is_value_within_dim( file_lats, grid%latitude )
    if ( .not. status ) then
      write(message,*) "Requested latitude lies outside of bounds of input file latitudes!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    status = is_value_within_dim( file_lons, grid%longitude )
    if ( .not. status ) then
      write(message,*)"Requested longitude lies outside of bounds of input file longitudes!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! If we are still going, the required lat/lon must lie within the
    ! bounds of the dimensions, so find where..

    ! Check if the user has picked a point on the longitude grid...
    ! Find the point closest..
    tmp = minloc( file_lons - grid%longitude ) 
    grid_index = tmp(1)
    value = file_lons( grid_index )
    spacing = file_lons(2) - file_lons(1) ! assume constant grid spacing
    if ( value .eq. grid%longitude ) then
       ! set this as our index..
       file_grid%indices(1) = grid_index
       file_grid%no_bilinear = .true.
    else
       ! if not, find the nearest point just less than it...
       tmp = minloc( abs( file_lons - grid%longitude + 0.5d0*spacing ) )
       file_grid%indices(1) = tmp(1)
    end if

    ! Check if the user has picked a point on the latitude grid...
    tmp = minloc( file_lats - grid%latitude )
    grid_index = tmp(1)
    value = file_lats( grid_index )
    spacing = file_lats(2) - file_lats(1) ! assume constant grid spacing
    if ( ( value - grid%latitude ) .eq. grid%latitude ) then
       file_grid%indices(2) = grid_index
    else
       ! we only avoid bi-linear if both are true...
       file_grid%no_bilinear = .false.
       tmp = minloc( abs( file_lats - grid%latitude + 0.5d0*spacing ) )
       file_grid%indices(2) = tmp(1)
    end if

  end subroutine check_LatLon
  !
  !----------------------------------------------------------------------
  !
  logical function is_value_within_dim( dim , point )

    ! check whether a point lies within the bounds of a dimension !

    ! arguments..
    double precision,intent(in) :: dim(:), point

    ! check if value lies within it...
    if ( ( minval(dim) .lt. point ) .and. ( maxval(dim) .gt. point ) ) then
       is_value_within_dim = .true.
    else 
       is_value_within_dim = .false.
    end if

  end function is_value_within_dim
  !
  !----------------------------------------------------------------------
  !
  subroutine setup_spa_lists()

    ! Sets up a common structure which can be used by !
    ! all of SPA's input/output routines, containing  !
    ! information (and pointers to) SPA variables.    !

    use gv_hourscale,          only: canopy_store, hourts
    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid, time
    use gv_snow_info,          only: snowheight, snowweight
    use gv_soil_structure       ! loads!!
    use gv_veg                  ! loads!!
    use gv_carbon_model
    use linked_lists,       only: append, new_item, spa_dims, spa_vars, spaitem, print_item
    use log_tools
    use netcdf_tools

    implicit none

    ! Build a list of the dimensions in spa..
    if ( .not. associated( spa_dims ) ) allocate( spa_dims )
    call append( spa_dims , new_item( 'grid%core'          , grid%core          ) )
    call append( spa_dims , new_item( 'nos_canopy_layers'  , grid%canopy        ) )
    call append( spa_dims , new_item( 'nos_soil_layers'    , grid%soil          ) )
    call append( spa_dims , new_item( 'nos_wetting_laters' , grid%wetting       ) )
    call append( spa_dims , new_item( 'steps_per_day'      , time%steps_per_day ) )
    call append( spa_dims , new_item( 'time (years)'       , time%year          ) )
    call append( spa_dims , new_item( 'nopools'            , 1,grid%nopools  ) )
    call append( spa_dims , new_item( 'nofluxes'           , grid%steps_per_day,grid%nofluxes ) )
    call append( spa_dims , new_item( 'nopars'             , grid%nopars   ) )

    ! Build a list of the variables in spa..
    if ( .not. associated( spa_vars ) ) allocate( spa_vars )  
    call append( spa_vars , new_item( 'avN'                               , avn                               ) )
    call append( spa_vars , new_item( 'canopy_height'                     , canopy_height                     ) )
    call append( spa_vars , new_item( 'canopy_store'                      , canopy_store                      ) )
    call append( spa_vars , new_item( 'capac'                             , capac                             ) )
    call append( spa_vars , new_item( 'conductivity'                      , conductivity                      ) )
    call append( spa_vars , new_item( 'avtemp    '                        , avtemp                            ) )
    call append( spa_vars , new_item( 'dimen'                             , dimen                             ) )
    call append( spa_vars , new_item( 'drythick'                          , drythick                          ) )
    call append( spa_vars , new_item( 'GPP'                               , GPP                               ) )
    call append( spa_vars , new_item( 'gplant'                            , gplant                            ) )
    call append( spa_vars , new_item( 'hourts'                            , hourts                            ) )
    call append( spa_vars , new_item( 'iota'                              , iota                              ) )
    call append( spa_vars , new_item( 'kappaC'                            , kappaC                            ) )
    call append( spa_vars , new_item( 'kappaJ'                            , kappaJ                            ) )
    call append( spa_vars , new_item( 'latitude_radians'                  , latitude_radians                  ) )
    call append( spa_vars , new_item( 'max_depth'                         , max_depth                         ) )
    call append( spa_vars , new_item( 'max_storage'                       , max_storage                       ) )
    call append( spa_vars , new_item( 'minlwp'                            , minlwp                            ) )
    call append( spa_vars , new_item( 'resp_auto'                         , resp_auto                         ) )
    call append( spa_vars , new_item( 'resp_h_litter'                     , resp_h_litter                     ) )
    call append( spa_vars , new_item( 'resp_h_soilOrgMatter'              , resp_h_soilOrgMatter              ) )
    call append( spa_vars , new_item( 'rooted_layers'                     , rooted_layers                     ) )
    call append( spa_vars , new_item( 'root_k'                            , root_k                            ) )
    call append( spa_vars , new_item( 'root_radius'                       , root_radius                       ) )
    call append( spa_vars , new_item( 'rootresist'                        , rootresist                        ) )
    call append( spa_vars , new_item( 'snowheight'                        , snowheight                        ) )
    call append( spa_vars , new_item( 'snowweight'                        , snowweight                        ) )
    call append( spa_vars , new_item( 'through_fall'                      , through_fall                      ) )
    call append( spa_vars , new_item( 'tower_height'                      , tower_height                      ) )

    ! vector valued items..
    call append( spa_vars , new_item( 'pars'           , pars           , 'nopars'        ) )
    call append( spa_vars , new_item( 'conduc'         , conduc         , 'grid%core'     ) )
    call append( spa_vars , new_item( 'iceprop'        , iceprop        , 'grid%core'     ) )
    call append( spa_vars , new_item( 'lafrac'         , lafrac         , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'lafrac_dead'    , lafrac_dead    , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'lai'            , lai            , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'psil_pd'        , psil_pd        , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'Cfol_profile'   , Cfol_profile   , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'LCA'            , LCA            , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'layer_depth'    , layer_depth    , 'grid%core'     ) )
    call append( spa_vars , new_item( 'layer_height'   , layer_height   , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'LWPstore'       , LWPstore       , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'mineralfrac'    , mineralfrac    , 'grid%core'     ) )
    call append( spa_vars , new_item( 'nfrac'          , nfrac          , 'grid%canopy'   ) )
    call append( spa_vars , new_item( 'organicfrac'    , organicfrac    , 'grid%core'     ) )
    call append( spa_vars , new_item( 'root_length'    , root_length    , 'grid%core'     ) )
    call append( spa_vars , new_item( 'root_mass'      , root_mass      , 'grid%core'     ) )
    call append( spa_vars , new_item( 'thickness'      , thickness      , 'grid%core'     ) )
    call append( spa_vars , new_item( 'soil_frac_clay' , soil_frac_clay , 'grid%core'     ) )
    call append( spa_vars , new_item( 'soil_frac_sand' , soil_frac_sand , 'grid%core'     ) )
    call append( spa_vars , new_item( 'soil_temp'      , soil_temp      , 'grid%core'     ) )
    call append( spa_vars , new_item( 'waterfrac'      , waterfrac      , 'grid%core'     ) )
    call append( spa_vars , new_item( 'watericemm'     , watericemm     , 'grid%soil'     ) )
    call append( spa_vars , new_item( 'wettingbot'     , wettingbot     , 'grid%wetting'  ) )
    call append( spa_vars , new_item( 'wettingtop'     , wettingtop     , 'grid%wetting'  ) )

    ! array valued items..
    call append( spa_vars , new_item( 'POOLS'          , POOLS          , 'nopools'       ) )
    call append( spa_vars , new_item( 'FLUXES'         , FLUXES         , 'nofluxes'      ) )

 
  end subroutine setup_spa_lists
#endif
  !
  !----------------------------------------------------------------------
  !
end module spa_io
!
!------------------------------------------------------------------------
!

