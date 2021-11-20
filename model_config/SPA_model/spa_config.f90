! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_config

  !! This module reads the user-provided  !!
  !! runtime-configuration file for SPA.  !!

  use gv_scale_declarations, only: fname_length
#ifdef USE_NETCDF
  use spa_io_netcdf ! for the nc_met/nc_soils/nc_veg type definitions
#endif

  implicit none


  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: read_user_config, update_parameters_from_user_config &
           ,allocate_arrays

  interface check_allocated
    module procedure check_allocated_1d, check_allocated_2d
  end interface check_allocated

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine read_user_config( config_filename , section_names , met_nc )
#else
  subroutine read_user_config( config_filename , section_names )
#endif

    ! reads the user config, looking for the filenames !
    ! output dir and SPA configuration flags.          !

    use config_tools
    use gv_scale_declarations, only: fname_length, grid, met, time, user_opts, deg_to_rad
    use gv_veg,                only: latitude_radians
    use gv_carbon_model,       only: nopars, nofluxes, nopools, &
                                     nopars_default, nopars_crops, &
                                     nofluxes_default, nofluxes_crops, &
                                     nopools_default, nopools_crops, &
                                     pars_names,pars_names_default,pars_names_crops, &
                                     fluxes_names,fluxes_names_default,fluxes_names_crops, &
                                     pools_names,pools_names_default,pools_names_crops
    use log_tools
#ifdef USE_NETCDF
    use netcdf_tools
#endif

    ! arguments..
    character(fname_length),intent(in) :: config_filename !(input)  where to look for the config file
    type(ConfigSection), pointer       :: section_names   !(output) structure holding sections of the configuration file   
#ifdef USE_NETCDF
    type(nc_met_file),pointer          :: met_nc          !(output)
#endif

    ! local variables..
    type(ConfigSection), pointer :: section
    integer                      :: ios, dummy_integer
    character(len=100)           :: outdirtestfile

    ! Get the section names..
    call ConfigRead( config_filename, section_names )

    ! print the user's config to the logfile
    call write_log( "The contents of the configuration file are.." , msg_info , __FILE__ , __LINE__ )
    call PrintConfig( section_names, get_logunit() )
    call write_log_div

    ! Get the General options...
    call GetSection(section_names,section,'Options')
    if ( associated(section) ) then
      call GetValue(section,'canopy_MDecay_opts', user_opts%canopy_MDecay)
      if (user_opts%canopy_MDecay == 1) then
          call write_log( "Canopy momentum decay as in Williams et al 1996" , msg_info , __FILE__ , __LINE__ )
      elseif (user_opts%canopy_MDecay == 2) then
          call write_log( "Canopy momentum decay as in Smallman et al 2013" , msg_info , __FILE__ , __LINE__ )
      else
          call write_log( "Have set 'canopy_MDecay_opts' to invalid number" , msg_fatal , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'use_iWUE', user_opts%iWUE)
      if (user_opts%iWUE == 1 ) then
          call write_log( "Stomatal optimisation based on intrinsic water use efficiency (see Williams et al 1996)" &
                        , msg_info , __FILE__ , __LINE__ )
      elseif (user_opts%iWUE == 2) then
          call write_log( "Stomatal optimisation based on water use efficiency (see Bonan et al., 2014)" &
                        , msg_info , __FILE__ , __LINE__ )
      elseif (user_opts%iWUE == 3) then
          call write_log( "Stomatal optimization based on Ball-Berry" &
                        , msg_info , __FILE__ , __LINE__ )
      elseif (user_opts%iWUE == 4) then
          call write_log( "Stomatal optimization based on Medlyn" &
                        , msg_info , __FILE__ , __LINE__ )
      else
          call write_log( "Have set 'use_iWUE' to invalid state (.true. / .false.)" , msg_fatal , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'soil_hydrology_opts', user_opts%soil_hydrology_opts)
      if (user_opts%soil_hydrology_opts == 1) then
          call write_log( "Using Saxton 1986 soil equations" , msg_info , __FILE__ , __LINE__ )
      elseif (user_opts%soil_hydrology_opts == 2) then
          call write_log( "Using Hawkins",msg_info , __FILE__ , __LINE__ )
      elseif (user_opts%soil_hydrology_opts == 3) then
          call write_log( "Using Sandy Oregon soil equations" , msg_info , __FILE__ , __LINE__ )
      else
          call write_log( "Have set 'soil_hydrology_opts' to invalid number" , msg_fatal , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'prescribed_swc', user_opts%prescribed_swc)
      if (user_opts%prescribed_swc ) then
          call write_log( "using prescribed soil water content from met file" , msg_info , __FILE__ , __LINE__ )
      else
          call write_log( "Prognostic soil water content" , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'use_dalec'  , user_opts%use_dalec )
      if (user_opts%use_dalec) then
        call write_log ("dalec will be used to determine ecosystem phenology", msg_info, __FILE__ , __LINE__ )
      else
        call write_log ("ecosystem phenology will be prescribed", msg_info, __FILE__, __LINE__ )
      end if
      call GetValue(section,'dead_lai_energy_balance'  ,     user_opts%dead_lai_energy_balance )
      if (user_opts%dead_lai_energy_balance) then
        call write_log ("Including dead lai in energy balance", msg_info, __FILE__ , __LINE__ )
      else
        call write_log ("Not including dead lai in energy balance", msg_info, __FILE__, __LINE__ )
      end if
      call GetValue(section,'iterate_canopy'  ,     user_opts%iterative_canopy )
      if (user_opts%iterative_canopy) then
        call write_log ("Using iterative canopy option", msg_info, __FILE__ , __LINE__ )
      else
        call write_log ("Non-iterative canopy option in use", msg_info, __FILE__, __LINE__ )
      end if
      call GetValue(section, 'solve_canopy_temperature', user_opts%solve_canopy_temperature)
      if (user_opts%solve_canopy_temperature == 1) then
        call write_log ("Using steady state solution to canopy temperature", msg_info, __FILE__ , __LINE__ )
      else if (user_opts%solve_canopy_temperature == 2) then
        call write_log ("Using iterative solution to canopy temperature", msg_info, __FILE__ , __LINE__ )
      else
        call write_log ("Invalid option for solving canopy temperature", msg_fatal, __FILE__ , __LINE__ )
      end if
      call GetValue(section,'start_from_restart',user_opts%load_restart     )
      call GetValue(section,'plant_func_type',user_opts%plant_func_type     )
      select case (user_opts%plant_func_type)
      case (0) 
        call write_log( "Evergreen simulation." , msg_info , __FILE__ , __LINE__ )
        ! assign carbon model dimensions
        nopars = nopars_default ; nopools = nopools_default ; nofluxes = nofluxes_default
        allocate(fluxes_names(nofluxes),pools_names(nopools),pars_names(nopars))
        fluxes_names = fluxes_names_default ; pools_names = pools_names_default ; pars_names = pars_names_default
      case (1)
        call write_log( "Deciduous simulation." , msg_info , __FILE__ , __LINE__ )
        ! assign carbon model dimensions
        nopars = nopars_default ; nopools = nopools_default ; nofluxes = nofluxes_default
        allocate(fluxes_names(nofluxes),pools_names(nopools),pars_names(nopars))
        fluxes_names = fluxes_names_default ; pools_names = pools_names_default ; pars_names = pars_names_default
      case (2)
        call write_log( "Cereal crops/C4-grasses simulation." , msg_info , __FILE__ , __LINE__ )
        ! assign carbon model dimensions
        nopars = nopars_crops ; nopools = nopools_crops ; nofluxes = nofluxes_crops
        allocate(fluxes_names(nofluxes),pools_names(nopools),pars_names(nopars))
        fluxes_names = fluxes_names_crops ; pools_names = pools_names_crops ; pars_names = pars_names_crops
      case (3)
        call write_log( "Tuber crop." , msg_info , __FILE__ , __LINE__ )
        ! assign carbon model dimensions
        nopars = nopars_crops ; nopools = nopools_crops ; nofluxes = nofluxes_crops
        allocate(fluxes_names(nofluxes),pools_names(nopools),pars_names(nopars))
        fluxes_names = fluxes_names_crops ; pools_names = pools_names_crops ; pars_names = pars_names_crops
      case default
        write(message,*)"Plant func type not recognised - should be,"&
                      //" 1(evergreen/deciduous),  2(cereal crops)"&
                      //" or 3 (tuber crops)"
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      end select
      call GetValue(section,'layers_in_canopy',       grid%canopy         )
      ! sanity check..
      if ( grid%canopy .lt. 1 ) then
        call write_log( "Must be at least one canopy layer!" , msg_fatal , __FILE__ , __LINE__ )
      elseif ( grid%canopy .gt. 20 ) then
        call write_log( "Use of more than 20 canopy layers is not recommended!" , msg_warning , __FILE__ , __LINE__ )
      end if
      call GetValue(section,'layers_in_soil',         grid%soil           )
      ! sanity check..
      if ( grid%soil .lt. 3 ) then
        call write_log( "Must be at least 3 soil layers!" , msg_fatal , __FILE__ , __LINE__ )
      elseif ( grid%soil .gt. 40 ) then
        call write_log( "Use of more than 40 soil layers is not recommended!" , msg_warning , __FILE__ , __LINE__ )
      end if
      grid%core = grid%soil + 1
      call GetValue(section,'latitude_deg_north',     grid%latitude                  ) 
      ! convert latitude to module variable in radians
      latitude_radians=grid%latitude * deg_to_rad 
      call GetValue(section,'longitude_deg_east',     grid%longitude                 )
      call GetValue(section,'number_of_years',        time%nos_of_years              )
      call GetValue(section,'days_per_year',          dummy_integer                  )
      ! year cannot be less than 365 days if we want more than 1 year of
      ! simulation
      if (time%nos_of_years > 1 .and. dummy_integer < 365) then
         print*,"cannot have less than 365 days days_per_year if you want more than 1 year of simulation" 
         stop
      endif
      call GetValue(section,'leap_years',     user_opts%include_leap_years           )
      call GetValue(section,'start_date',          time%start_date                   )
      if (time%start_date == '') then
         print*,"start_date needed in format YYYY-MM-DD_HH:MM:SS"
         stop
      endif
      call GetValue(section,'steps_per_day',          time%steps_per_day             )
      ! now allocate the number of years in the simulation for the days in each
      ! of those years
      allocate(time%days_in_year(time%nos_of_years)) 
      if (.not.user_opts%include_leap_years) then
          time%days_in_year = dummy_integer 
      endif

      ! sanity check..
      if ( time%steps_per_day .lt. 24 ) then
        call write_log( "Use of less than 24 steps per day is not recommended (resolving the dirunal cycle is important!)" , &
                         msg_warning , __FILE__ , __LINE__ )
      end if
      time%seconds_per_step = ( 24d0 * 3600d0 ) / time%steps_per_day
      call GetValue(section,'loop_met_driver',       user_opts%loop_met_driver         )
      call GetValue(section,'loop_phenology_driver', user_opts%loop_pheno_driver       )
      call GetValue(section,'log_messages_to_screen:', user_opts%print_msg_at_severity )
      call set_msg_level( user_opts%print_msg_at_severity )
    else
      write(message,*)"Could not find [Options] section in the config file.  Cannot continue without it!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! Get the initialisation input files...
    call GetSection(section_names,section,'Initialisation')
    if ( associated(section) ) then
       if ( user_opts%load_restart ) then
         ! Restart model..
         call GetValue(section,'restart',user_opts%restart_filename )
         if (user_opts%restart_filename=='') then
            write(message,*)'Restart filename must be supplied in the user-config file.'
            call write_log( message , msg_fatal, __FILE__ , __LINE__ )
         end if
       else ! for restart
         ! Standard initialisation..
         ! soils
         call GetValue(section,'soils',user_opts%soils_filename)
         if (user_opts%soils_filename=='') then
            write(message,*)'Soils filename must be supplied in the user-config file.'
            call write_log( message , msg_fatal, __FILE__ , __LINE__ )
         end if
         ! vegetation
         call GetValue(section,'vegetation',user_opts%veg_filename)
         if (user_opts%veg_filename=='') then
            write(message,*)'Vegetation filename must be supplied in the user-config file.'
            call write_log( message , msg_fatal , __FILE__ , __LINE__ )
         end if
         ! carbon model
         call GetValue(section,'carbon',user_opts%carbon_filename)
         if (user_opts%carbon_filename=='' .and. user_opts%use_dalec) then
            write(message,*)'Carbon filename must be supplied in the user-config file.'
            call write_log( message , msg_fatal , __FILE__ , __LINE__ )
         end if
       end if ! end restart
       ! crops
       if (user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3) then
          call GetValue(section,'crops',user_opts%crops_filename)
          if (user_opts%crops_filename=='') then
             write(message,*) 'Crops filename must be supplied in the user-config file.'
             call write_log( message , msg_fatal , __FILE__ , __LINE__ )
          end if
       end if
       ! ecosystem phenology file if needed
       if (.not.user_opts%use_dalec) then
          call GetValue(section, 'phenology', user_opts%phenology_filename)
          if (user_opts%phenology_filename=='') then
             write(message,*) 'Phenology filename must be supplied in the user-config file.'
             call write_log( message , msg_fatal , __FILE__ , __LINE__ )
          end if
       end if 
    else ! for associated
       write(message,*)"Could not find [Initialisation] section in the config file."
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! Meteorology driver
    call GetSection( section_names , section , 'Meteorology driver' )
    if ( associated(section) ) then
      call GetValue( section , 'filename' , user_opts%met_filename )
      if (user_opts%met_filename=='') then
        write(message,*) "Meteorology filename must be supplied in the user-config file."
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      end if
      call GetValue( section , 'temp_in_kelvin' , user_opts%temp_in_kelvin )
      call GetValue( section , 'file_is_netcdf' , user_opts%met_file_is_nc )
      call GetValue( section , 'file_has_co2' ,  user_opts%use_co2_from_met_file  )
      if ( .not. user_opts%use_co2_from_met_file ) then 
        call write_log( "Using a constant CO2 concentration." , msg_info , __FILE__ , __LINE__ )
        call GetValue( section , 'constant_CO2_concentration' , met%const_ambient_co2 )
      end if
      call GetValue( section , 'file_has_par' ,  user_opts%use_ppfd_from_met_file  )
      if ( .not. user_opts%use_ppfd_from_met_file ) then
        call write_log( "PPFD will be calculated from SW-radiation." , msg_info , __FILE__ , __LINE__ )
      end if
      call GetValue(section,'file_has_incoming_longwave', user_opts%met_file_has_lw_rad )
      if ( .not. user_opts%met_file_has_lw_rad ) then
        call write_log( "Incoming lw_rad will be calculated from surface temperature." , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue( section , 'file_has_sfc_pressure' , user_opts%met_file_has_sfc_press )
      call GetValue( section , 'file_has_sw_diffuse' , user_opts%met_file_has_sw_diffuse )
      call GetValue( section,'  file_has_snowfall', user_opts%met_file_has_snowfall )
      if ( .not. user_opts%met_file_has_snowfall ) then
        call write_log( "Snowfall will be estimates from precipitation < 0oC." , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'file_has_disturbance', user_opts%met_file_has_disturbance)

#ifdef USE_NETCDF
      if ( user_opts%met_file_is_nc ) then
        allocate( met_nc%header )
        met_nc%header%name = user_opts%met_filename
        call get_met_nc_varnames_from_config( section , met_nc )
      end if
#endif
      call GetValue( section , 'precip_is_rate' , user_opts%precip_is_rate )
      call GetValue( section , 'vpd_or_specific_humidity', user_opts%vpd_or_specific_humidity)
      if (trim(user_opts%vpd_or_specific_humidity) /= "vpd" .and. & 
          trim(user_opts%vpd_or_specific_humidity) /= "specific_humidity") then
          write(message,*)'Invalid vpd_or_specific_humidity option provided, valid options are "vpd" and "specific_humidity"'
          print*,"you wrote ",trim(user_opts%vpd_or_specific_humidity)
          call write_log( message , msg_fatal, __FILE__ , __LINE__ )
      end if
    else
      write(message,*)"Could not find [Meteorology driver] section in the config file."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! Get the general output detail..
    call GetSection(section_names,section,'Output')
    if (associated(section)) then
      call GetValue( section , 'directory'               , user_opts%output_directory        )
      ! check the directory exists by trying to create a file in it..
      outdirtestfile = trim(user_opts%output_directory)//'spa_test_dir_exists'
      open( iostat=ios , unit=999 , file=outdirtestfile , status='new' )
      close( 999 , status='delete' )
      if ( ios .ne. 0 ) then
        write(message,*)"Cannot find output directory: ",trim(user_opts%output_directory)
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      end if
      call GetValue( section , 'standard_csv_output'     , user_opts%std_csv_output          )
      call GetValue( section , 'write_restart_file'      , user_opts%make_restart_file       )
      call GetValue( section , 'restart_write_frequency' , user_opts%restart_write_frequency )
      call GetValue( section , 'Ecofluxes_csv_output'    , user_opts%Ecofluxes_csv_output)
      call GetValue( section , 'canopy_csv_output'       , user_opts%canopy_csv_output)
      call GetValue( section , 'Cstock_csv_output'       , user_opts%Cstock_csv_output)
      call GetValue( section , 'soils_csv_output'        , user_opts%soils_csv_output)

    end if

    call allocate_arrays

  end subroutine read_user_config
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine update_parameters_from_user_config( section_names )

    ! make sure any parameters provided by the user overwrite !
    ! the values taken from the input files.                  !

    use config_tools
    use carbon_model_crop_mod, only: harvest_day, plough_day, sow_day, vernal_calcs
    use gv_scale_declarations, only: met
    use gv_snow_info,          only: snowheight, snowweight

    implicit none

    ! arguments..
    type(ConfigSection), pointer :: section_names ! holds structure of the configuration file   

    ! local variables..
    type(ConfigSection), pointer :: section

    ! For each section, we look for possible parameters.
    ! GetValue only changes the third argument IF it
    ! finds something.

    call GetSection( section_names , section , 'Meteorology Parameters' )
    if (associated(section)) then
       call GetValue(section,'constant_CO2_concentration',met%const_ambient_co2)
       call GetValue(section,'constant_sfc_pressure',met%const_sfc_pressure)
    end if
 
    call GetSection(section_names,section,'Crop Parameters')
    if (associated(section)) then
       call GetValue(section,'plough_day_number' ,plough_day  )
       call GetValue(section,'sow_day_number'    ,sow_day     )
       call GetValue(section,'harvest_day_number',harvest_day )
       call GetValue(section,'vernalization'     ,vernal_calcs)
    end if

    call GetSection(section_names,section,'Soil Parameters')
    if (associated(section)) then
       call GetValue(section,'snowweight',snowweight)
       call GetValue(section,'snowheight',snowheight)
    end if

  end subroutine update_parameters_from_user_config
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine allocate_arrays

    ! The size of some arrays depends on grid-size, which !
    ! can be set at runtime. These arrays must then be    !
    ! allocated at the start of a run, rather than being  !
    ! hard-coded. This routine allocates them.            !

    use gv_clim,               only: gbw, wetev, leaf_heat_conductance
    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid, time
    use gv_soil_structure
    use gv_veg
    use gv_snow_info
    use gv_carbon_model,       only: FLUXES, POOLS, pars, &
                                     nopars, nofluxes, nopools

    implicit none

    ! clim..
    call check_allocated( gbw   , grid%canopy )
    call check_allocated( wetev , time%steps_per_day )
    call check_allocated( leaf_heat_conductance , grid%canopy )

    ! hydrol..
    call check_allocated( soil_frac_clay , grid%core )
    call check_allocated( soil_frac_sand , grid%core )

    ! snow..
    if ( .not. allocated(            Dsfix ) )  allocate(            Dsfix( grid%snow ) )
    if ( .not. allocated(            Dsnow ) )  allocate(            Dsnow( grid%snow ) )
    if ( .not. allocated(             Sice ) )  allocate(             Sice( grid%snow ) )
    if ( .not. allocated(             Sliq ) )  allocate(             Sliq( grid%snow ) )
    if ( .not. allocated(            Tsnow ) )  allocate(            Tsnow( grid%snow ) )

    ! soil..
    if ( .not. allocated(cond1) ) then
      ! (it's fair to assume that if one hasn't been, none of them have)..
      allocate( cond1(grid%core), cond2(grid%core), cond3(grid%core), potA(grid%core), potB(grid%core) )
    end if
    if ( .not. allocated(           conduc ) )  allocate(           conduc( grid%core ) )
    if ( .not. allocated(   field_capacity ) )  allocate(   field_capacity( grid%core ) )
    if ( .not. allocated(  fraction_uptake ) )  allocate(  fraction_uptake( grid%core ) )
    if ( .not. allocated(          iceprop ) )  allocate(          iceprop( grid%core ) )
    if ( .not. allocated(      layer_depth ) )  allocate(      layer_depth( grid%core ) )
    if ( .not. allocated(      mineralfrac ) )  allocate(      mineralfrac( grid%core ) )
    if ( .not. allocated(      organicfrac ) )  allocate(      organicfrac( grid%core ) )
    if ( .not. allocated(         porosity ) )  allocate(         porosity( grid%core ) )
    if ( .not. allocated(          pptgain ) )  allocate(          pptgain( grid%core ) )
    if ( .not. allocated(      root_length ) )  allocate(      root_length( grid%core ) )
    if ( .not. allocated(        root_mass ) )  allocate(        root_mass( grid%core ) )
    if ( .not. allocated(        soil_temp ) )  allocate(        soil_temp( grid%core ) )
    if ( .not. allocated( soil_temp_nplus1 ) )  allocate( soil_temp_nplus1( grid%core ) )
    if ( .not. allocated(            soilR ) )  allocate(            soilR( grid%core ) )
    if ( .not. allocated(           soilR1 ) )  allocate(           soilR1( grid%core ) )
    if ( .not. allocated(           soilR2 ) )  allocate(           soilR2( grid%core ) )
    if ( .not. allocated(              SWP ) )  allocate(              SWP( grid%core ) )
    if ( .not. allocated(        thickness ) )  allocate(        thickness( grid%core ) )
    if ( .not. allocated(        waterfrac ) )  allocate(        waterfrac( grid%core ) )
    if ( .not. allocated(        watergain ) )  allocate(        watergain( grid%core ) )
    if ( .not. allocated(       watericemm ) )  allocate(       watericemm( grid%soil ) )
    if ( .not. allocated(        waterloss ) )  allocate(        waterloss( grid%core ) )
    if ( .not. allocated( wettingbot ) )  then
      allocate( wettingbot( grid%wetting ) ) 
      allocate( wettingtop( grid%wetting ) )
    end if

    ! veg..
    if ( .not. allocated( c3 ) ) allocate( c3( grid%canopy ) ) ; c3 = .False.
    if ( .not. allocated( canopy_soil_resistance ) )  allocate( canopy_soil_resistance( grid%canopy )        )
    if ( .not. allocated(                    ess ) )  allocate( ess( time%steps_per_day )                    )
    if ( .not. allocated(                   gppt ) )  allocate( gppt( time%steps_per_day )                   )
    if ( .not. allocated(  system_energy_balance ) )  allocate( system_energy_balance ( time%steps_per_day ) )
    if ( .not. allocated(                 lafrac ) )  allocate( lafrac( grid%canopy )                        )
    if ( .not. allocated(            lafrac_dead ) )  allocate( lafrac_dead( grid%canopy )                   )
    if ( .not. allocated(                    lai ) )  allocate( lai( grid%canopy )                           )
    if ( .not. allocated(                psil_pd ) )  allocate( psil_pd( grid%canopy )                       )
    if ( .not. allocated(               dead_lai ) )  allocate( dead_lai( grid%canopy )                      )
!    if ( .not. allocated(           Cfol_profile ) )  allocate( Cfol_profile( grid%canopy )                  )
    if ( .not. allocated(                    LCA ) )  allocate( LCA( grid%canopy )                           )
    if ( .not. allocated(     can_sensible_layer ) )  allocate( can_sensible_layer( grid%canopy )            )
    if ( .not. allocated( leaf_temp_intermediate ) )  allocate( leaf_temp_intermediate( grid%canopy )        )
    if ( .not. allocated(              leaf_temp ) )  allocate( leaf_temp( grid%canopy )                     )
    if ( .not. allocated(           layer_height ) )  allocate( layer_height( grid%core )                    )
    if ( .not. allocated(               LWPstore ) )  allocate( LWPstore( grid%canopy )                      )
    if ( .not. allocated(            LWPprevious ) )  allocate( LWPprevious( grid%canopy )                   )
    if ( .not. allocated(                  nfrac ) )  allocate( nfrac( grid%canopy )                         )
    if ( .not. allocated(                    nla ) )  allocate( nla( grid%canopy )                           )
    if ( .not. allocated(                  respt ) )  allocate( respt( time%steps_per_day )                  )
    if ( .not. allocated(               soiletmm ) )  allocate( soiletmm( time%steps_per_day )               )
    if ( .not. allocated(            transt_kg_s ) )  allocate( transt_kg_s( time%steps_per_day )            )
    if ( .not. allocated(                 transt ) )  allocate( transt( time%steps_per_day )                 )
    if ( .not. allocated(          sensible_heat ) )  allocate( sensible_heat(time%steps_per_day)            )

    ! carbon...
    if ( .not. allocated(            FLUXES ) )  allocate( FLUXES(time%steps_per_day,nofluxes)  )
    if ( .not. allocated(             POOLS ) )  allocate( POOLS(1,nopools)                     )
    if ( .not. allocated(              pars ) )  allocate( pars(nopars)                         )

    ! time...
    if ( .not. allocated(    lemod ) ) allocate(    lemod(time%steps_per_day) )
    if ( .not. allocated(    mmmod ) ) allocate(    mmmod(time%steps_per_day) )
    if ( .not. allocated(   neemod ) ) allocate(   neemod(time%steps_per_day) )
    if ( .not. allocated( timeflux ) ) allocate( timeflux(time%steps_per_day) )
    if ( .not. allocated(    wetle ) ) allocate(    wetle(time%steps_per_day) )
    if ( .not. allocated(potential_evaporation) ) allocate(potential_evaporation(time%steps_per_day))
    if ( .not. allocated(netrad_day) ) allocate(netrad_day(time%steps_per_day) )

    ! set some variables as zero at initialisation
    potA = 0d0 ; potB = 0d0 ; cond1 = 0d0 ; cond2 = 0d0 ; cond3 = 0d0
    leaf_temp=0d0 ; lai = 0d0 ; psil_pd = 0d0 ; dead_lai = 0d0 ; lafrac = 0d0 ; lafrac_dead = 0d0 ; LCA = 0d0
    system_energy_balance = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; pars = 0d0 !; Cfol_profile = 0d0
    call check_allocated( flux , varsize=(/time%steps_per_day,grid%canopy/) )

  end subroutine allocate_arrays
  !
  !----------------------------------------------------------------------
  !
  subroutine check_allocated_1d( var, varsize )

    use log_tools

    implicit none

    ! arguments..
    double precision,allocatable,intent(inout) :: var(:)
    integer,intent(in) :: varsize

    if ( .not. allocated( var ) ) then
        ! allocate it..
        allocate( var( varsize ) )
    else
        ! check it is actually the desired size..
         if ( size(var) .ne. varsize ) then
             write(message,*)"Variable not allocated to same size as stated varsize!"
             call write_log( message , msg_error , __FILE__ , __LINE__ )
         end if
     end if

  end subroutine
  !
  !----------------------------------------------------------------------
  !
  subroutine check_allocated_2d( var, varsize )

    use log_tools

    implicit none

    ! arguments..
    double precision,allocatable,intent(inout)  :: var(:,:)
    integer,dimension(:),intent(in) :: varsize

    ! local variables..
    integer :: varshape(2)

    if ( .not. allocated( var ) ) then
        ! allocate it..
        allocate( var( varsize(1),varsize(2) ) )
    else
        ! check it is actually the desired size..
         varshape = shape(var)
         if ( ( varshape(1) .ne. varsize(1) ) .or. ( varshape(2) .ne. varsize(2) ) ) then
             write(message,*)"Variable not allocated to same size as stated varsize!"
             call write_log( message , msg_error , __FILE__ , __LINE__ )
         end if
     end if

  end subroutine check_allocated_2d
  !
  !----------------------------------------------------------------------
  !
#ifdef USE_NETCDF
  subroutine get_met_nc_varnames_from_config( section , met_nc )

    ! If the input met file is netcdf we need to check the !
    !  [Meteorology driver] section for details of what    !
    !  the variables are called in the met input file.     !

    use config_tools
    use gv_scale_declarations, only: met, user_opts
    use log_tools
    use spa_io_netcdf,         only: nc_met_file

    implicit none

    ! arguments..
    type(ConfigSection),pointer :: section
    type(nc_met_file),pointer   :: met_nc

    ! define default names for dimensions, and then go see if the user
    ! has given different ones...
    ! (I'm sure there should be a better way of doing this...)
    allocate( met_nc%time   ) ; met_nc%time%name = 'time'
    allocate( met_nc%lat    ) ; met_nc%lat%name  = 'lat'
    allocate( met_nc%lon    ) ; met_nc%lon%name  = 'lon'
    call GetValue( section, 'time',      met_nc%time%name   )
    call GetValue( section, 'latitude',  met_nc%lat%name    )
    call GetValue( section, 'longitude', met_nc%lon%name    )

    ! define default names for variables.. and then go see if the user
    ! has given different ones...
    ! (I'm sure there should be a better way of doing this...)
    allocate( met_nc%co2    ) ; met_nc%co2%name    = 'co2'
    allocate( met_nc%par    ) ; met_nc%par%name    = 'par'
    allocate( met_nc%precip    ) ; met_nc%precip%name    = 'precip'
    allocate( met_nc%airt    ) ; met_nc%airt%name    = 'temp'
    allocate( met_nc%sfc_pressure  ) ; met_nc%sfc_pressure%name  = 'sfc_pressure'
    allocate( met_nc%sw_rad  ) ; met_nc%sw_rad%name  = 'sw_rad'
    allocate( met_nc%sw_diffuse  ) ; met_nc%sw_diffuse%name  = 'sw_diffuse'
    allocate( met_nc%lw_rad  ) ; met_nc%lw_rad%name  = 'lw_rad'
    allocate( met_nc%vpd    ) ; met_nc%vpd%name    = 'vpd'
    allocate( met_nc%sp_moist) ; met_nc%sp_moist%name = 'sp_moist'
    allocate( met_nc%wind_spd ) ; met_nc%wind_spd%name = 'wind_spd'

    call GetValue( section, 'air_temperature',         met_nc%airt%name      )
    call GetValue( section, 'temp_in_kelvin',          met_nc%airt_in_kelvin )
    if ( user_opts%use_co2_from_met_file ) then
      call GetValue( section, 'carbon_dioxide',        met_nc%co2%name      )
    else
      met_nc%co2%name = ''
    end if
    if ( user_opts%use_ppfd_from_met_file ) then
      call GetValue( section, 'photo_active_rad',       met_nc%par%name     )
    else
      met_nc%par%name = ''
    end if
    call GetValue( section, 'precipitation',           met_nc%precip%name      )
    if ( user_opts%met_file_has_sfc_press ) then
      call GetValue( section, 'surface_pressure',      met_nc%sfc_pressure%name    )
    else
      met_nc%sfc_pressure%name = ''
      call GetValue(section,'constant_sfc_pressure',   met%const_sfc_pressure )
    end if
    call GetValue( section, 'sw_radiation',            met_nc%sw_rad%name    )
    call GetValue( section, 'vapour_pressure_deficit', met_nc%vpd%name      )
    call GetValue( section, 'wind_speed',              met_nc%wind_spd%name   )

    call write_log("The variable names SPA is expecting to find in met nc file are:", &
                     msg_info, __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%time%name),' ',trim(met_nc%lat%name),' ',trim(met_nc%lon%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%co2%name),' ',trim(met_nc%par%name),' ',  &
                     trim(met_nc%precip%name),' ',trim(met_nc%airt%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%sfc_pressure%name),' ',trim(met_nc%sw_rad%name),' ', &
                     trim(met_nc%vpd%name),' ',trim(met_nc%wind_spd%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )

  end subroutine get_met_nc_varnames_from_config
#endif
  !
  !----------------------------------------------------------------------
  !
end module spa_config
!
!----------------------------------------------------------------------
!

