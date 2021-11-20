! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module spa_io_netcdf

  !! This module declares derived types describing SPA-specific  !!
  !! input/output files, and procedures for reading those files. !!

  use gv_scale_declarations, only: fname_length
  use netcdf_tools,          only: nc_header, nc_dimension, nc_variable

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: load_met_nc_data, open_met_nc
  ! Variables..
  public :: gridcalc, nc_met_file


  ! Required for finding particular spot on grid that we will use..
  type gridcalc
    integer,dimension(2) :: indices  = (1,1)   ! bottom-left corner of box of grid-pts
    ! that encompass user-desired location.
    logical           :: no_bilinear = .false. ! bilinear flag (assume false=>do calcs)
  end type gridcalc

  ! Input Meteorolgy file structure..
  type nc_met_file
    type(nc_header),pointer    :: header       => null() ! information about the file, such as name & dims
    type(nc_dimension),pointer :: time         => null() ! time
    type(nc_dimension),pointer :: lat          => null() ! latitude
    type(nc_dimension),pointer :: lon          => null() ! longitude
    type(nc_variable),pointer  :: airt         => null() ! surface air temperature
    logical                    :: airt_in_kelvin = .False. ! assume SAT in Celcius
    type(nc_variable),pointer  :: co2          => null() ! carbon dioxide atmospheric concentration
    type(nc_variable),pointer  :: par          => null() ! photosyntheticaly active radiation
    type(nc_variable),pointer  :: precip       => null() ! precipitation
    type(nc_variable),pointer  :: sw_rad       => null() ! shortwave radiation
    type(nc_variable),pointer  :: sw_diffuse   => null() ! diffuse short wave
    type(nc_variable),pointer  :: sfc_pressure => null() ! surface atmospheric pressure
    type(nc_variable),pointer  :: vpd          => null() ! vapour pressure deficit
    type(nc_variable),pointer  :: sp_moist     => null() ! specific moisture
    type(nc_variable),pointer  :: lw_rad       => null() ! long wave radiation
    type(nc_variable),pointer  :: wind_spd     => null() ! wind speed
!    type(nc_variable),pointer  :: snowfall => null() ! snow fall
    type(gridcalc),pointer     :: grid         => null() ! point on grid that we will use
  end type nc_met_file

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine load_met_nc_data( ncfile )

    ! Read the met input file and load data to pointers !
    ! Unlike the csv file, which is read line-by-line as!
    ! it is needed, the netcdf file contents are read   !
    ! all in one go.                                    !

    use gv_scale_declarations, only: met, time, user_opts
    use log_tools
    use netcdf_tools,          only: get_nc_var

    implicit none

    ! arguments..
    type(nc_met_file)  :: ncfile

    ! local variables..
    integer :: time_length

    time_length = size(ncfile%time%values)

    ! Carbon dioxide..
    allocate( met%ambient_co2(time_length) )
    if ( user_opts%use_co2_from_met_file ) then   
      call write_log("Try to load Ambient CO2 concentration from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%co2 )
      met%ambient_co2 = reshape( ncfile%co2%data%grid_3d , (/ time_length /) )
      call write_log("..successful.")
    else
      ! use a default value..
      met%ambient_co2 = met%const_ambient_co2
      call write_log("Ambient CO2 concentration set to default value (330ppm)")
    end if

    ! Precipitation
    allocate( met%precip(time_length) )
    call write_log("Try to load Precipitation from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%precip )
    call write_log ("..successful.")
    met%precip = reshape( ncfile%precip%data%grid_3d , (/ time_length /) )
    if ( user_opts%precip_is_rate ) then
      call write_log("Converting precipitation from rate to volume by "&
                   //"multiplying by nos of seconds per timestep.")
      met%precip = met%precip * time%seconds_per_step
    end if

    ! Short wave radiation
    allocate( met%sw_rad(time_length) )
    call write_log("Try to load SW radiation from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%sw_rad )
    call write_log("..successful.")
    met%sw_rad = reshape( ncfile%sw_rad%data%grid_3d  , (/ time_length /) )

    if ( user_opts%met_file_has_sw_diffuse ) then
       ! Diffuse short wave radiation
       allocate( met%sw_diffuse(time_length) )
       call write_log("Try to load diffuse SW radiation from met netcdf file..")
       call get_nc_var( ncfile%header%id, ncfile%sw_diffuse )
       call write_log("..successful.")
       met%sw_diffuse = reshape( ncfile%sw_diffuse%data%grid_3d  , (/ time_length /) )
    endif

    allocate( met%lw_rad(time_length) )
    if ( user_opts%met_file_has_lw_rad ) then
       ! Long wave radiation
       call write_log("Try to load LW radiation from met netcdf file..")
       call get_nc_var( ncfile%header%id, ncfile%lw_rad )
       call write_log("..successful.")
       met%lw_rad = reshape( ncfile%lw_rad%data%grid_3d  , (/ time_length /) )
    endif

    ! Surface air temperature
    allocate( met%temp(time_length) )
    call write_log("Try to load Surface air temperature from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%airt )
    ! Adjust to Kelvin, if needed, and also check that units seem sane..
    if ( ncfile%airt_in_kelvin ) then
      if ( minval( ncfile%airt%data%grid_3d ) .lt. 100. ) then
        write(message,*)"Units of temperature appear to be wrong. You "//&
             "specified Kelvin, but minimum is below 100K; is this correct??"
        call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
      end if
      ncfile%airt%data%grid_3d = ncfile%airt%data%grid_3d - 273.15 ! convert T to Celcius
    else
      if ( maxval( ncfile%airt%data%grid_3d ) .gt. 150. ) then
        write(message,*)"Units of temperature appear to be wrong. You "//&
             "specified Celcius, but maximum is above 150C; is this correct??"
        call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
      end if
    end if
    call write_log("..successful.")
    met%temp = reshape( ncfile%airt%data%grid_3d , (/ time_length /) )

    ! Surface pressure
    allocate( met%sfc_pressure(time_length) )
    if ( user_opts%met_file_has_sfc_press ) then
      call write_log("Try to load Surface pressure from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%sfc_pressure )
      call write_log("..successful.")
      met%sfc_pressure = reshape( ncfile%sfc_pressure%data%grid_3d , (/ time_length /) )
    else
      ! use a default value..
      met%sfc_pressure = met%const_sfc_pressure
      call write_log("Surface pressure set to default value (100000 Pa)")
    end if

    ! Vapour pressure deficit / specific humidity
    allocate( met%vpd(time_length) )
    if (user_opts%vpd_or_specific_humidity == "specific_humidity") then
       call write_log("Try to load specific humidity from met netcdf file..")
       call get_nc_var( ncfile%header%id, ncfile%vpd )
       call write_log("..successful.")
       met%vpd = reshape( ncfile%vpd%data%grid_3d , (/ time_length /) )
    else if (user_opts%vpd_or_specific_humidity == "vpd") then
       call write_log("Try to load vapour pressure deficit from met netcdf file..")
       call get_nc_var( ncfile%header%id, ncfile%vpd )
       call write_log("..successful.")
       met%vpd = reshape( ncfile%sp_moist%data%grid_3d , (/ time_length /) )
    else ! not specified
       write(message,*)"Have not specified vpd or specfic humidity"
       call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Wind speed
    allocate( met%wind_spd(time_length) )
    call write_log("Try to load wind speed from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%wind_spd )
    call write_log("..successful.")
    met%wind_spd = reshape( ncfile%wind_spd%data%grid_3d , (/ time_length /) )
    ! adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2 )  met%wind_spd = 0.2

    ! Photosynthetically active radiation
    allocate( met%ppfd(time_length) )
    if ( user_opts%use_ppfd_from_met_file ) then
      call write_log("Try to load photosynthetically active radiation from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%par )
      call write_log("..successful.")
      met%ppfd = reshape( ncfile%par%data%grid_3d , (/ time_length /) )
    else
      ! PAR isn't available in most climate models, so
      ! instead we use a loose approximation...
      met%ppfd = 2.3 * met%sw_rad
    end if

    

  end subroutine load_met_nc_data
  !
  !----------------------------------------------------------------------
  !
  subroutine open_met_nc( header , time , lat , lon )

    ! Open a netcdf input file, and fill out the header with !
    ! information on handles to the file, and the dimensions !
    ! of time, latitude and longitude.                       !

    use netcdf_tools, only: get_dim_info, nc_header, nc_dimension, open_nc_file

    implicit none

    ! arguments..
    type(nc_header),pointer    :: header
    type(nc_dimension),pointer :: time, lat, lon

    ! open netCDF file..
    call open_nc_file( header )

    ! populate each dimension with their dim & var handle
    !  ids, their length and their actual values..
    call get_dim_info( header%id , time%name , time%dim_id , &
                                      time%var_id , time%values )
    call get_dim_info( header%id, lat%name , lat%dim_id , &
                                      lat%var_id , lat%values )
    call get_dim_info( header%id, lon%name , lon%dim_id , &
                                      lon%var_id , lon%values  )

  end subroutine open_met_nc
  !
  !----------------------------------------------------------------------
  !
end module spa_io_netcdf
!
!------------------------------------------------------------------------
!
