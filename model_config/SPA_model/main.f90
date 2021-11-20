! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


program main_spa

  !! Soil Plant Atmosphere model coupled to DALEC !!

  use canopy,                  only: timestep_calcs
  use gv_scale_declarations,   only: time, user_opts
  use gv_soil_structure,       only: rooted_layers, root_mass, root_length, &
                                     potB, watericemm
  use gv_veg,                  only: stock_roots
  use soil_air,                only: roots
  use log_tools
  use spa_io,                  only: handle_output, start_spa, & 
                                     update_met_drivers, update_phenology_drivers

  implicit none

  ! local variables
  integer            :: yr, day, step, start_day
  double precision :: start, finish

  ! get some time information
  call cpu_time(start)

  ! start logging..
  call open_log( unit=100 , fname="spa.log" )

  ! read user config, open files & initialise (spa_io.f90)
  call start_spa
 
  ! starting point and adjust the time%day back one day. This is cause time% day
  ! will be incremented to the correct day again in "increment_day" 
  start_day = time%day ; time%day = time%day - 1

  ! for each year..
  do yr = 1 , time%nos_of_years

     ! in first year don't do this to avoid resetting user inputted start day
     if (yr == 1) then
        call increment_time( 'first_year' , time , start_day)
     else
        call increment_time( 'year' , time , start_day)
     end if

     ! for each day...
     do day = start_day , time%days_in_year(yr)

        ! increment day component of model timing information
        call increment_time( 'day' , time, start_day)

        if (.not. user_opts%use_dalec) then
           ! update phenology for lai and root C
           call update_phenology_drivers
        end if ! use_dalec

        ! Update root structure..
        call write_log( 'Updating root structure.' , msg_info , __FILE__ , __LINE__ )
        call roots(stock_roots,rooted_layers,root_length,root_mass)

        ! for each sub-daily slice...
        do step = 1 , time%steps_per_day
if (maxval(watericemm) /= maxval(watericemm)) then
   print*,"error", step, day, yr
   print*,"water",watericemm
   stop
end if
           call increment_time( 'step' , time , start_day)

           ! get met driver for step (spa_io.f90)
!           call write_log( 'Get next chunk of met data' , msg_info , __FILE__ , __LINE__ )
           call update_met_drivers( time )

           ! transfer carbon and water (canopy.f90)
!           call write_log('Entering the timestep calculations' , msg_info , __FILE__ , __LINE__ )
           call timestep_calcs( time )

           ! write output if needed (spa_io.f90)
!           call write_log('Dealing with any output (if needed)' , msg_info , __FILE__ , __LINE__ )
           if (user_opts%Ecofluxes_csv_output) call handle_output( 1 , time )
           if (user_opts%Cstock_csv_output) call handle_output(2, time)

           call write_log_div

        enddo ! step of day

        call write_log('              ' , msg_info , __FILE__ , __LINE__ )
        call write_log_div

     enddo ! day of year

     ! End of year, sync the log..
     call sync_log

#ifdef USE_NETCDF
     ! check if it is time to write a restart file...
     call check_restart
#endif

    ! write the crops output..
    if (user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3) call handle_output( 8 , time )

  enddo ! year of run

#ifdef USE_NETCDF
  ! Check that, if restarts were wanted, we wrote at least one!..
  call check_restart( end_of_run = .true. )
#endif

  ! finish logging..
  call close_log( 0 )

  ! final shout of the model
  call cpu_time(finish)
  print '("   SPA time to completation = ",f9.5," seconds.")',finish-start

contains
  !
  !-------------------------------------------
  !
  subroutine increment_time( type , time , start_day)

    ! Increment parts of the time-holder variable. !
  
    use gv_scale_declarations, only: time_holder
    use log_tools,             only: write_log

    implicit none

    ! arguments..
    character(*),intent(in)         :: type
    type(time_holder),intent(inout) :: time
    integer, intent(inout) :: start_day

    select case (type)
    case ('first_year')
      time%year = time%year + 1
      time%step = 0
      write(message,*)'year is: ',time%year
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('year')
      ! next year, reset day & step..
      time%year = time%year + 1
      time%day  = 0
      time%step = 0
      start_day = 1 ! reset counter
      write(message,*)'year is: ',time%year
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('day')
      ! next day, reset step...
      time%day  = time%day + 1
      time%run_day = time%run_day + 1
      time%step = 0
      write(message,*)'day is: ',time%day
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('step')
      ! next step.
      ! Also update the count of total-steps, and the daytime..
      time%step = time%step + 1
      write(message,*)'step-of-day is: ',time%step
      call write_log( trim(message), msg_info  , __FILE__ , __LINE__ )
      time%steps_count = time%steps_count + 1
      write(message,*)'number of steps so far is: ',time%steps_count
      call write_log( trim(message), msg_info , __FILE__ , __LINE__ )

      ! Re-calculate the current time in units of days..
      ! Consider step 1 to be at 00:00 on the day, and the
      ! last step of the day to be at or before 23.59...
      time%daytime = ( time%year -1 ) * time%days_in_year(time%year)    &
                      + time%day                                  &
                       + ( dble(time%step) - 1 ) / dble(time%steps_per_day)
      write(message,*)'daytime is: ',time%daytime
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case default
      write(message,*)'Do not recognise type of update: ',type
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
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
end program main_spa

