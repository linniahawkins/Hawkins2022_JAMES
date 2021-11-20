! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module light

  !!> description needed <!!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: solar, solar_update_for_iterative_canopy, &
            solar_update_for_marginal_return
  ! Variables..
  public :: checkpar, leaf_sun, longem, nirabs, parabs_shade, parabs_sun &
           ,nirrefl_crop,parrefl_crop,soilpar,soilpar_coef,soilnir,soilnir_coef

  double precision :: fdiff,              & ! ?
          declination,        & !
          sky_absorbed_ppfd,  & ! ?
          soilabsn,           & ! ?
          soil_absorbed_ppfd, & ! ?
          soilabsl,           & ! ?
          soilnir_coef,       & ! surface litter adjustment to soilnir
          soilpar_coef          ! surface litter adjustment to soilpar

  double precision,dimension(:),allocatable :: &
         checkpar,     & ! check PAR radiation balances
         leaf_sun,     & ! fraction of leaf in a given canopy layer exposed to sun
         longem,       & ! long wave emission from each canopy layer
         nirabs,       & ! NIR abosrbed by each canopy layer
         parabs_shade, & ! PAR absorbed by each canopy layer in sun exposed leaves
         parabs_sun      ! PAR absorbed by each canopy layer in shade covered leaves

  double precision,parameter :: So = 1360d0,  & ! solar constant (Wm-2)
                           nirrefl = 0.43d0,  & ! NIR leaf reflectance (Sitka Spruce 0.16, grass/crop ~ 0.38)
                          nirtrans = 0.26d0,  & ! NIR transmittance from leaf layer to next
                           parrefl = 0.16d0,  & ! Baldocchi 0.11 ! Sitka Spruce 0.07, Grass / crop 0.11
                          partrans = 0.16d0,  & ! ?
                           soilpar = 0.033d0, & ! soil reflectance for PAR
                           soilnir = 0.023d0, & ! soil reflectance for NIR
                      nirrefl_crop = 0.50d0,  & ! NIR reflectance for dead crop (Nagler et al., 2003)
                      parrefl_crop = 0.30d0,  & ! PAR reflectance for dead crop (Nagler et al., 2003)
                         spherical = 2d0        ! spherical leaf angle distribution has a value 2.

  save

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine solar

    use gv_clim,                only: ppfd_to_par, par_top, sw_rad, sw_diffuse, lw_rad &
                                     ,temp_bot, temp_top 
    use gv_irradiance_sunshade, only: check, soilnet, skyabsl, daylength
    use gv_scale_declarations,  only: boltz, grid, pi, time, user_opts, deg_to_rad
    use gv_snow_info,           only: snowalb_nir, snowalb_par, snowheight, fsnow
    use gv_veg,                 only: lai, dead_lai, latitude_radians, modrnet, lafrac,lafrac_dead 
    use math_tools,             only: calculate_daylength_hours
    use spa_io,                 only: handle_output
    use log_tools

    implicit none

    ! local variables..
    integer :: counter, i
!    double precision :: albedo = 0.0, sw_absorbed = 0.0, sw_total = 0.0, fdiff_total = 0.0, kbm_total = 0.0
    double precision    :: ab, beam, clump, diff, em, estfdiff, ff, kbm, lnet, long,  &
               nirbeam, nirdiff, period, PIb, PId, radnet, rtime, skyabsn, soem,     &
               soilt, suml, sumn, sump, sun, SWR, totlong, totnir, totpar, xhour,    &
               soilnir_corrected, soilpar_corrected,      &
               reflsurf

    double precision,dimension(grid%canopy) :: abslong, abspar_shade, abspar_sun, absnir, sunfrac, temp, totparabs, lai_in_rad, &
                                   parrefl_vector, nirrefl_vector

    ! include an extra value each side bracketing grid%canopy..
    double precision,dimension(0:grid%canopy+1) :: downlong, downnir, downpar, uplong, upnir, uppar

!if (time%steps_count == 1) open( unit=999 , file='./bob_evap.txt' , status='replace' ,action='readwrite' )
    ! Check stuff is allocated before we get started..
    if (.not. allocated(    checkpar) )  allocate(    checkpar(grid%canopy) )
    if (.not. allocated(    leaf_sun) )  allocate(    leaf_sun(grid%canopy) )
    if (.not. allocated(      longem) )  allocate(      longem(grid%canopy) )
    if (.not. allocated(      nirabs) )  allocate(      nirabs(grid%canopy) )
    if (.not. allocated(parabs_shade) )  allocate(parabs_shade(grid%canopy) )
    if (.not. allocated(  parabs_sun) )  allocate(  parabs_sun(grid%canopy) )

    ! calculations..
    period    = dble(time%steps_per_day)
    clump     = 1d0    ! clumping factor - =1 means unclumped, random distribution
    ff        = 1d0   ! correction factor for stems etc on attenuation, currently unused and untested, set at 1.0
    totparabs = 0d0
    rtime     = dble( time%step )

    ! calculate declination 
    declination  = - asin ( sin ( 23.45d0 * deg_to_rad ) &
                 * cos ( 2d0 * pi * ( dble(time%day) + 10d0 ) / 365d0 ) )

    ! calculate daylength in hours..
    daylength=calculate_daylength_hours(dble(time%day),latitude_radians)

    ! Estimate diffuse fraction of radiation from measured vs maximum radiation..
    ! (diffuse fraction of radiation from Amthor, 1994) if diffuse fraction has
    ! not been provided
    if (.not.user_opts%met_file_has_sw_diffuse) then
        estfdiff = frac_diffuse_rad( declination , par_top , rtime )
        fdiff = estfdiff
    else
        fdiff = sw_diffuse 
    endif

    ! reset arrays
    leaf_sun = 0d0
    downnir = 0d0 ; upnir  = 0d0 ; downpar    = 0d0 ; uppar        = 0d0 ; downlong = 0d0
    uplong  = 0d0 ; absnir = 0d0 ; abspar_sun = 0d0 ; abspar_shade = 0d0 ; abslong  = 0d0
    do i = 1 , grid%canopy
       temp( i ) = temp_top - 0.1111d0 * dble(i-1) * ( temp_top - temp_bot )
    enddo

    ! detemine light extinction coefficent from sun angle..
    xhour = deg_to_rad * 15d0 * ( rtime - 0.5d0 * period ) * 24d0 / period

    ! sin(beta) - solar geometry affects attenuation
    sun = cos(latitude_radians) * cos(xhour) * cos(declination) + sin(latitude_radians) * sin(declination)

    ! determine extinction coefficient for direct beam radiation..
    ! (but only for light levels above a threshold)
    if ( sun .ge. 0.06d0 ) then
        kbm = 1d0 / ( spherical * sun )
    else
        ! At low sun angles extinction coefficient is zero..
        ! ( Only occurs near sunrise and sunset. Without this )
        ! ( correction we get very unrealistic estimates.     )
        sun = 0d0
        kbm = 0d0
    end if

    if ( sun .gt. 0d0 ) then
        ! If sun has risen then partition par between beam and diffuse
        beam = ( 1d0 - fdiff ) * par_top    ! light attentuation(PPFD) initialise layer above canopy
        diff =         fdiff   * par_top
    else 
        beam = 0d0
        diff = par_top
    end if

    ! ??
    SWR     = sw_rad
    PIb     = 0.48d0 * ( 1d0 - fdiff ) * SWR
    PId     = 0.5d0 * SWR - PIb
    nirbeam = ( 1d0 - fdiff ) * SWR - PIb
    nirdiff = fdiff * SWR - PId

    ! load long wave and change units from W.m-2 -> kW.m-2
    long = lw_rad * 1d-3

    totnir  = nirbeam + nirdiff
    totpar  = beam + diff
    totlong = long

    ! to prevent very low lai values causing precision errors
    if ((sum(lai) + sum(dead_lai)) < 1d-30) then
       lai_in_rad = 0d0
    else
       lai_in_rad = lai + dead_lai
    end if

    ! reset (module) values to zero
    soilabsn = 0d0 ; soil_absorbed_ppfd = 0d0 ; sky_absorbed_ppfd = 0d0 ; skyabsn = 0d0
    soilabsl = 0d0 ; skyabsl = 0d0 ; counter = 0 ; em = 0d0 ; ab = 0d0

    ! determine layer specific reflectances canopy
    ! this impacts crops only
    if (( user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3) .and. user_opts%dead_lai_energy_balance) then
       do i=1,grid%canopy
          if (lafrac(i) > 0d0 .and. lafrac_dead(i) > 0d0) then
             ! if both living and dead fractions are greater than zero then we have
             ! problem
             write(message,*)"Both living and dead leaf area fractions are", &
                             "greater than zero, this is a problem.",        & 
                             "Please check the veg file" 
             call write_log( trim(message) , msg_error , __FILE__ , __LINE__ ) ; stop
          else if (lafrac(i) > 0d0) then
             ! if living leaf area frac present then add living reflectances 
             nirrefl_vector(i) = nirrefl
             parrefl_vector(i) = parrefl
          else if (lafrac_dead(i) > 0d0) then
             ! if dead leaf area frac present then add dead foliage reflectances
             nirrefl_vector(i) = nirrefl_crop
             parrefl_vector(i) = parrefl_crop
          else 
             ! we assume that both must be zero
             nirrefl_vector(i) = 0d0
             parrefl_vector(i) = 0d0
          end if
       enddo ! canopy loop
       ! adjust soil reflectance based on surface litter
       soilnir_corrected=soilnir*soilnir_coef
       soilpar_corrected=soilpar*soilpar_coef
    else
       ! we are not a crop so
       nirrefl_vector=nirrefl
       parrefl_vector=parrefl
       soilnir_corrected=soilnir
       soilpar_corrected=soilpar
    end if ! is this crop pft

    ! start the attenuation routines..
    do while ( counter .lt. 3 )

      ! multiple passes through the canopy - 3 is enough to account for every photon in general
      if ( totpar .gt. 1d0 ) then ! if the sun is up, then call the routines....

        fsnow = 1d0 - exp( - snowheight / 0.1d0)  ! fraction of snow cover on the ground

        ! firstly calculate PAR attenuation - only do sunlit and shaded calculations for PAR..
        reflsurf = (1d0 - fsnow) * soilpar_corrected + fsnow * snowalb_par
        call attenuate_PAR( ff , lai_in_rad , beam , diff , parrefl_vector , partrans , reflsurf , &
                            kbm , clump , abspar_sun , abspar_shade , uppar ,       &
                            downpar , soil_absorbed_ppfd , sky_absorbed_ppfd , sunfrac , sump           )

       ! set leaf_sun after only the first pass, when beam radiation is incident
        if ( counter .eq. 0 ) then
          do i = 1 , grid%canopy
             leaf_sun( i ) = sunfrac( i )
          enddo
        end if

        ! next do Near Infra Red..
        reflsurf = (1d0 - fsnow) * soilnir_corrected + fsnow * snowalb_nir
        call attenuate_NIR( ff , lai_in_rad , nirbeam , nirdiff , nirrefl_vector , nirtrans , reflsurf , kbm , &
                            clump , absnir , upnir , downnir , soilabsn , skyabsn , sumn )

      else
        abspar_sun = 0d0 ; abspar_shade = 0d0
        absnir     = 0d0 ; sunfrac      = 0d0
      end if

      soilt = temp_top    ! use air temp from current timestep

      ! finally calculate the longwave..
      call longwave( ff , lai_in_rad , long , suml , abslong , uplong , downlong , soilabsl , &
                     skyabsl , counter , temp , soilt , em , soem , ab , clump )

      ! reset everything..
      beam = 0d0 ; diff   = 0d0 ; nirbeam = 0d0 ; nirdiff = 0d0
      long = 0d0 ; radnet = 0d0 ; lnet    = 0d0

      ! increment the counter..
      counter = counter + 1 

    enddo
!print*,sum(lai),skyabsl/(sum(lai)*2*0.001*0.96*boltz*(temp_top+273.15)**4)
!print*,sum(lai),skyabsl/(lw_rad*0.001)
!print*,sum(lai),sum(abslong)/(lw_rad*0.001)
!print*,sum(lai),(skyabsn - (sky_absorbed_ppfd / ppfd_to_par))/sw_rad
    ! calculate output for use in other routines
    do i = 1 , grid%canopy
        parabs_sun(i) =   abspar_sun(grid%canopy+1-i)    ! beam par absorbed in each layer
      parabs_shade(i) = abspar_shade(grid%canopy+1-i)    ! diffuse par absorbed in each layer

      ! total shortwave radiation for heating = par + nir , per m2 ground area per layer
      nirabs(i)    = absnir(grid%canopy+1-i)
      longem(i)    = abslong(grid%canopy+1-i) * 1000d0
      lnet         = lnet   + longem(i)
      radnet       = radnet + nirabs(i) 
      totparabs(i) = abspar_sun(grid%canopy+1-i) + abspar_shade(grid%canopy+1-i)
    enddo
!if (sum(lai) > 0d0 .and. sw_rad > 0d0) then
!print*,sum(lai),skyabsn/sw_rad,sum(nirabs)/sw_rad
!print*,sum(lai),(sky_absorbed_ppfd/ppfd_to_par)/sw_rad,sum(totparabs/ppfd_to_par)/sw_rad
!endif
!if (sum(lai) > 0d0)print*,sum(lai),skyabsl/lw_rad,sum(longem)/lw_rad
!if ((sum(nirabs(1:4))/sw_rad) > (sum(parabs_sun(1:4)/ppfd_to_par)+sum(parabs_shade(1:4)/ppfd_to_par))/sw_rad ) print*,"FALSE"
!if (kbm > 0.0) then
!    sw_absorbed = sw_absorbed + (sum(nirabs(1:4))+sum(parabs_sun(1:4)/ppfd_to_par)+sum(parabs_shade(1:4)/ppfd_to_par) ) 
!    sw_total = sw_total + sw_rad 
!    fdiff_total = fdiff_total + fdiff
!    kbm_total = kbm_total + kbm
!endif
!if (time%step == 24) then
!    xhour = deg_to_rad * 15.0 * ( 12.0 - 0.5 * 24.0 ) * 24.0 / 24.0
!    sun = cos(latitude_radians) * cos(xhour) * cos(declination) + sin(latitude_radians) * sin(declination)
!    kbm = 1 / ( spherical * sun )
!    estfdiff = frac_diffuse_rad( declination , par_top , 12d0 )
!write(999,'(f,3(",",f))'),sw_absorbed / sw_total,fdiff_total/daylength,daylength,kbm!,kbm_total/daylength
!    write(999,'(f,4(",",f))'),sw_absorbed / sw_total,estfdiff,daylength,kbm,sum(lai)!,kbm_total/daylength
!    print*,sw_total
!    sw_absorbed = 0d0 ; sw_total = 0d0 ; fdiff_total = 0d0 ; kbm_total = 0d0
!endif
!if (sw_rad > 1)print*,(sum(nirabs)+sum(totparabs/ppfd_to_par))/sw_rad, sum(lai),daylength
!print*,lnet / lw_rad, sum(lai), lnet
!print*,(skyabsl*1000)/lw_rad,sum(lai)
    !  Radiation absorbed by soil (W.m-2); no long emission, isothermal or otherwise
    soilnet = soilabsn + 1000d0 * soilabsl + soil_absorbed_ppfd / ppfd_to_par

    ! generate canopy energy balance output
    check = sum(abspar_sun) + sum(abspar_shade) + soil_absorbed_ppfd + sky_absorbed_ppfd
    if ( user_opts%canopy_csv_output .and. par_top .gt. 0d0 ) then
      call handle_output( 5 , output_data =  (/ soilnet, soilabsn, soilabsl,  &
                        soil_absorbed_ppfd, par_top, fdiff, check, sky_absorbed_ppfd, abspar_sun, &
                        parabs_sun, abspar_shade, parabs_shade, leaf_sun /)   )
    end if
    modRnet  = SWR - skyabsn - sky_absorbed_ppfd / ppfd_to_par - 1d3 * ( skyabsl - totlong )
    checkpar = 0d0

  end subroutine solar
  !
  !----------------------------------------------------------------------
  !
  subroutine solar_update_for_marginal_return(lai_in_rad,leaf_sun_in_rad &
                                             ,apar_sun,apar_shade &
                                             ,anir,along)

    use gv_clim,                only: ppfd_to_par, par_top, sw_rad, sw_diffuse, lw_rad &
                                     ,temp_bot, temp_top 
    use gv_irradiance_sunshade, only: skyabsl
    use gv_scale_declarations,  only: boltz, grid, pi, time, user_opts, deg_to_rad
    use gv_snow_info,           only: snowalb_nir, snowalb_par, snowheight,fsnow
    use gv_veg,                 only: latitude_radians, lafrac,lafrac_dead 
    use spa_io,                 only: handle_output
    use log_tools

    implicit none

    ! arguments
    double precision, dimension(grid%canopy), intent(inout) :: lai_in_rad &
                                                  ,leaf_sun_in_rad,along,anir  &
                                                  ,apar_sun,apar_shade
    ! local variables..
    integer :: counter, i
    double precision    :: ab, beam, clump, diff, em, estfdiff, ff, kbm, long,        &
               nirbeam, nirdiff, period, PIb, PId, rtime, skyabsn, soem,  &
               soilt, suml, sumn, sump, sun, SWR,totlong, totnir, totpar, xhour, &
               soilnir_corrected, soilpar_corrected,   &
               reflsurf

    double precision,dimension(grid%canopy) :: abslong, abspar_shade, abspar_sun, absnir, sunfrac, temp, &
                                   parrefl_vector, nirrefl_vector

    ! include an extra value each side bracketing grid%canopy..
    double precision,dimension(0:grid%canopy+1) :: downlong, downnir, downpar, uplong, upnir, uppar

    ! calculations..
    period    = dble(time%steps_per_day)
    clump     = 1d0   ! clumping factor - =1 means unclumped, random distribution
    ff        = 1d0   ! correction factor for stems etc on attenuation, currently unused and untested, set at 1.0
    rtime     = dble( time%step )

    ! Estimate diffuse fraction of radiation from measured vs maximum radiation..
    ! (diffuse fraction of radiation from Amthor, 1994) if diffuse fraction has
    ! not been provided
    if (.not.user_opts%met_file_has_sw_diffuse) then
       estfdiff = frac_diffuse_rad( declination , par_top , rtime )
       fdiff = estfdiff
    else
       fdiff = sw_diffuse 
    endif

    ! reset arrays
    leaf_sun = 0d0
    downnir = 0d0 ; upnir  = 0d0 ; downpar    = 0d0 ; uppar        = 0d0 ; downlong = 0d0
    uplong  = 0d0 ; absnir = 0d0 ; abspar_sun = 0d0 ; abspar_shade = 0d0 ; abslong  = 0d0
    do i = 1 , grid%canopy
       temp( i ) = temp_top - 0.1111d0 * dble(i-1) * ( temp_top - temp_bot )
    enddo

    ! detemine light extinction coefficent from sun angle..
    xhour = deg_to_rad * 15d0 * ( rtime - 0.5d0 * period ) * 24d0 / period

    ! sin(beta) - solar geometry affects attenuation
    sun = cos(latitude_radians) * cos(xhour) * cos(declination) + sin(latitude_radians) * sin(declination)

    ! determine extinction coefficient for direct beam radiation..
    ! (but only for light levels above a threshold)
    if ( sun .ge. 0.06d0 ) then
        kbm = 1 / ( spherical * sun )
    else
        ! At low sun angles extinction coefficient is zero..
        ! ( Only occurs near sunrise and sunset. Without this )
        ! ( correction we get very unrealistic estimates.     )
        sun = 0d0
        kbm = 0d0
    end if

    if ( sun .gt. 0d0 ) then
        ! If sun has risen then partition par between beam and diffuse
        beam = ( 1d0 - fdiff ) * par_top    ! light attentuation(PPFD) initialise layer above canopy
        diff =         fdiff   * par_top
    else 
        beam = 0d0
        diff = par_top
    end if

    ! ??
    SWR     = sw_rad
    PIb     = 0.48d0 * ( 1d0 - fdiff ) * SWR
    PId     = 0.5d0 * SWR - PIb
    nirbeam = ( 1d0 - fdiff ) * SWR - PIb
    nirdiff = fdiff * SWR - PId

    ! load long wave and change units from W.m-2 -> kW.m-2
    long = lw_rad * 1d-3

    totnir  = nirbeam + nirdiff
    totpar  = beam + diff
    totlong = long

    ! reset (module) values to zero
    soilabsn = 0d0 ; soil_absorbed_ppfd = 0d0 ; sky_absorbed_ppfd = 0d0 ; skyabsn = 0d0
    soilabsl = 0d0 ; skyabsl = 0d0 ; counter = 0 ; em = 0d0 ; ab = 0d0

    ! determine layer specific reflectances canopy
    ! this impacts crops only
    if (( user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3) .and. user_opts%dead_lai_energy_balance) then
       do i=1,grid%canopy
          if (lafrac(i) > 0d0 .and. lafrac_dead(i) > 0d0) then
             ! if both living and dead fractions are greater than zero then we have
             ! problem
             write(message,*)"Both living and dead leaf area fractions are", &
                             "greater than zero, this is a problem.",        & 
                             "Please check the veg file" 
             call write_log( trim(message) , msg_error , __FILE__ , __LINE__ ) ; stop
          else if (lafrac(i) > 0d0) then
             ! if living leaf area frac present then add living reflectances 
             nirrefl_vector(i) = nirrefl
             parrefl_vector(i) = parrefl
          else if (lafrac_dead(i) > 0d0) then
             ! if dead leaf area frac present then add dead foliage reflectances
             nirrefl_vector(i) = nirrefl_crop
             parrefl_vector(i) = parrefl_crop
          else 
             ! we assume that both must be zero
             nirrefl_vector(i) = 0d0
             parrefl_vector(i) = 0d0
          end if
       enddo ! canopy loop
       ! adjust soil reflectance based on surface litter
       soilnir_corrected = soilnir*soilnir_coef
       soilpar_corrected = soilpar*soilpar_coef
    else
       ! we are not a crop so
       nirrefl_vector = nirrefl
       parrefl_vector = parrefl
       soilnir_corrected = soilnir
       soilpar_corrected = soilpar
    end if ! is this crop pft

    ! start the attenuation routines..
    do while ( counter .lt. 3 )

      ! multiple passes through the canopy - 3 is enough to account for every photon in general
      if ( (totpar) .gt. 1d0 ) then    ! if the sun is up, then call the routines....

        ! firstly calculate PAR attenuation - only do sunlit and shaded calculations for PAR..
        reflsurf = (1d0 - fsnow) * soilpar_corrected + fsnow * snowalb_par
        call attenuate_PAR( ff , lai_in_rad , beam , diff , parrefl_vector , partrans , reflsurf , &
                            kbm , clump , abspar_sun , abspar_shade , uppar ,       &
                            downpar , soil_absorbed_ppfd , sky_absorbed_ppfd , sunfrac , sump           )

       ! set leaf_sun after only the first pass, when beam radiation is incident
        if ( counter .eq. 0 ) then
            do i = 1 , grid%canopy
               leaf_sun_in_rad( i ) = sunfrac( i )
            enddo
        end if

        ! next do Near Infra Red..
        reflsurf = (1d0 - fsnow) * soilnir_corrected + fsnow * snowalb_nir
        call attenuate_NIR( ff , lai_in_rad , nirbeam , nirdiff , nirrefl_vector , nirtrans , reflsurf , kbm , &
                            clump , absnir , upnir , downnir , soilabsn , skyabsn , sumn       )
      else
        abspar_sun = 0d0 ; abspar_shade = 0d0
        absnir     = 0d0 ; sunfrac      = 0d0
      end if

      ! use current air temperature for soil temperature in calculating longwave
      ! balance
      soilt = temp_top  

      ! finally calculate the longwave..
      call longwave( ff , lai_in_rad , long , suml , abslong , uplong , downlong , soilabsl , &
                     skyabsl , counter , temp , soilt , em , soem , ab , clump         )

      ! reset everything..
      beam = 0d0 ; diff = 0d0 ; nirbeam = 0d0 ; nirdiff = 0d0
      long = 0d0

      ! increment the counter..
      counter = counter + 1 

    enddo ! for do while counter

    ! calculate output for use in other routines
    do i = 1 , grid%canopy
       apar_sun(i)   = abspar_sun(grid%canopy+1-i)    ! beam par absorbed in each layer
       apar_shade(i) = abspar_shade(grid%canopy+1-i)    ! diffuse par absorbed in each layer
       ! total shortwave radiation for heating = par + nir , per m2 ground area per layer
       anir(i)    = absnir(grid%canopy+1-i)
       along(i)    = abslong(grid%canopy+1-i) * 1000d0
    enddo

  end subroutine solar_update_for_marginal_return
  !
  !----------------------------------------------------------------------
  !
  subroutine solar_update_for_iterative_canopy

  ! Redistributed long wave radiation through the canopy at the end of each
  ! canopy iteration

    use gv_clim,                only: lw_rad
    use gv_irradiance_sunshade, only: skyabsl
    use gv_scale_declarations,  only: boltz, grid, pi, time_holder
    use gv_veg,                 only: lai, dead_lai, modrnet, leaf_temp_intermediate
    use gv_hourscale,           only: hourts
    use spa_io,                 only: handle_output

    implicit none

    ! local variables..
    integer :: counter, i
    double precision    :: ab, clump, em, ff, lnet, long, &
               radnet, soem, &
               soilt, suml, totlong

    double precision,dimension(grid%canopy) :: abslong, temp, lai_in_rad

    ! include an extra value each side bracketing grid%canopy..
    double precision,dimension(0:grid%canopy+1) :: downlong, uplong

    ! calculations..
    clump     = 1d0    ! clumping factor - =1 means unclumped, random distribution
    ff        = 1d0   ! correction factor for stems etc on attenuation, currently unused and untested, set at 1.0.

    ! reset diffuse arrays
    downlong = 0d0 ; uplong = 0d0 ; abslong = 0d0
    do i = 1 , grid%canopy
       temp( i ) = leaf_temp_intermediate(i)
    enddo

    ! determine longwave incident from sky (Jones p.27) kW m-2 
    ! load long wave and change units from W.m-2 -> kW.m-2
    long = lw_rad * 1d-3
    totlong = long

    ! hack to prevent very low lai values creating NaN values during
    ! calculations
    if ( (sum(lai) + sum(dead_lai) ) < 1d-30 ) then
        lai_in_rad = 0d0
    else
        lai_in_rad = lai + dead_lai
    end if

    ! remove isothermal long wave from system net radiation
    modrnet=modrnet+1e3*(skyabsl-totlong)

    ! reset (module) values to zero
    soilabsl = 0d0 ; skyabsl = 0d0 ; counter = 0 ; em = 0d0 ; ab = 0d0

    ! start the attenuation routines..
    do while ( counter .lt. 3 )
      soilt = hourts-273.15d0    ! soil temperature from current step (K->oC)
      ! finally calculate the longwave..
      call longwave( ff , lai_in_rad , long , suml , abslong , uplong , downlong , soilabsl , &
                     skyabsl , counter , temp , soilt , em , soem , ab , clump )
      ! reset everything..
      long = 0d0 ; radnet = 0d0 ; lnet = 0d0
      ! increment the counter..
      counter = counter + 1
    enddo ! while counter

    ! update model radiation net with new long wave variables
    modrnet = modrnet-1d3*(skyabsl-totlong)

    ! calculate output for use in other routines
    do i = 1 , grid%canopy
       ! total shortwave radiation for heating = par + nir , per m2 ground area
       ! per layer
       longem(i) = abslong(grid%canopy+1-i) * 1000d0
       lnet      = lnet + longem(i)
    enddo

  end subroutine solar_update_for_iterative_canopy
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  subroutine attenuate_NIR( ff , lai , bm , df , refl , trans , reflsoil , kbm , clump , &
                            absorb , updf , downdf , soilrad , skyrad , sum_rad )

    use gv_scale_declarations, only: grid

    implicit none

    ! arguments..
    double precision,intent(in)    :: bm, clump, df, ff, kbm, lai(grid%canopy), refl(grid%canopy), reflsoil, trans
    double precision,intent(inout) :: absorb(grid%canopy), downdf(0:grid%canopy+1), skyrad, soilrad, updf(0:grid%canopy+1)
    double precision,intent(out)   :: sum_rad

    ! local variables..
    integer :: i
    double precision :: decay, intbm, intdf, kdf, leafrad
    double precision    :: beam(0:grid%canopy+1)


    ! calculations..
      beam(grid%canopy+1) = bm
    downdf(grid%canopy+1) = df

    ! diffuse radiation approximated by beta=30 degrees
    kdf = 0.5d0

    do i = grid%canopy , 1 , -1

       ! determine unintercepted radiation
       beam(i) = beam(i+1) * exp( -kbm * clump * lai(grid%canopy+1-i) )
       ! and thus intercepted radiation
       intbm = beam(i+1) - beam(i)
       ! now for diffuse
       decay = exp( -kdf * clump * lai(grid%canopy+1-i) )
       downdf(i) = downdf(i) + downdf(i+1) * decay
       intdf = downdf(i+1) * ( 1. - decay )
       ! correct for transmittance (trans)
       absorb(i) = absorb(i) + intbm * ( 1d0 - trans - refl(i) ) / ff
       ! correct for reflectance (refl)
       absorb(i) = absorb(i) + intdf * ( 1d0 - trans - refl(i) ) / ff
       ! add transmitted beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + trans * ( intbm + intdf )
       ! add reflected beam & diffuse to upward diffuse
       updf(i) = updf(i) + refl(i) * ( intbm + intdf )

    enddo

    ! reflectance by soil of radiation that penetrates all vegetation
    updf(0) = reflsoil * ( beam(1) + downdf(1) )
    soilrad = soilrad + ( 1d0 - reflsoil ) * ( beam(1) + downdf(1) )

    ! reset the downwelling diffuse radiation..
    downdf = 0d0

    do i = 1 , grid%canopy

       ! now return upwards through the canopy, dealing with reflected radiation
       ! unintercepted
       decay = exp( -kdf * clump * lai(grid%canopy+1-i) )
       updf(i) = updf(i) + updf(i-1) * decay
       ! intercepted
       intdf = updf(i-1) * ( 1. - decay )
       ! absorbed determined from transmittance/reflectance
       absorb(i) = absorb(i) + intdf * ( 1. - trans - refl(i) ) / ff
       ! add reflected beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + refl(i) * intdf
       ! add transmitted beam & diffuse to upward diffuse
       updf(i)   = updf(i) + trans * intdf
       updf(i-1) = 0d0

    enddo

    ! accumulate radiation back to the sky
    skyrad = skyrad + updf(grid%canopy)
    ! reset upwards moving radation
    updf(grid%canopy) = 0d0
    ! estimate total absorbed by canopy
    leafrad = sum( absorb(1:grid%canopy) )
    ! sum components for later checking
    sum_rad = soilrad + skyrad + leafrad

  end subroutine attenuate_NIR
  !
  !----------------------------------------------------------------------
  !
  subroutine attenuate_PAR( ff , lai , bm , df , refl , trans , reflsoil ,  &
                            kbm , clump , absorb_sun, absorb_shade , updf , &
                            downdf , soilrad , skyrad , sunfrac , sum_rad   )

    use gv_scale_declarations, only: grid

    implicit none

    ! arguments..
    double precision,intent(in)    :: clump, bm, df, ff, kbm, lai(grid%canopy), refl(grid%canopy), reflsoil, trans
    double precision,intent(inout) :: absorb_sun(grid%canopy), absorb_shade(grid%canopy), downdf(0:grid%canopy+1), &
                          skyrad, soilrad, updf(0:grid%canopy+1)
    double precision,intent(out)   :: sum_rad, sunfrac(grid%canopy)

    ! local variables..
    integer :: i
    double precision    :: beam(0:grid%canopy+1), cumlai, decay, intbm, intdf, kdf, leafrad, suncum, sunla(grid%canopy), sunprev


    ! calculations..
    beam(grid%canopy+1) = bm
    downdf(grid%canopy+1) = df

    ! initial values..
    cumlai = 0d0 ; sunprev = 0d0

    ! For ease, we approximate diffuse radiation by beta=30 degrees,
    ! which gives an extinction coefficient of..
    kdf = 0.5d0

    ! Calculate how the radiation is diffused as it travels down
    ! through the canopy layers..
    do i = grid%canopy , 1 , -1

      ! If there is sunlight and some extinction in this layer then..
      if ( ( kbm .gt. 0d0 ) .and. ( bm .gt. 0d0 ) ) then

         ! first calculate cumulative leaf area..
         cumlai = cumlai + lai( grid%canopy+1 - i )
         ! then total sunlit for cumlai. ( kbm = 1/(2sinBeta) )
         suncum = ( 1d0 - exp( -kbm * cumlai ) ) / kbm
         ! sunlit area in this layer is cumulative sunlit minus previous
         sunla(grid%canopy+1-i) = suncum - sunprev
         ! save previous
         sunprev = suncum  

         ! determine sunlight fraction
         if ( sunla( grid%canopy+1 - i ) .gt. 0d0 ) then
             sunfrac( grid%canopy+1 - i ) = sunla( grid%canopy+1 - i ) / lai( grid%canopy+1 - i )
         else
             sunfrac( grid%canopy+1 - i ) = 0d0
         end if
         ! Beam radiation on sunlit leaves is constant
         beam( i ) = bm  

         ! Intercepted radiation dI = Io.k.Ls
         ! (dI=intercepted rad, Io=downwelling beam, Ls =sunlit leaf area, k=extinction)
         intbm = bm * kbm * sunla( grid%canopy+1 - i )

      else

         ! no sunlight..
         intbm      = 0d0 ! intercepted radiation
         beam(i)    = 0d0 !
         sunfrac(i) = 0d0 ! fraction of sunlit leaf

      end if

      ! now for diffuse
      decay     = exp( -kdf * clump * lai(grid%canopy+1-i) ) ! attenuation factor
      downdf(i) = downdf(i) + downdf(i+1) * decay            ! diffuse radiation that passes through layer
      intdf     = downdf(i+1) * ( 1d0 - decay )              ! intercepted diffuse radiation in each layer

      ! Absorption of direct beam (correct interception for transmittance and reflectance)..
      absorb_sun(i)   = absorb_sun(i)   + intbm * ( 1d0 - trans - refl(i) ) / ff

      ! Absorption of diffuse..
      absorb_shade(i) = absorb_shade(i) + intdf * ( 1d0 - trans - refl(i) ) / ff

      ! add transmitted beam & diffuse to downward diffuse
      downdf( i ) = downdf( i ) + trans * ( intbm + intdf )

      ! add reflected beam & diffuse to upward diffuse
      updf( i )   = updf( i ) + refl(i) * ( intbm + intdf )


    enddo

    if ( bm .gt. 0d0 ) then
        ! Direct beam radiation reaching soil surface is directly calculated by this eqn...
        beam(1) = ( 1d0 / kbm - suncum ) * bm * kbm
    end if

    ! reflectance by soil of radiation that penetrates all vegetation
    updf(0) = reflsoil * ( beam(1) + downdf(1) )
    soilrad = soilrad + ( 1d0 - reflsoil ) * ( beam(1) + downdf(1) )

    ! reset downwelling diffuse radiation..
    downdf = 0d0

    do i = 1 , grid%canopy

      ! now return upwards through the canopy, dealing with reflected radiation
      ! unintercepted
      decay = exp( -kdf * clump * lai(grid%canopy+1-i) )
      updf(i) = updf(i) + updf(i-1) * decay

      ! intercepted
      intdf = updf( i - 1 ) * ( 1d0 - decay )

      ! absorbed determined from transmittance/reflectance
      absorb_shade(i) = absorb_shade(i) + intdf * ( 1d0 - trans - refl(i) ) / ff

      ! add reflected beam & diffuse to downward diffuse
      downdf(i) = downdf(i) + refl(i) * intdf


      ! add transmitted beam & diffuse to upward diffuse
      updf(i)   = updf(i) + trans * intdf
      updf(i-1) = 0d0

    enddo
    
    ! accumulate radiation lost back to the sky
    skyrad  = skyrad + updf(grid%canopy)
    ! reset upward radiation 
    updf(grid%canopy) = 0d0
    ! determine radiation absorbed by canopy
    leafrad = 0d0
    leafrad = sum( absorb_sun(1:grid%canopy) ) + sum( absorb_shade(1:grid%canopy) )
    ! accumulate radiation components for later energy balance checking
    sum_rad = soilrad + skyrad + leafrad

  end subroutine attenuate_PAR
  !
  !----------------------------------------------------------------------
  !
  double precision function frac_diffuse_rad( dec , parsteps , rtime )
    
    ! Determines the ratio of actual-to-potential radiation (ksko)   !
    ! for the day, from the total potential daily PAR, and then uses !
    ! a relationship from Erbs et al (1982) to estimate fraction of  !
    ! incoming radiation that is diffuse.                            !

    use gv_scale_declarations, only: pi, time, deg_to_rad
    use gv_veg, only: latitude_radians

    implicit none

    ! arguments..
    double precision,intent(in) :: dec, parsteps, rtime

    ! local variables..
    double precision :: ff, hourangle, ko, ks, ksko

    ! calculations..
    ks        = parsteps / 2.3d0       ! (MJ m-2 d-1)
    hourangle = deg_to_rad * 15d0 * ( rtime - 0.5d0 * dble(time%steps_per_day) ) * 24d0 / dble(time%steps_per_day)

    ! extraterrestrial radiation (J m-2)
    ko        = So * ( sin( dec ) * sin( latitude_radians ) + cos( latitude_radians ) * cos( dec ) * cos( hourangle ) )
    ksko      = ks / ko  ! hourly ksko

    ! Diffuse fraction
    if (ksko .lt. 0d0) then
        ff = 1d0
    else if ( ksko .lt. 0.22d0 ) then
        ff = 1d0 - 0.09d0 * ksko
    else if ( ksko .lt. 0.8 ) then
        ff = 0.9511d0                  &
            - 0.1604d0 * ksko          &
             + 4.388d0 * ksko**2       &
              - 16.638d0 * ksko**3     &
               + 12.336d0 * ksko**4
    else 
        ff = 0.165d0
    end if
    frac_diffuse_rad = ff 

  end function frac_diffuse_rad
  !
  !----------------------------------------------------------------------
  !
  subroutine longwave( ff, lai, df, sum_rad, absorb, updf, downdf, soilrad, skyrad, count, & 
                       Ta, Ts, totemit, soilem, totab, clump)

    use gv_scale_declarations, only: boltz, grid
    use gv_veg,                only: emiss

    implicit none

    ! arguments..
    integer,intent(in) :: count
    double precision,intent(in)    :: clump, df, ff, lai(grid%canopy), Ta(grid%canopy), Ts
    double precision,intent(inout) :: absorb(grid%canopy), downdf(0:grid%canopy+1), skyrad, soilrad, totemit, updf(0:grid%canopy+1)
    double precision,intent(out)   :: soilem, sum_rad, totab

    ! local variables..
    integer :: i
    double precision    :: decay, eup, edown, intdf, kdf, leafrad

    downdf(grid%canopy+1) = df

    ! diffuse radiation approximated by beta=30 degrees
    kdf = 0.5d0
    do i = grid%canopy , 1 , -1

       decay = exp( -kdf * clump * lai(grid%canopy+1-i) )

       if ( count .eq. 0 ) then

         ! longwave radiative heat loss from top side of leaf (KW m-2)
         ! The '*(1-decay)' corrects emissions into the 1-D, vertical
         eup = 0.001d0 * emiss * boltz * ( Ta(grid%canopy+1-i) + 273.15d0 )**4 * ( 1d0 - decay )
         totemit = totemit + 2d0 * eup * lai(grid%canopy+1-i)

       else

         eup = 0d0

       end if

       ! and bottom side..
       edown     = eup
       downdf(i) = downdf(i) + downdf(i+1) * decay + edown
       intdf     = downdf(i+1) * ( 1d0 - decay )
       ! correct for transmittance (trans) and reflectance (refl) & PAI/LAI ratio
       absorb(i) = absorb(i) + intdf * emiss / ff - eup - edown
       ! add transmitted diffuse to downward diffuse
       downdf(i) = downdf(i) + 0.5d0 * ( 1d0 - emiss ) * intdf
       ! add reflected diffuse to upward diffuse
       updf(i) = updf(i) + 0.5d0 * ( 1d0 - emiss ) * intdf + eup
       totab = totab + intdf * emiss
    enddo

    ! reflectance by soil of radiation that penetrates all vegetation
    ! Note that 1-emissivity = 0.04
    soilem = emiss * boltz * ( Ts + 273.15d0 )**4 * 0.001d0
    if ( count .eq. 0 ) then
        updf(0) = 0.04d0 * downdf(1) + soilem
    else
        updf(0) = 0.04d0 * downdf(1)
    end if
    soilrad = soilrad + emiss * downdf(1)

    ! reset downwelling diffuse radiation..
    downdf = 0d0

    do i = 1 , grid%canopy

       ! now return upwards through the canopy, dealing with reflected radiation
       ! unintercepted
       decay   = exp( -kdf * clump * lai(grid%canopy+1-i) )
       updf(i) = updf(i) + updf(i-1) * decay
       ! intercepted
       intdf = updf(i-1) * ( 1d0 - decay )
       ! absorbed determined from transmittance/reflectance
       absorb(i) = absorb(i) + intdf * emiss / ff
       ! add reflected beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + 0.02d0 * intdf
       ! add transmitted beam & diffuse to upward diffuse
       updf(i)   = updf(i) + 0.02d0 * intdf
       updf(i-1) = 0d0
       totab     = totab + intdf * emiss

    enddo

    skyrad = skyrad + updf(grid%canopy)
    updf(grid%canopy) = 0d0
    leafrad = sum( absorb(1:grid%canopy) )
    sum_rad = soilrad + skyrad + leafrad

  end subroutine longwave
  !
  !----------------------------------------------------------------------
  !
end module light
!
!------------------------------------------------------------------------
!
