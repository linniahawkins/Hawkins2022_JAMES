! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module soil_air

  !! > this module needs a summary < !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: saxton_parameters, soil_conductivity, soil_porosity, soil_resistance, &
            soil_water_potential, water_retention, water_uptake_layer, &
            calculate_canopy_soil_resistance, roots
  ! Variables..
  public :: drainlayer, liquid, slayer, unsat


  integer :: slayer,    &  ! Soil layer variable used when passing from a subroutine to a function which is acting
                           !  on a specific soil layer.
     water_retention_pass  ! Count of number of iterations gone through for soil layer in the Saxton equations

  double precision :: drainlayer, & ! The field capacity of the specific layer being worked on in the s/r soil_balance
                          liquid, & ! Liquid fraction of the water present within a given soil layer, used in s/r soil_balance
                           unsat    ! Unsaturated volume of the soil layer below the current one being modelled (m3 m-2).
                                    !  Calculated each timestep.

  save

contains
  !
  !----------------------------------------------------------------------
  !
  subroutine saxton_parameters()

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: cond1, cond2, cond3, potA, potB

    implicit none

    integer :: i
    double precision, parameter :: A = -4.396d0,  B = -0.0715d0, CC = -4.880d-4, D = -4.285d-5, &
                                   E = -3.140d0,  F = -2.22d-3,   G = -3.484d-5, H = 0.332d0,   & 
                                   J = -7.251d-4, K = 0.1276d0,   P = 12.012d0,  Q = -7.551d-2, &
                                   R = -3.895d0,  T = 3.671d-2,   U = -0.1103d0, V = 8.7546d-4, &
                                   mult1 = 100d0, mult2 = 2.778d-6, mult3 = 1000d0

    ! layed out in this manor to avoid memory management issues in module
    ! variables
    potA = A + (B * soil_frac_clay) + (CC * soil_frac_sand * soil_frac_sand) + &
               (D * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    do i = 1, grid%core
       potA(i)  = exp(potA(i))
    end do
    potA = potA * mult1
    potB  = E + (F * soil_frac_clay * soil_frac_clay) + (G * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    cond1 = mult2
    cond2 = P + (Q * soil_frac_sand)
    cond3 = R + (T * soil_frac_sand) + (U * soil_frac_clay) + (V * soil_frac_clay * soil_frac_clay)

  end subroutine saxton_parameters
  !
  !----------------------------------------------------------------------
  !
  double precision function soil_conductivity( wf )

    ! Used in the soil drainage integrator. !
    ! Returns a single-point value.         !
    ! 'slayer' is a module variable that    !
    !  provides the soil-layer number.      !

    use gv_soil_structure, only: cond1, cond2, cond3

    implicit none

    ! arguments..
    double precision, intent(in) :: wf ! fraction of water in soils

    if ( wf .lt. 0.05d0 ) then    ! Avoid floating-underflow fortran error
      soil_conductivity = 1d-30
    else
      soil_conductivity = cond1(slayer) * exp( cond2(slayer) + cond3(slayer) / wf )
      ! Soil conductivity (m s-1 )
    end if

  end function soil_conductivity
  !
  ! ---------------------------------------------------------------------
  !
  subroutine soil_porosity

   ! Porosity, the fraction of empty space in the soil,
   ! is estimated from Saxton equations. 

    use gv_hydrol,             only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: porosity

    implicit none

    ! local variables..
    integer :: i
    double precision    :: H, J, K

    ! saxton params relevant to porosity..
    H = 0.332d0 ; J = -7.251d-4 ; K = 0.1276d0

    ! loop over soil layers..
    do i = 1 , grid%core
        if (i .lt. 2) then
          porosity(i) = 0.2d0
        elseif (i .lt. 4) then
          porosity(i) = 0.25d0
        elseif (i .lt. 9) then
          porosity(i) = 0.3d0 
        else
          porosity(i) = 0.35d0
          !porosity(i) = H + J * soil_frac_sand(i) + K * log10( soil_frac_clay(i) )
        end if
    enddo
  end subroutine soil_porosity
  !
  ! ---------------------------------------------------------------------
  !
  subroutine soil_resistance(soilR_out,soilR1_out,soilR2_out &
                            ,root_length_in,root_mass_in,rooted_layers_in)

    !Calculates soil and root root hydraulic resistances based on root
    !distribution information, This information is combined to give the total
    !soil root hydraulic resistance used in photosynthetic calculations.

    use gv_meteo,              only: head
    use gv_scale_declarations, only: pi, mol_to_g_water, grid
    use gv_soil_structure,     only: abovebelow, conduc, thickness
    use gv_veg,                only: rootresist, root_radius  ! some parameters

    implicit none

    ! arguments
    integer, intent(in) :: rooted_layers_in
    double precision, dimension(grid%core), intent(in) :: root_length_in, &
                                                          root_mass_in
    double precision, dimension(grid%core), intent(out) :: soilR_out, & !
                                                          soilR1_out, & !
                                                          soilR2_out    !

    ! local variables..
    integer :: i
    double precision    :: Lsoil, rs, rs2

    soilR_out = 0d0

    ! Calculate soil-root hydraulic resistance
    do i = 1 , rooted_layers_in
      Lsoil = conduc(i) / head    !converts from ms-1 to m2 s-1 MPa-1
      if ( Lsoil .lt. 1d-35 ) then    !prevent floating point error
        soilR_out(i) = 1d35
      else 
        rs  = sqrt( 1d0 / (root_length_in(i) * pi ) )
        rs2 = log( rs / root_radius ) / ( 2d0 * pi * root_length_in(i) * thickness(i) * Lsoil )
        ! soil water resistance
        soilR1_out(i) = rs2 * 1d-6 * mol_to_g_water * 0.001d0    ! convert from MPa s m2 m-3 to MPa s m2 mmol-1
        !second component of below ground resistance related to root hydraulics
        !rrcheck=rootresist/(root_mass(i)*thickness(i)/abovebelow)
        ! root reistance
        soilR2_out(i) = rootresist / ( root_mass_in(i) * thickness(i) / abovebelow )
        soilR_out(i)  = soilR1_out(i) + soilR2_out(i) ! MPa s m2 mmol-1
      end if
    enddo

  end subroutine soil_resistance
  !
  ! ---------------------------------------------------------------------
  !
  subroutine soil_water_potential()

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    use gv_scale_declarations, only: user_opts
    use gv_soil_structure,     only: potA, potB, rooted_layers, SWP, waterfrac &
                                     , swp_params

    implicit none

    ! local variables..
    integer :: i
    double precision    :: soil_wp,a,b,c,d

    a = swp_params( 1 )
    b = swp_params( 2 )
    c = swp_params( 3 )
    d = swp_params( 4 )

    do i = 1 , rooted_layers
       if ( waterfrac(i) .gt. 0d0 ) then
           if (user_opts%soil_hydrology_opts ==2) then ! user defined
               soil_wp = a - b/(1d0+exp((waterfrac(i)-c)/d))
           elseif (user_opts%soil_hydrology_opts ==3) then ! Rheur et al., 2014
               soil_wp = -0.4d0 - 1.6d0*1d0/(1d0+exp((waterfrac(i)-0.096d0)/0.0184d0))
           else
               soil_wp = -0.001d0 * potA(i) * waterfrac(i)**potB(i)   !  Soil water potential (MPa)
           end if
       else
           soil_wp = -9999d0
       end if
       SWP(i) = soil_wp
    enddo

  end subroutine soil_water_potential
  !
  !----------------------------------------------------------------------
  !
  subroutine water_retention

    ! field capacity calculations for saxton eqns !

    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: field_capacity,potB
    use math_tools,            only: zbrent

    implicit none

    ! local variables..
    integer        :: i
    double precision           :: x1, x2
    double precision,parameter :: xacc = 0.0001d0

    do i = 1 , grid%core
      x1   = 0.1d0    ! low guess
      x2   = 0.7d0    ! high guess
      water_retention_pass = i
      ! field capacity is water content at which SWP = -10 kPa
      field_capacity( i ) = zbrent( 'water_retention:water_retention_saxton_eqns' , water_retention_saxton_eqns , x1 , x2 , xacc )
    enddo

  end subroutine water_retention
  !
  !----------------------------------------------------------------------
  !
  subroutine water_uptake_layer( total_est_evap )

    ! Find from which layer water is withdrawn !

    use gv_scale_declarations, only: grid, time
    use gv_soil_structure,     only: fraction_uptake,iceprop,rooted_layers,soilR,SWP,weighted_SWP,potB!,soilR1
    use gv_veg,                only: canopy_soil_resistance, lai, minlwp, gplant
    use gv_daily_averages
    use log_tools

    implicit none

    ! arguments..
    double precision,intent(out) :: total_est_evap

    ! local variables..
    integer          :: i, p
    double precision :: est_evap(grid%core), frac

    ! -- calculations begin below --

    frac = 0d0 ! local var so should be initialised to zero anyway, but just to be sure
    total_est_evap = 0d0 ; weighted_SWP    = 0d0
    est_evap       = 0d0 ; fraction_uptake = 0d0

    do i = 1 , rooted_layers
       ! estimate max transpiration (mmol.m-2.s-1) from gradient-gravity / soil resistance.
       est_evap(i) = ( SWP(i) - minlwp ) / soilR(i)
       est_evap(i) = max( 0d0 , est_evap(i) )       ! no negative 
       if ( iceprop(i) .gt. 0d0 ) est_evap(i) = 0d0 ! no uptake from frozen soils
    enddo
    total_est_evap = sum( est_evap )

    ! weighted soil water potential
    if ( total_est_evap .gt. 0d0 ) then
       ! Water was evaporated from some layers..
       do i = 1 , rooted_layers
          weighted_SWP = weighted_SWP + SWP(i) * est_evap(i)
          ! fraction of total et taken from layer i...
          fraction_uptake(i) = est_evap(i) / total_est_evap
       enddo
       weighted_SWP = weighted_SWP / total_est_evap
    else
       ! No water was evaporated..
       fraction_uptake(:) = 1d0 / dble(rooted_layers)
    end if

    ! prepare some output for training
    if (allocated(daily_weighted_SWP)) then
        daily_weighted_SWP(time%step)=weighted_SWP
!        daily_weighted_soilR(time%step)=sum(soilR1(1:rooted_layers)*fraction_uptake(1:rooted_layers))
        daily_weighted_soilR(time%step)=((weighted_SWP-minlwp)/total_est_evap) & 
                                       + (1d0/(gplant*sum(lai)))
    endif

    if ( nint ( sum( fraction_uptake ) ) .ne. 1d0 .and. time%steps_count > 1 .and. total_est_evap > 0d0) then
        call write_log( 'The sum of uptake fraction is not (nearly) equal to 1 '&
                      //' in water_uptake_layer' , msg_warning , __FILE__ , __LINE__ )
    end if
    if ( ( fraction_uptake(1) .gt. 1d0 ) .or. ( fraction_uptake(1) .lt. 0d0 ) ) then
        call write_log( 'Problem with the uptake fraction (either >1 or 0<)' , &
                        msg_warning , __FILE__ , __LINE__ )
    end if

    ! calculate combined canopy_soil_resistance
    call calculate_canopy_soil_resistance(canopy_soil_resistance,soilR,lai,rooted_layers)

  end subroutine water_uptake_layer
  !
  !---------------------------------------------------------------------- 
  !
  subroutine calculate_canopy_soil_resistance(canopy_soil_resistance_out,soilR_in,lai_in,rooted_layers_in)

    use gv_scale_declarations, only: time, grid

    ! calculate the canopy layer specific soil+root hydraulic resistances

    ! arguments
    integer, intent(in) :: rooted_layers_in
    double precision, dimension(grid%core), intent(in):: soilR_in
    double precision, dimension(grid%canopy), intent(in) :: lai_in
    double precision, dimension(grid%canopy), intent(out) :: canopy_soil_resistance_out

    ! local variables 
    integer :: p, i
    double precision :: frac

    canopy_soil_resistance_out = 0d0    ! reset
    do p = 1 , grid%canopy
      frac = 0d0
      if (lai_in(p) > 0d0) frac = lai_in(p) / sum(lai_in)
      if (frac > 0d0) then
          do i = 1, rooted_layers_in
             ! soil resistance for each canopy layer is related to leaf area
             ! NOTE: summed as conductances
             canopy_soil_resistance_out(p) = canopy_soil_resistance_out(p) + 1d0 / ( soilR_in(i) / frac )
          end do
       else
          canopy_soil_resistance_out(p) = 0.001d0
       end if ! frac > 0
    end do ! p through canopy
    ! then convert to resistance
    canopy_soil_resistance_out = 1d0 / canopy_soil_resistance_out

    ! explicit return
    return

  end subroutine calculate_canopy_soil_resistance
  !
  !----------------------------------------------------------------------
  !
  subroutine roots(stock_roots_in,rooted_layers_out,root_length_out,root_mass_out)

    ! determines root distribution given root_biomass !

    use gv_scale_declarations, only: pi, grid
    use gv_soil_structure,     only: layer_depth, max_depth, root_biomass, &
                                     root_reach, root_k, surf_biomass, thickness
    use gv_veg,                only: root_radius, root_density
    use math_tools,            only: zbrent

    implicit none

    ! arguments
    integer, intent(out) :: rooted_layers_out
    double precision, intent(in) :: stock_roots_in
    double precision, dimension(grid%core), intent(out) :: root_length_out, &
                                                             root_mass_out 

    ! local variables..
    integer        :: i
    double precision           :: cumdepth, curr, depth, mult, preb, prev, &
                                  root_depth, root_cross_sec_area, slpa, x1, x2, xx(1)
    double precision,parameter :: xacc = 0.0001d0, min_root_biomass = 5d0

    depth = 0d0
    preb  = 0d0

    ! convert from gC to g biomass
    root_biomass = stock_roots_in * 2d0

    ! always provide a minimum root biomass
    root_biomass = max( min_root_biomass , root_biomass )
    root_cross_sec_area = pi * root_radius * root_radius    ! root X-sectional area (m2)
    root_depth = max_depth * root_biomass / ( root_k + root_biomass )    ! rmass=fine root C

    i = size( shape( layer_depth ) )  
    xx = minloc( layer_depth , MASK=layer_depth .gt. root_depth )
    rooted_layers_out = int( xx(1) )
    root_reach = layer_depth( rooted_layers_out )    ! to what soil layer do the roots pentrate?

    ! ensure 50% of root mass is in top 25% of rooted layers 
    ! see pipochronosequence.xls for derivation of this relationship
    mult = min( 1d0 / thickness(1) , max( 2d0 , 11d0 * exp( -0.006d0 * root_biomass ) ) )

    ! estimate the amount of root at the surface as input to the exponential
    ! decay model
    surf_biomass = root_biomass * mult

    ! conduct integration of the exponential equation.
    if (rooted_layers_out > 1) then
        x1 = 0.1d0
        x2 = 10.0d0
        ! determine slope of root distribution given rooting depth 
        ! and ratio of root mass to surface root density
        slpa = zbrent( 'roots:root_dist' , root_dist , x1 , x2 , xacc )
        prev = 1d0 / slpa
        cumdepth = 0d0
        do i = 1 , rooted_layers_out
           cumdepth       = cumdepth + thickness(i)
           curr           = 1d0 / slpa * exp( -slpa * cumdepth )
           ! root_mass is g biomass, i.e. twice the C content
           root_mass_out(i)   = ( prev - curr ) * surf_biomass
           root_length_out(i) = root_mass_out( i ) / ( root_density * root_cross_sec_area )    ! (m m-3 soil)
           prev               = curr
        enddo
    else
        root_mass_out(1)  = root_biomass
        root_length_out(1)= root_biomass / root_density * root_cross_sec_area
    end if

  end subroutine roots
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
  double precision function root_dist( xin )    ! parameters

    use gv_soil_structure

    ! the integral of the exponential decay (eq 1) should be equal to the ratio
    ! of the total stock and the max value in profile (eq 2). Therefore, here we
    ! are solving for the minimum difference. Search integration of exponential
    ! function. NOTE: here for simplicity here we neglect application of log to
    ! both sides of the equation.

    implicit none

    ! arguments..
    double precision, intent(in) :: xin       ! slope parameter

    ! see pipo chronosequence.xls for this relationship
    root_dist = ( 1d0 - exp( -xin * root_reach ) ) / xin ! eq 1
    root_dist = root_dist - root_biomass / surf_biomass  ! eq 2
    
    return

  end function root_dist
  !
  !----------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    use gv_scale_declarations, only: user_opts
    use gv_soil_structure, only: potA, potB

    implicit none

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision ::soil_wp

!    ! calculate the soil water potential (kPa)..
!    soil_WP = -0.001 * potA( water_retention_pass ) * xin**potB( water_retention_pass )
!    water_retention_saxton_eqns = 1000. * soil_wp + 10.    ! 10 kPa represents air-entry swp
    ! calculate the soil water potential (kPa)..
    soil_wp = -potA( water_retention_pass ) * xin**potB( water_retention_pass )
    water_retention_saxton_eqns = soil_wp + 10d0    ! 10 kPa represents air-entry swp

  end function water_retention_saxton_eqns
  !
  !----------------------------------------------------------------------
  !
end module soil_air
!
!------------------------------------------------------------------------
