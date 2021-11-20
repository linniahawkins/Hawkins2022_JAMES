! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_initialise

  !! This module gives initial values to variables at the !!
  !!  start of a 'from-scratch' SPA simulation.           !!
  !! The alternative is to start SPA from some previous   !!
  !!  state, which does not require these steps.          !!

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: initialise_crops, initialise_soils, initialise_veg

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first...
  !
  !----------------------------------------------------------------------
  !
  subroutine initialise_crops

    ! Anything required to initialise crops.. !

    use carbon_model_crop_mod

    implicit none

  end subroutine initialise_crops
  !
  !----------------------------------------------------------------------
  !
  subroutine initialise_soils

    ! Go through the steps necessary to initialise the soils from scratch !

    use gv_hourscale,          only: hourts
    use gv_scale_declarations, only: grid
    use gv_soil_structure,     only: soil_temp, thickness, waterfrac, watericemm, wettingbot, wettingtop, field_capacity &
                                    ,prevwater, potB, potA
    use gv_snow_Info,          only: snowheight, snowweight, snow_watermm
    use soil_air,              only: saxton_parameters, soil_porosity, water_retention

    implicit none

    ! local variables..
    integer :: i
    logical :: saturated

    ! calculate the conductivity parameters
    call saxton_parameters  
    ! calculate field capacity
    call water_retention

    watericemm = 0d0 ; wettingbot = 0d0 ; wettingtop = 0d0
    snow_watermm = 0d0 ; snowheight = 0d0 ; snowweight = 0d0
    hourts = soil_temp(1) ! initial soil temp estimate for light routine
    ! loop through to determine relative saturation
    i = 0 ; saturated = .false.
    do while (saturated .eqv. .false.) 
       ! first where is the top of the saturated layer (well actually where is
       ! it not dry)
       i = i + 1
       wettingtop(i) = thickness(i)*(1d0-(waterfrac(i)/field_capacity(i)))
       if (wettingtop(i) < thickness(i) .or. i == (grid%wetting-1)) saturated = .true.
    end do
    saturated = .false.
    do while (saturated .eqv. .false.) 
       ! then where does the saturated layer end (well actually saturated or one
       ! layer below the wettingtop)
       i = i + 1
       wettingbot(i) = thickness(i)*(1d0-(waterfrac(i)/field_capacity(i)))
       if (wettingbot(i) > 0d0 .or. i == grid%wetting) saturated = .true.
    end do

    ! calculate soil porosity   
    call soil_porosity   

    ! v large mesophyll conductance - ACi curves have not been Epron adjusted
    prevwater  = 0d0
    do i = 1 , grid%soil
      prevwater = prevwater + 1d3 * ( waterfrac(i) * thickness(i) )
    enddo

  end subroutine initialise_soils
  !
  !----------------------------------------------------------------------
  !
  subroutine initialise_veg

    ! Calculate the initial leaf-water potential, which requires !
    ! knowing the soil water potential, which in turn requires   !
    ! knowing about the conductance and resistance of the soil.  !

    use gv_meteo,              only: head
    use gv_scale_declarations, only: grid, user_opts, dble_zero, dble_one
    use gv_soil_structure,     only: conduc, rooted_layers, root_length, root_mass &
                                    ,waterfrac, weighted_SWP, soilR, soilR1, soilR2
    use gv_veg,                only: canopy_height, lai, layer_height, LWPstore, &
                                     LCA, avN, nfrac, lafrac, Nla
    use gv_carbon_model,       only: pars, POOLS
    use soil_air,              only: slayer, soil_conductivity, soil_resistance, &
                                     soil_water_potential, water_uptake_layer

    implicit none

    ! local variables..
    integer :: i
    double precision    :: dummy_output
  
    canopy_height = layer_height(1)

    ! Calculate soil-conductance. This requires the saxton parameters
    ! to have been calculated already (which happens in init_soils)
    do slayer = 1 , grid%core   ! loop over layers..
      conduc(slayer) = soil_conductivity( waterfrac(slayer) )
    enddo

    ! need some initial values, not sure why crops is excluded here, 
    ! possible should test removing this condition and seeing what happens
    if (user_opts%plant_func_type == 2 .or. user_opts%plant_func_type == 3) then
        ! initial lai estimate (made up)
        lai = lafrac * 1d0
        ! estimate initial LCA based on maintaining C:N ratio through canopy
        LCA = 0.d0 ; Nla = (nfrac * avN * sum( lai ))
        where (lai > 0d0) Nla = Nla / lai ! avN per canopy layer (gN.m-2)
        where (lai > 0d0) LCA = pars(17) * (Nla / avN)
        ! now update with actual initial condtions to be consistent with C cycle
        ! models
        lai = lafrac * POOLS(1,2) / pars(17)

        ! And finally used to get the initial leaf-water potential..
        do i = 1 , grid%canopy
           LWPstore(i) = - head * layer_height(i)
        enddo
    else
        ! Need to provide some initial values..
        root_length   = 0.1d0
        root_mass     = 0.1d0
        rooted_layers = grid%core
        ! initial lai estimate (made up)
        lai = lafrac * dble_one 
        ! estimate initial LCA based on maintaining C:N ratio through canopy
        LCA = 0d0 ; Nla = (nfrac * avN * sum( lai )) 
        where (lai > 0d0) Nla = Nla / lai ! avN per canopy layer (gN.m-2)
        where (lai > 0d0) LCA = pars(17) * (Nla / avN)
        ! now update with actual initial condtions to be consistent with C cycle
        ! models (note that POOLS are zeroed initially)
        lai = lafrac * POOLS(1,2) / pars(17)

        ! Calculate resistance to flow (we need soilR)
        call soil_resistance(soilR,soilR1,soilR2,root_length,root_mass,rooted_layers)
        ! and the soil water potential..
        call soil_water_potential
        ! ..use to calculate the weighted soil-water-potential..
        call water_uptake_layer( dummy_output )
        ! And finally used to get the initial leaf-water potential..
        do i = 1 , grid%canopy
           LWPstore(i) = weighted_SWP - head * layer_height(i)
        enddo

    end if ! pft

  end subroutine initialise_veg
  !
  !----------------------------------------------------------------------
  !
end module spa_initialise
!
!------------------------------------------------------------------------
!

