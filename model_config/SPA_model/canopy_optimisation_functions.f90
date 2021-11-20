
module canopy_optimisation_mod

   ! description

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental
  ! overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: calculate_lai_marginal_return, & ! marginal return on additional LAI
            calculate_root_marginal_return

  ! publica variables
  public :: deltaNfrac, deltaLCA

  ! needed variables
  double precision :: Ntop ! N at the very top of the canopy (gN)

  ! some model parameters
  double precision :: deltaNfrac = 0.05d0  & ! fraction of N per canopy layer 
                     ,deltaLCA = 5d0   ! gC.m-2 leaf area

contains

  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_lai_marginal_return(lai_in)

    use gv_scale_declarations, only: grid,time,user_opts,umol_to_g_carbon
    use leaf, only: assimilate, set_leaf
    use light, only: solar_update_for_marginal_return,longem
    use gv_veg, only: Nla,lafrac,LWPprevious,LCA
    use gv_clim, only: leaf_heat_conductance,gbw,ppfd_to_par &
                      ,temp_top,temp_bot,wdbot,wdtop
    use gv_meteo, only: la,layer_LCA,nit,dew_evap,live_frac &
                       ,rad,temp,wdef,gbb,gbh,par,psil

    ! Based on leaf_processes subroutine to calculate the potential GPP (Agross)
    ! (gC.m-2.day-1) if potential leaf growth is given go ahead. 
    ! Currently the radiative transfer is updated and current soil mosture demands, 
    ! however no feedback on water is carried through the marginal calculation over the day.

    implicit none

    ! arguments
    double precision, intent(in) :: lai_in 

    ! local variables
    integer :: i
    double precision, dimension(grid%canopy) :: frac_sun,apar_sun,apar_shade &
                                   ,anir,along,lai_for_solar
    double precision    :: agr,     & !
               dew_rad, & ! 
               gsm,     & !
               modet,   & !
          modet_kg_s,   & !
               res        !

    ! call updated radiative transfer
    apar_sun = 0d0 ; apar_shade = 0d0 ; frac_sun = 0d0
    anir = 0d0 ; along = 0d0 ; lai_for_solar = lai_in * lafrac
    call solar_update_for_marginal_return(lai_for_solar,frac_sun &
                                         ,apar_sun,apar_shade &
                                         ,anir,along)

    ! loop through canopy to estimate GPP return
    calculate_lai_marginal_return = 0d0
    do i = 1, grid%canopy

       ! do sun leaf areas first
       la = frac_sun(i) * (lai_in * lafrac(i))
       if (la > 0d0) then
           nit  = Nla(i)
           layer_LCA = LCA(i)
           psil = LWPprevious(i) 
           par  = ( apar_sun(i) + apar_shade(i) * frac_sun(i) ) / la
           if (user_opts%iterative_canopy) then
              ! calcalate impact on energy balance of dew to this canopy layer
              ! (W.m-2)
              dew_rad = (dew_evap*live_frac*lafrac(i)*frac_sun(i)) / la
           else
              dew_rad = 0.
           end if

           ! net radiation = shortwave + longwave radiation
           rad = 0.001d0 * ( anir(i) * frac_sun(i) / la  &
!                         + along(i) * frac_sun(i) / la &
                         + longem(i) * frac_sun(i) / la &
                         + (par / ppfd_to_par) - dew_rad )
           temp = temp_top - 0.1111d0 * dble( i - 1 ) * ( temp_top - temp_bot )
           wdef = wdtop - 0.1111d0 * dble( i - 1 ) * ( wdtop - wdbot )
           gbb  = gbw(i)                    ! boundary layer for water vapour (m.s-1)
           gbh  = leaf_heat_conductance(i)  ! boundary layer for heat (m.s-1)
           call set_leaf(i,frac_sun(i))    ! see leaf.f90
           call assimilate(i,time,.true.,modet,modet_kg_s,agr,res,gsm,.true.)    ! see leaf.f90
           ! update marginal calculation
           calculate_lai_marginal_return = calculate_lai_marginal_return &
                                         + (agr*time%seconds_per_day*umol_to_g_carbon)
       endif ! sun leaves

       ! do shaded leaf areas second
       la = (1d0-frac_sun(i)) * (lai_in*lafrac(i))
       if (la > 0d0) then
           nit  = Nla(i)
           layer_LCA = LCA(i)
           psil = LWPprevious(i) 
           par  = apar_shade(i) * (1d0-frac_sun(i)) / la
           if (user_opts%iterative_canopy) then
              ! calcalate impact on energy balance of dew to this canopy layer
              ! (W.m-2)
              dew_rad = (dew_evap*live_frac*lafrac(i)*(1d0-frac_sun(i))) / la
           else
              dew_rad = 0d0
           end if
           ! net radiation = shortwave + longwave radiation
           rad = 0.001d0 * ( anir(i) * (1d0-frac_sun(i)) / la  &
                         + longem(i) * (1d0-frac_sun(i)) / la &
!                         + along(i) * (1d0-frac_sun(i)) / la &
                         + (par / ppfd_to_par) - dew_rad )
           temp = temp_top - 0.1111d0 * dble( i - 1 ) * ( temp_top - temp_bot )
           wdef = wdtop - 0.1111d0 * dble( i - 1 ) * ( wdtop - wdbot )
           gbb  = gbw(i)                    ! boundary layer for water vapour (m.s-1)
           gbh  = leaf_heat_conductance(i)  ! boundary layer for heat (m.s-1)
           call set_leaf(i,(1d0-frac_sun(i)))    ! see leaf.f90
           call assimilate(i,time,.false.,modet,modet_kg_s,agr,res,gsm,.true.)    ! see leaf.f90

           ! update marginal calculation (gC.m-2.day-1)
           calculate_lai_marginal_return = calculate_lai_marginal_return &
                                         + (agr*time%seconds_per_day*umol_to_g_carbon)

       endif ! shade leaves

    enddo ! end canopy loop

    ! explicit return command
    return

  end function calculate_lai_marginal_return
  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_root_marginal_return(root_in)

    use gv_scale_declarations, only: grid,time,user_opts,umol_to_g_carbon &
                                    ,dble_one, dble_zero
    use leaf, only: assimilate, set_leaf
    use light, only: solar_update_for_marginal_return,longem
    use gv_veg, only: Nla,lafrac,LWPprevious,LCA,lai,canopy_soil_resistance
    use gv_clim, only: leaf_heat_conductance,gbw,ppfd_to_par &
                      ,temp_top,temp_bot,wdbot,wdtop
    use gv_meteo, only: la,layer_LCA,nit,dew_evap,live_frac &
                       ,rad,temp,wdef,gbb,gbh,par,psil
    use soil_air, only: soil_resistance, calculate_canopy_soil_resistance, roots
    use gv_soil_structure, only: rooted_layers, root_mass, root_length, &
                                 soilR, soilR1, soilR2

    ! Based on leaf_processes subroutine to calculate the potential GPP (Agross)
    ! (gC.m-2.day-1) if potential root growth is given go ahead. 
    ! Currently the radiative transfer is updated and current soil mosture demands, 
    ! however no feedback on water is carried through the marginal calculation over the day.

    implicit none

    ! arguments
    double precision, intent(in) :: root_in 

    ! local variables
    integer :: i, rooted_layers_hold
    double precision, dimension(grid%canopy) :: frac_sun,apar_sun,apar_shade &
                                   ,anir,along,lai_for_solar,canopy_soil_resistance_hold
    double precision, dimension(grid%core) :: root_mass_hold, &
                                            root_length_hold, &
                        soilR_hold, soilR1_hold, soilR2_hold
    double precision    :: agr,     & !
               dew_rad, & ! 
               gsm,     & !
               modet,   & !
          modet_kg_s,   & !
               res        !

    ! call updated radiative transfer
    apar_sun = 0d0 ; apar_shade = 0d0 ; frac_sun = 0d0
    anir = 0d0 ; along = 0d0 ; lai_for_solar = lai
    call solar_update_for_marginal_return(lai_for_solar,frac_sun &
                                         ,apar_sun,apar_shade &
                                         ,anir,along)

    ! hold root hydraulic related variables in memory locally 
    root_mass_hold = root_mass ; root_length_hold = root_length
    soilR_hold = soilR ; soilR1_hold = soilR1 ; soilR2_hold = soilR2 
    rooted_layers_hold = rooted_layers ; canopy_soil_resistance_hold = canopy_soil_resistance

    ! update root hydraulic parameters
    call roots(root_in,rooted_layers,root_length,root_mass)
    call soil_resistance(soilR,soilR1,soilR2,root_length,root_mass,rooted_layers) 
    call calculate_canopy_soil_resistance(canopy_soil_resistance,soilR,lai,rooted_layers)

    ! loop through canopy to estimate GPP return
    calculate_root_marginal_return = 0d0
    do i = 1, grid%canopy

       ! do sun leaf areas first
       la = frac_sun(i) * lai(i)
       if (la > dble_zero) then
           nit  = Nla(i)
           layer_LCA = LCA(i)
           psil = LWPprevious(i) 
           par  = ( apar_sun(i) + apar_shade(i) * frac_sun(i) ) / la
           if (user_opts%iterative_canopy) then
              ! calcalate impact on energy balance of dew to this canopy layer
              ! (W.m-2)
              dew_rad = (dew_evap*live_frac*lafrac(i)*frac_sun(i)) / la
           else
              dew_rad = dble_zero
           end if

           ! net radiation = shortwave + longwave radiation
           rad = 0.001d0 * ( anir(i) * frac_sun(i) / la  &
!                         + along(i) * frac_sun(i) / la &
                         + longem(i) * frac_sun(i) / la &
                         + (par / ppfd_to_par) - dew_rad )
           temp = temp_top - 0.1111d0 * dble( i - 1 ) * ( temp_top - temp_bot )
           wdef = wdtop - 0.1111d0 * dble( i - 1 ) * ( wdtop - wdbot )
           gbb  = gbw(i)                    ! boundary layer for water vapour (m.s-1)
           gbh  = leaf_heat_conductance(i)  ! boundary layer for heat (m.s-1)
           call set_leaf(i,frac_sun(i))    ! see leaf.f90
           call assimilate(i,time,.true.,modet,modet_kg_s,agr,res,gsm,.true.)    ! see leaf.f90
           ! update marginal calculation
           calculate_root_marginal_return = calculate_root_marginal_return &
                                         + (agr*time%seconds_per_day*umol_to_g_carbon)
       endif ! sun leaves

       ! do shaded leaf areas second
       la = (dble_one-frac_sun(i)) * lai(i)
       if (la > 0d0) then
           nit  = Nla(i)
           layer_LCA = LCA(i)
           psil = LWPprevious(i) 
           par  = apar_shade(i) * (dble_one-frac_sun(i)) / la
           if (user_opts%iterative_canopy) then
              ! calcalate impact on energy balance of dew to this canopy layer
              ! (W.m-2)
              dew_rad = (dew_evap*live_frac*lafrac(i)*(dble_one-frac_sun(i))) / la
           else
              dew_rad = dble_zero
           end if
           ! net radiation = shortwave + longwave radiation
           rad = 0.001d0 * ( anir(i) * (dble_one-frac_sun(i)) / la  &
                         + longem(i) * (dble_one-frac_sun(i)) / la &
!                         + along(i) * (dble_one-frac_sun(i)) / la &
                         + (par / ppfd_to_par) - dew_rad )
           temp = temp_top - 0.1111d0 * dble( i - 1 ) * ( temp_top - temp_bot )
           wdef = wdtop - 0.1111d0 * dble( i - 1 ) * ( wdtop - wdbot )
           gbb  = gbw(i)                    ! boundary layer for water vapour (m.s-1)
           gbh  = leaf_heat_conductance(i)  ! boundary layer for heat (m.s-1)
           call set_leaf(i,(dble_one-frac_sun(i)))    ! see leaf.f90
           call assimilate(i,time,.false.,modet,modet_kg_s,agr,res,gsm,.true.)    ! see leaf.f90

           ! update marginal calculation (gC.m-2.day-1)
           calculate_root_marginal_return = calculate_root_marginal_return &
                                         + (agr*time%seconds_per_day*umol_to_g_carbon)

       endif ! shade leaves

    enddo ! end canopy loop

    ! revert root hydraulic related variables back to original values
    root_mass = root_mass_hold ; root_length = root_length_hold
    soilR = soilR_hold ; soilR1 = soilR1_hold ; soilR2 = soilR2_hold
    rooted_layers = rooted_layers_hold ; canopy_soil_resistance = canopy_soil_resistance_hold

    ! explicit return command
    return

  end function calculate_root_marginal_return
  !
  !----------------------------------------------------------------------
  !
!
!----------------------------------------------------------------------
!
end module canopy_optimisation_mod
