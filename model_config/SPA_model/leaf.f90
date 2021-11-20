! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module leaf

  !!           > module descriptor here <             !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA...
  public :: assimilate, gs2, gsx, g1_med, leaf_balance, leaf_temperature, set_leaf, leaf_temp_precision

  ! variables only used within this module..

  double precision, parameter :: Vc_kurtosis = 0.143d0, & ! kurtosis of Vcmax temp response
                                 Vj_kurtosis = 0.172d0, & ! kurtosis of Jmax temp response
                                 leaf_temp_precision = 0.01d0

  ! parameters for Arrhensis adjustments to photosynthetic parameters
  double precision, parameter :: kc_saturation = 310d0,    &
                              kc_half_sat_conc = 23.956d0, &
                                 ko_saturation = 155d0,    &
                              ko_half_sat_conc = 14.509d0, &
                            co2comp_saturation = 36.5d0,   &
                         co2comp_half_sat_conc = 9.46d0

  logical :: c3path ! use C3 or C4 photosynthetic models 

  double precision :: et, & ! Evapotranspiration through stomata,for a given canopy layer.
                            !  Both (g.m-2.s-1) and (mol.m-2.s-1) are used
          gamma1, & ! gamma arrhensis curves for modification of limitation on photosynthesis
             gs2, & ! common block stomata conductance
             gsx, & ! testing medlyn gsx umol/m2/s
        vpd_term, & ! vpd term in medlyn model (kPa)
          g1_med, & ! medlyn g1 parameter
              lt, & ! leaf temperature for photosynthesis (oC)
              kc, & ! half saturation concentration for carboxylation rate (umol.m-2)
              ko, & ! half saturation concentration for oxygenation rate (umol.m-2)
            resp, & ! Canopy layer respiration (umolC.m-2)
           vcmax, & ! Temperature modified vcmax rates
           vjmax    ! Temperature modifed vjmax rates

  save

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first..
  !
  !----------------------------------------------------------------------
  !
  subroutine assimilate(clayer,time,in_sun_flag,etr,etr_kg_s,agr,res,gsm,marginal)

    ! Assimilation is determined through equalising assimilation !
    !  according to the Farquhar model and a diffusion model,    !
    !  through varying internal CO2 concentration (Ci).          !
    ! Rsoil calculations have been moved to soil_air.            !
    ! switch variables.. ded to specify                          !
    !    -conductance vs conductivity                            !
    !    -Farquhar parameters are N-linked or not                !

    use gv_clim,               only: atmos_press
    use gv_metab,              only: an, ci, layer_capac, rn
    use gv_meteo,              only: lwp_pd, gbb, la, layer_LCA, nit, psil, rad, rad_pass, Rcon, temp, wdef
    use gv_scale_declarations, only: boltz_kW, time_holder, user_opts, g_to_mol_water, mol_to_g_water, &
                                     umol_to_g_water
    use gv_veg,                only: co2amb, c3, emiss, flux, psil_pd
    use gv_carbon_model,       only: twq
    use gv_soil_structure,     only: SWP
    use math_tools,            only: zbrent
    use spa_io,                only: handle_output

    implicit none

    integer,intent(in)           :: clayer
    type(time_holder),intent(in) :: time
    logical,intent(in)           :: in_sun_flag, marginal
    double precision,intent(out)             :: agr, etr, etr_kg_s, gsm, res

    ! local variables..
    double precision :: ad, asn, check1, check2, lt_diff, conv, darkresp, dpsildt, g1, g2, gs, &
            lambda, lt, maxg, netrad, prevpsil, xacc, tmp, d13c, gsmed

    ! reset some veg-module variables..
    an = 0d0 ; ci = 0d0 ; et = 0d0
    
    ! initialise local variables..
    ad = 0d0 ; asn = 0d0 ; gs = 0.00001d0 ; xacc = 0.0001d0
    c3path = c3( clayer )  ! determine PS pathway for layer

    if (time%step .eq. time%steps_per_day/6 ) then ! 4 am

        psil_pd(clayer) = SWP(4) ! psil

    end if

    lwp_pd = psil_pd(clayer)

    ! Brent's method to determine gs (m s-1) using Ball-Berry
    ! set maximum and minimum values
    g1 = 0.00005d0 ; g2 = 0.05d0

    if ( user_opts%iWUE .gt. 2 ) then
      !check if conditions allow photosynthesis
      !check1 = stomdiff( g1 )
      !check2 = stomdiff( g2 )
      check1 = -1
      check2 = 1
      maxg   = 0.002d0
      if ( check1*check2 .lt. 0 ) then
        if ( maxg .gt. 0.00002d0 ) then
            lt_diff = zbrent('assimilate:TleafFunc',TleafFunc,temp-20d0,temp+20d0,0.05d0)
        else
            gs = maxg
        end if 

        prevpsil = psil
        call deltapotential( psil )    ! change in leaf water potential
        dpsildt = psil - prevpsil

        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
           lt = leaf_temperature( gs2 )
        elseif (user_opts%solve_canopy_temperature == 2) then
           lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if

        !  determine if radiation is isothermal or net radiation
        if (rad_pass > 1) then
          ! radiation already netrad
          netrad = rad
        else
          ! determine net radiation from isothermal, longwave correction (kW.m-2)
          netrad = rad-4d0 * emiss * boltz_kW * ( temp + 273.15d0 )**3 * ( lt - temp )
        end if

        et = trans( gs2 , lt , netrad , wdef , gbb )
        !et = evap( gs2 , lt , netrad , wdef , gbb )
        res = resp ; agr = ( an + resp )
        d13c = 4.4d0+22.6d0*(ci/co2amb)

      else ! dark
        gs = 0.00004d0
        ci = co2amb

        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
            gs2 = gs
            lt = leaf_temperature( gs )
        elseif (user_opts%solve_canopy_temperature == 2) then
            gs2 = gs
            lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if

        !  determine if radiation is isothermal or net radiation
        if (rad_pass > 1) then
            ! radiation already netrad
            netrad=rad
        else
            ! determine net radiation from isothermal, longwave correction
            ! (kW.m-2)
            netrad=rad - 4d0 * emiss * boltz_kW * ( temp + 273.15d0 )**3 * ( lt - temp )
        end if

        !et = evap( gs , lt , netrad , wdef , gbb )
        et = trans( gs2 , lt , netrad , wdef , gbb )
!       darkresp = Rm_reich(lt,nit,twq,layer_LCA)
        darkresp = rn * nit * exp( log(2d0) * ( lt - 10d0 ) / 10d0 )
        an  = -darkresp
        res = darkresp
        agr = 0.0
        d13c = 4.4d0+22.6d0*(ci/co2amb)
        prevpsil = psil
        call deltapotential( psil )  ! change in leaf water potential
        dpsildt = psil - prevpsil

      end if

    elseif (user_opts%iWUE .lt. 3) then
      check1 = stomdiff( g1 )
      check2 = stomdiff( g2 )
      maxg   = 0.002d0
      g1_med = 0d0

      if ( check1 * check2 .lt. 0 ) then

        if ( maxg .gt. 0.00002d0 ) then 
          gs = zbrent( 'assimilate:stomdiff' , stomdiff , g1 , g2 , xacc )
        else        ! drought limitation - psil has fallen to minpsil (=minlwp)
          gs = maxg
        end if

        prevpsil = psil
        call deltapotential( psil )    ! change in leaf water potential
        dpsildt = psil - prevpsil

        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
           gs2 = gs
           lt = leaf_temperature( gs )
        elseif (user_opts%solve_canopy_temperature == 2) then
           gs2 = gs
           lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if

        !  determine if radiation is isothermal or net radiation
        if (rad_pass > 1) then
          ! radiation already netrad
          netrad = rad
        else
          ! determine net radiation from isothermal, longwave correction
          ! (kW.m-2)
          netrad = rad-4d0 * emiss * boltz_kW * ( temp + 273.15d0 )**3 * ( lt - temp )
        end if

        !et = evap( gs2 , lt , netrad , wdef , gbb )
        et = trans(gs2 , lt , netrad , wdef , gbb ) 
        res = resp ; agr = ( an + resp )
        d13c = 4.4d0+22.6d0*(ci/co2amb) 

      else ! dark
        gs = 0.00004d0
        ci = co2amb

        ! select method for calculating surface energy balance
        if (user_opts%solve_canopy_temperature == 1) then
            gs2 = gs
            lt = leaf_temperature( gs )
        elseif (user_opts%solve_canopy_temperature == 2) then
            gs2 = gs
            lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
        end if

        !  determine if radiation is isothermal or net radiation
        if (rad_pass > 1) then
            ! radiation already netrad
            netrad=rad
        else
            ! determine net radiation from isothermal, longwave correction
            ! (kW.m-2)
            netrad=rad - 4d0 * emiss * boltz_kW * ( temp + 273.15d0 )**3 * ( lt - temp )
        end if

        !et = evap( gs , lt , netrad , wdef , gbb )
        et = trans(gs2 , lt , netrad , wdef , gbb )
!       darkresp = Rm_reich(lt,nit,twq,layer_LCA)
        darkresp = rn * nit * exp( log(2d0) * ( lt - 10d0 ) / 10d0 )
        an  = -darkresp
        res = darkresp
        agr = 0.0
        d13c = 4.4d0+22.6d0*(ci/co2amb)
        prevpsil = psil
        call deltapotential( psil )  ! change in leaf water potential
        dpsildt = psil - prevpsil

      end if
    end if

    lambda = 1000d0*(2501d0 - 2.364d0 * lt) ! (J.kg-1)
    ! apply consistent energy balance restriction to transpiration
    et = max(0d0, et)
    ! the convert from g/m2/s to mol/m2/s
    et  = et * g_to_mol_water

    ! correct leaf respiration, Anet, Agross and ET to leaf area
    res = res * la ; an = an * la ; agr = agr * la
    etr = et * mol_to_g_water * 1d-3 * lambda * la ! convert et (mol m-2 s-1) to latent heat (Wm-2)
    etr_kg_s = etr / lambda ! transpiration but now converted back to kg.m-2.s-1 or mm

    ! water flux at base of trunk..
    flux( time%step , clayer ) = flux( time%step , clayer ) + la * et * 1000d0 &
                                  + layer_capac * dpsildt / time%seconds_per_step
    conv = 1000d0 * atmos_press / ( ( lt + 273.15d0 ) * Rcon )    !convert from m s-1 to mmol m-2 s-1
    gsm  = gs2 * conv * la

    gsmed = gsx * 1000d0 * la

    ! Produce output..
    if (user_opts%canopy_csv_output .and. .not.marginal) then
        if ( in_sun_flag ) then ! sun calculations
          call handle_output( 3 , time , output_data = (/ dble(clayer), gsm, agr, &
                                                           res, etr, lt, conv, nit, layer_LCA, d13c, gsmed, vpd_term, g1_med, lwp_pd /) )
        else
          call handle_output( 4 , time , output_data = (/ dble(clayer), gsm, agr, &
                                                           res, etr, lt, conv, nit, layer_LCA, d13c, gsmed, vpd_term, g1_med, lwp_pd /) )
        end if ! sun / shade    layer_LCA
    endif ! generate canopy output

  end subroutine assimilate
  !
  !----------------------------------------------------------------------
  !
  double precision function leaf_balance( leaf_temp_in )
    ! determines leaf temperature by balancing leaf energy balance.
    ! Achieved through bisection proceedure

    use gv_clim,               only: cp_air
    use gv_meteo,              only: abs_pot_conv, air_density_kg, gbb, gbh, temp, rad,rad_pass, wdef
    use gv_scale_declarations, only: boltz
    use gv_veg,                only: emiss

    implicit none

    double precision, intent(in) :: leaf_temp_in
    double precision net, LH, HFX, lambda, rho

    ! calculate netrad (W.m-2)
    if (rad_pass == 1 ) then
        net = ( rad * 1d3 ) - 4d0 * emiss * boltz * ( temp + 273.15d0 )**3 * ( leaf_temp_in - temp )
    else
        net = rad * 1d3
    end if

    ! calculate latent heat of vaporisation (J.kg-1)
    lambda  = ( 2501000d0 - 2364d0 * leaf_temp_in ) 
    ! load air density (kg.m-3)
    rho = air_density_kg

    ! calculate evaporation converted g.m-2.s-1 -> kg.m-2.s-1
    LH = evap( gs2 , leaf_temp_in , (net*1d-3) , wdef , gbb ) * 1d-3 
    ! convert kg.m-2.s-1 -> W.m-2.
    ! restrict to near positive always values to prevent un-realistic
    ! transpiration values 
!    LH = max(-0.5d0, LH*lambda)
    LH = LH * lambda

    ! calculate sensible heat (W.m-2)
    HFX = (gbh*2d0)*cp_air*rho*(((leaf_temp_in+273.15d0)*abs_pot_conv)-((temp+273.15d0)*abs_pot_conv))

    ! calculate leaf level energy balance closure (W/m2)
    leaf_balance=net-LH-HFX

  end function leaf_balance
  !
  !----------------------------------------------------------------------
  !
  double precision function leaf_temperature( gs )

    ! determines leaf temperature. !

    use gv_clim,               only: cp_air
    use gv_meteo,              only: air_density_kg, gbb, gbh, rad, temp, wdef
    use gv_scale_declarations, only: boltz_kW
    use gv_veg,                only: emiss
    use log_tools

    implicit none

    ! arguments..
    double precision,intent(in) :: gs

    ! local variables..
    double precision :: de, denom, diff, ghr, gr, lambda, psych, q, rho, rhr, rt, s, slope, ta, tt, cp

    tt    = temp
    q     = rad
    ta    = tt + 273.15d0
    rho   = air_density_kg * 1d3   ! density of air g m-3 (t-dependent)
    cp    = cp_air * 1d-6 ! convert to correct units

    ! slope of saturation vapour pressure curve (t-dependent)..
    s     = 6.1078d0 * 17.269d0 * 237.3d0 * exp( 17.269d0 * tt / ( 237.3d0 + tt ) )
    slope = 0.1d0 * ( s / ( 237.3d0 + tt )**2 )
    psych = 0.1d0 * ( 0.646 * exp( 0.00097 * tt ) )       ! psych is temp-dependent 
    lambda = 0.001d0 * ( 2501d0 - 2.364d0 * tt )            ! latent heat of vapourisation (KJ g-1)

   ! convert water deficit to vapour pressure deficit..
    de    = wdef * lambda * psych / ( rho * cp )
    gr    = 4d0 * emiss * boltz_kW * ta**3 / ( rho * cp )   ! remove leaf area sensitivity from gr
    ghr   = gr + 2d0 * gbh                               ! total thermal conductance
    rhr   = 1d0 / ghr
    rt    = 1d0 / gs + 1d0 / gbb                        ! combined leaf resistance
    denom = psych * rt + slope * rhr
    diff  = rhr * rt * psych * q / ( rho * cp * denom ) - rhr * de / denom  ! temperature difference

    leaf_temperature = diff + tt          ! leaf temp

  end function leaf_temperature
  !
  !---------------------------------------------------------------------
  !
  subroutine set_leaf( clayer , frac )

    ! sets Farquhar parameters and hydraulic parameters for each leaf layer

    use gv_metab, only: ht, layer_capac, rplant, rsoil, vcm, vjm, rn
    use gv_meteo, only: la, nit
    use gv_soil_structure, only: weighted_SWP,SWP
    use gv_veg,   only: canopy_soil_resistance, capac, conductivity, &
                        gplant, kappac, kappaj, layer_height

    implicit none

    ! arguments..
    double precision:: fgplant
    integer, intent(in) :: clayer
    double precision,    intent(in) :: frac

    ! local variables..
    real, parameter :: leaf_gpp_to_rleaf = 0.105d0 
    real, parameter :: propresp = 1d0 ! proportional respiration constant based on N content
                                      !  and assumed temperature base of 10 oC.

    ! calculate Vcmax given current nitrogen (umolC.m-2.s-1)
    vcm = kappac * nit
    ! calculate Jmax given current nitrogen (umolC.m-2.s-1) 
    vjm = kappaj * nit 
    ! coefficient for Ryan (1991) Leaf maintenance respiration model
    rn = leaf_gpp_to_rleaf * propresp  ! respiration constant umol CO2/g N at 10 deg C
    rn = rn / la           ! convert to resp per m2 

    ! extract height of the current canopy layer
    ht = layer_height( clayer )

    ! plant hydraulic resistances are determined by the amount of leaf area in
    ! the sunlit or shaded fraction
    if ( conductivity .eq. 1 ) then
       !fgplant = gplant
       fgplant = gplant*(0.2d0 + 0.8d0/(1d0+exp(-(SWP(4)+0.784d0)/0.163d0)))
       rplant = ht / ( fgplant * la )  ! MPa s mmol-1 (per layer)
    else
       !fgplant = gplant
       fgplant = gplant*(0.2d0 + 0.8d0/(1d0+exp(-(SWP(4)+0.784d0)/0.163d0)))
       rplant = 1d0 / ( fgplant * la )  ! conductance is constant with height
    end if
    ! calculate canopy layer specific water capacitance
    layer_capac = capac * la
    ! load soil+root resistance of each canopy layer
    rsoil = canopy_soil_resistance( clayer )
    ! now correct rsoil according to the fraction of total leaf area in this run
    rsoil = rsoil / frac

  end subroutine set_leaf
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others.  To maximise precision, the !
    ! calculations have been split & d0 has been used. !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    ! local variables..
    double precision :: answer, denominator, numerator

    numerator   = t - 25d0
    denominator = t + 273.15d0
    answer      = a * exp( b * 1d0 * numerator / denominator )
    arrhenious  = answer

  end function arrhenious
  !
  ! ---------------------------------------------------------------------
  !
  double precision function cdiff( xmid )

    ! difference between metabolic assimilation !
    ! rate and diffusive assimilation rate      !

    use gv_meteo, only: gbb, par 

    implicit none

    ! arguments..
    double precision,intent(in) :: xmid

    ! local variables..
    double precision :: adx, anx

    ! metabolic based photosynthesis umolCO2.m-2.s-1
    if ( c3path ) then
        anx = farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , xmid )
    else
        anx = collatz( vcmax , resp , par , xmid )
    end if

    ! diffusion rate based photosynthesis umolCO2.m-2.s-1 
    adx = diffusion( gs2 , xmid , et , gbb , lt )

    ! difference between metabolic and diffusion based photosynthesis.
    cdiff = adx - max( anx , 0d0 )

  end function cdiff
  !
  !----------------------------------------------------------------------
  !
  double precision function collatz( vcmax , resp , par , ci )

    !  metabolic C4 assimilation rate (umol.C.m-2 ground area s-1)!

    implicit none

    ! arguments..
    double precision, intent (in) :: ci , par , resp , vcmax

    ! local variables..
    double precision :: alpharf, an, bee, beta, k, theta, v1, v2

    alpharf = 0.067d0  ! mol/mol
    k       = 0.7d0    ! mol/m2/s
    theta   = 0.83d0   ! curvature of quantum response curve, Collatz table2
    beta    = 0.93d0   ! collatz table 2
    bee     = alpharf * par + vcmax
    ! Rubisco and light limited capacity..
    v1      = ( bee - sqrt( bee**2 - 4d0 * theta * alpharf * par * vcmax ) ) / ( 2d0 * theta )
    bee     = v1 + k * ci
    ! v1 and CO2 limitation..
    v2      = ( bee - sqrt( bee**2 - 4d0 * beta * v1 * k * ci ) ) / ( 2d0 * beta )
    an      = v2 - resp
    collatz = an

  end function collatz
  !
  !----------------------------------------------------------------------
  !
  double precision function comp( ccc , g1 )

    ! determine light-limitation to photosynthesis !

    use gv_metab,              only: metabolic_opt_temp, rn, vcm, vjm, &
                                     vcmax_max_temp, jmax_max_temp
    use gv_meteo,              only: nit, layer_LCA, par, temp
    use gv_scale_declarations, only: user_opts
    use gv_carbon_model,       only: twq
    use math_tools,            only: zbrent 

    implicit none

    ! arguments..
    double precision,intent(in) :: ccc, g1

    ! local variables..
    double precision            :: gamma1, kc, ko, lt, resx, vcmax, vcmt, vjmax, vjmt

    !  {temperature modifications begin below: based on air - not leaf - temperature}
    gs2 = g1
    ! select method for calculating surface energy balance
    if (user_opts%solve_canopy_temperature == 1) then
       lt = leaf_temperature( gs2 )
    elseif (user_opts%solve_canopy_temperature == 2) then
       lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
    end if

    ! apply non-gaussian temperature modification on carboxylation and
    ! electron-transport rates 
    vcmt  = tempmet(vcmax_max_temp,metabolic_opt_temp,Vc_kurtosis,lt)    ! {temperature modifications begin below}
    vjmt  = tempmet(jmax_max_temp,metabolic_opt_temp,Vj_kurtosis,lt)
    vcmax  = vcmt * vcm
    vjmax  = vjmt * vjm
    ! Temperature adjustments for Michaelis-Menten coefficients 
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    kc     = arrhenious(kc_saturation,kc_half_sat_conc,lt)
    ko     = arrhenious(ko_saturation,ko_half_sat_conc,lt)

    gamma1 = arrhenious(co2comp_saturation,co2comp_half_sat_conc,lt)
!    resx = Rm_reich(lt,nit,twq,layer_LCA)
    resx = rn * nit * exp( log(2d0) * ( lt - 10d0 ) / 10d0 )

    if ( c3path ) then
      comp = farquhar( vcmax , vjmax , kc , ko , gamma1 , resx , par , ccc )
    else
      comp = collatz( vcmax , resx , par , ccc )
    end if

  end function comp
  !
  !----------------------------------------------------------------------
  !
  subroutine deltapotential( leafwp )

    ! change in leaf water potential !

    use math_tools, only: dxsav, kmax, ode_int

    implicit none

    ! arguments..
    double precision,intent(inout)  :: leafwp

    ! local variables..
    integer, parameter :: nvar = 1
    integer            :: nbad, nok
    double precision   :: eps, h1, hmin, x1, x2, ystart(nvar)

    eps  = 1.0d-4
    h1   = 0.01d0
    hmin = 0d0
    kmax = 100
    x1   = 1d0
    x2   = 2d0

    dxsav     = ( x2 - x1 ) * 0.05d0 !/ 20.0
    ystart(1) = leafwp    ! initial conditions

    call ode_int( 'deltapotential:lwp_diff_eqn' , ystart , nvar , x1 , x2 , eps , h1 , hmin , nok , nbad , lwp_diff_eqn )

    leafwp = ystart( 1 )

  end subroutine deltapotential
  !
  !----------------------------------------------------------------------
  !
  double precision function diffusion( gs , ci , e , gbbb , ttemp )

    ! diffusion limited assimilation rate (umol.C.m-2 ground area s-1)!

    use gv_clim,  only: atmos_press
    use gv_meteo, only: gi, Rcon
    use gv_veg,   only: co2amb

    implicit none

    ! arguments..
    double precision,intent(in) :: ci , gs , e , gbbb , ttemp

    ! local variables..
    double precision            :: convert , gt

    diffusion = 0d0

    ! see Jones Appendix 3 for conversion ms-1 to mol s-1)
    ! and Jones Appendix 2 for ratio CO2 diffusion to water diffusion through
    ! stomata (=1.65)
    ! see Jones Appendix 2 for ratio of CO2 diffusion to water diffusion through
    ! boundary layer (=1.37)
    convert = atmos_press / ( Rcon * ( 273.15d0 + ttemp ) )

    ! total leaf conductance (converts from ms-1 to mols-1)
    gt = convert / ( 1.65d0 / gs + 1.37d0 / gbbb + 1d0 / gi )

    ! interaction between water flux out of and CO2 flux into the stomata
    ! determines diffusion limited assimilation rate. Note that mol.m-2.s-1 of
    ! the conductance is converted to umol.m-2.s-1 due to the co2amb and ci
    ! being in ppm or umol/mol
    diffusion = ( gt - 0.5d0 * e ) * co2amb - ( gt + 0.5d0 * e ) * ci

  end function diffusion
  !
  !----------------------------------------------------------------------
  !
  double precision function evap( gs , tt , q , wdef , gbb )

    !  determine evapotranspiration rate (g m-2 s-1) !
    !  from q (kW m-2), tt (oC), wdef (g m-3) and    !
    !  gbb and gs (m s-1)                            !

    implicit none

    ! arguments..
    double precision,intent(in) :: gbb & ! canopy level boundary conductance for water vapour (m.s-1)
                                  ,gs  & ! incoming stomatal conductance (m.s-1)
                                  ,q   & ! canopy level net radiation (kW.m-2)
                                  ,tt  & ! incoming leaf temperature (oC)
                                  ,wdef  ! water deficit of air (kg.m-3)

    ! local variables..
    double precision :: eps    & ! coefficient
                       ,lambda & ! latent heat of vapourisation ((KJ.g-1))
                       ,psych  & ! Pyschrometer constant (kPa.K-1); coefficient between air
                                 ! temperature and wet bulb temp, therefore determining potential rate of latent
                                 ! energy exchange
                       ,s      & ! Straight line approximation of the true slope (below) with temperature
                       ,slope    ! rate of change of saturation vapour pressure with temperature
 
!    gcut   = 0.0005
!    gleaf  = gcut + gs
    ! slope of saturation vapour pressure curve (t-dependent)
    s      = 6.1078d0 * 17.269d0 * 237.3d0 * exp( 17.269d0 * tt / ( 237.3d0 + tt ) )
    slope  = 0.1d0 * ( s / ( 237.3d0 + tt )**2 )      ! (kPa K-1)
    psych  = 0.0646d0 * exp( 0.00097d0 * tt )         ! psych is temp-dependent (kPa K-1)
    eps    = slope / psych                            ! response of epsilon to temp
    lambda = ( 2.5010d0 - 0.002364d0 * tt )           ! latent heat of vapourisation (KJ g-1)
    evap   = ( eps * q / lambda + wdef * gbb ) / ( 1d0 + eps + gbb / gs ) ! (g m-2 s-1)

  end function evap

  !
  !----------------------------------------------------------------------
  !
  double precision function trans( gs , tt , q , wdef , gbb )

    !  determine evapotranspiration rate (g m-2 s-1) !
    !  from q (kW m-2), tt (oC), wdef (g m-3) and    !
    !  gbb and gs (m s-1)              

    use gv_scale_declarations, only: mol_to_g_water  !
    use gv_clim,               only: atmos_press, cp_air !
    use gv_meteo,              only: Rcon, air_density_kg !

    implicit none

    ! arguments..
    double precision,intent(in) :: gbb & ! canopy level boundary conductance for water vapour (m.s-1)
                                  ,gs  & ! incoming stomatal conductance (m.s-1)
                                  ,q   & ! canopy level net radiation (kW.m-2)
                                  ,tt  & ! incoming leaf temperature (oC)
                                  ,wdef  ! water deficit of air (kg.m-3)

    ! local variables..
    double precision :: eps    & ! coefficient
                       ,lambda & ! latent heat of vapourisation ((KJ.g-1))
                       ,psych  & ! Pyschrometer constant (kPa.K-1); coefficient between air
                                 ! temperature and wet bulb temp, therefore
                                 ! determining potential rate of latent
                                 ! energy exchange
                       ,s      & ! Straight line approximation of the true slope (below) with temperature
                       ,slope  & ! rate of change of saturation vapour pressure with temperature
                       ,rho    & ! density of air
                       ,cp     & ! specific heat capacity of air
                       ,vpd      ! vapor presure deficit (mol/mol)
 
    ! local variables..
    double precision :: convert, gsx , de 

    ! Convert gs (m/s) to  (mol/m2/s)
    convert = ( Rcon * ( 273.15d0 + tt ) ) / atmos_press
    gsx = gs / convert

    ! slope of saturation vapour pressure curve (t-dependent)
    s      = 6.1078d0 * 17.269d0 * 237.3d0 * exp( 17.269d0 * tt / ( 237.3d0 + tt ) )
    slope  = 0.1d0 * ( s / ( 237.3d0 + tt )**2 )      ! (kPa K-1)
    ! ---------------- relative humidity from mixing ratio ------------
    rho = air_density_kg * 1d3  ! density of air g/m3 (t-dependent)
    cp = cp_air * 1d-6 ! convert to correct units
    psych = 0.1d0 * ( 0.646 * exp( 0.00097 * tt ) )   ! phsychometric constant is temperature dependent
    lambda = 0.001d0 * ( 2501d0 - 2.364d0 * tt )  ! latent heat of vapor

    de = (wdef * lambda * psych / (rho * cp ) )*1000d0  ! VPD (Pa)
    vpd = de/atmos_press ! convert to mol/mol
    trans = vpd * gsx * mol_to_g_water ! (gm-2 s-1)

  end function trans

  !
  !----------------------------------------------------------------------
  !
  function farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , ci )

    !  metabolic C3 assimilation rate (umol.C.m-2 ground area s-1)!

    implicit none

    ! arguments..
    double precision,intent(in) :: ci, & ! CO2 conc at site of reaction
                               gamma1, & ! CO2 compensation points
                                   kc, & ! Michaelis-Menten constants for CO2
                                   ko, & ! and Oxygen
                                  par, & 
                                 resp, &
                                vcmax, &
                                vjmax

    ! function result..
    double precision :: farquhar

    ! local variables..
    double precision :: alphaj, bee, oi, theta, vc, vj, wc, wj

    oi       = 210.0d0 ! O2 partial pressure umol/mol
    alphaj   = 0.385d0 ! initial slope of quantum response curve
    theta    = 0.7d0   ! curvature of quantum response curve
    farquhar = 0d0
    ! determine Rubisco limited carboxylation rate..
    wc       = ( vcmax * ci ) / ( ci + kc * ( 1d0 + oi / ko ) )
    ! determine potential rate of RuBP regeneration..
    bee      = alphaj * par + vjmax
    vj       = ( bee - sqrt( bee**2 - 4d0 * theta * alphaj * par * vjmax ) ) / ( 2d0 * theta )
    ! determine RuBP regeneration limited carboxylation rate..
    wj       = vj * ci / ( 4.5d0 * ci + 10.5d0 * gamma1 )
    ! determine limiting rate, carboxylation or regeneration
    vc       = min( wc , wj )
    ! net photosynthetic rate, function of CO2 limitation and canopy autotrophic
    ! respiraiton
    farquhar = ( vc * ( 1d0 - gamma1 / ci ) - resp )

  end function farquhar
  !
  !----------------------------------------------------------------------
  !
  double precision function leaf_compensation_point( gs )

    ! leaf compensation point !

    use gv_metab, only: an, ci

    implicit none

    ! arguments..
    double precision,intent(in) :: gs

    ! local variables..
    double precision :: dummy  ! intent-out from stomata, but then never used

    call stomata( gs , ci , dummy )

    leaf_compensation_point = an

  end function leaf_compensation_point
  !
  !----------------------------------------------------------------------
  !
  subroutine lwp_diff_eqn( time_dummy, y, dydt )

    ! differential equation describing change in !
    ! leaf water potential given supply & demand !
    ! See Williams et al 1996 eqn A7 for description

    use gv_metab,              only: ht, layer_capac, rplant, rsoil
    use gv_meteo,              only: la, head, psis
    use gv_scale_declarations, only: max_nos_iterations, time

    implicit none

    ! arguments..
    double precision,intent(in)    :: y(max_nos_iterations)
    double precision,intent(in)    :: time_dummy ! dummy argument, provided for ode_int
    double precision,intent(out)   :: dydt(max_nos_iterations)
 
    dydt(1) = time%seconds_per_step * ( psis - ( head * ht ) - 1000d0 * la * et * &
         ( rplant + rsoil ) - y(1) ) / ( layer_capac * ( rplant + rsoil ) )

  end subroutine lwp_diff_eqn
  !
  !----------------------------------------------------------------------
  !
  subroutine minstom( lowg , maxg )

    ! minimum gs for net C fixation, !
    ! max gs in drought situation.   !

    use gv_metab,      only: an, ci
    use gv_meteo,      only: psil, temp
    use gv_veg,        only: minlwp
    use math_tools,    only: zbrent

    implicit none

    ! arguments..
    double precision,intent(out) :: lowg, maxg

    ! local variables..
    double precision :: low, lwp, high, mings, xacc

    an   = 0d0
    low  = 0.00002d0
    high = 0.05d0
    maxg = high
    xacc = 0.00000001d0
    lwp  = psil
    call stomata( low, ci, lwp )
    if ( (an .gt. 0d0) .or. (temp .lt. 5d0) ) then
       mings = 0.00005d0
    else
       mings = zbrent( 'minstom:leaf_compensation_point' , leaf_compensation_point, low, high, xacc )  !add delta
       mings = mings + 0.00004d0
    end if
    ! check lwp limit
    lwp = psil
    call stomata( mings, ci, lwp )
    if ( lwp .lt. minlwp ) then
       maxg = low
    end if

    lowg = mings

  end subroutine minstom
  !
  !----------------------------------------------------------------------
  !
!  double precision function Rm_reich(leaf_temperature,Narea,airt_warmest_quarter,LCA_in)
!  
!    ! Maintenance respiration (umolC.m-2.s-1) calculated based on modified
!    ! version of the Reich et al (2008) calculation. Temperature modification,
!    ! the mean temperate of the warmest quarter, introduced by Quinn et al.,
!    ! (submitted) applies to the N scaler modifying the baseline metabolic costs
!
!    ! arguments
!    double precision, intent(in) :: leaf_temperature, & ! input temperature of metabolising tissue (oC)
!                                               Narea, & ! N content gN.m-2 leaf area
!                                              LCA_in, & ! leaf C per unit leaf area (gC.m-2)
!                                airt_warmest_quarter    ! mean air temperature (oC) of the warmest quarter of year
!
!    ! local variables
!    double precision, parameter :: Q10 = 2d0,    & ! Q10 response of temperature (baseline = 20oC)
!                          Q10_baseline = 20d0,   & ! Baseline temperature for Q10
!                 N_exponetial_response = 1.411d0,  & ! 
!                    N_scaler_intercept = 1.0246d0, & !
!                     N_scaler_gradient = 0.0359d0, & !
!                           N_g_to_mmol = (1d0/14d0)*1d3    ! i.e. 14 = atomic weight of N
!    double precision :: LMA, N_scaler, Q10_adjustment, Nconc ! Nconc = mmol g-1
!
!    ! convert LCA (gC.m-2) -> LMA (gBiomass.m-2)
!    LMA = LCA_in * 2d0
!
!    ! calculate N concentration as function of biomass:N
!    Nconc = (Narea / LMA) * N_g_to_mmol
!
!    ! calculate instantaneous Q10 temperature response
!    Q10_adjustment = Q10**((leaf_temperature-Q10_baseline)*0.1d0)
!
!    ! calculate temperature sensitive N_scaler
!    N_scaler = N_scaler_intercept - N_scaler_gradient * airt_warmest_quarter
!    ! calculate leaf maintenance respiration (nmolC.g-1.s-1).
!    ! NOTE: that the coefficients in Reich et al., 2008 were calculated from
!    ! log10 linearised version of the model, thus N_scaler is already in log10()
!    ! scale. To remove the need of applying log10(Nconc) and 10**Rm_reich the
!    ! scaler is reverted instead to the correct scale for the exponential form
!    ! of the equations.
!    Rm_reich = Q10_adjustment * (10d0**N_scaler) * Nconc ** N_exponetial_response
!    ! convert nmolC.g-1.s-1 to umolC.m-2.s-1
!    Rm_reich = Rm_reich*1d-3*LMA
!
!    ! explicit return command
!    return
!
!  end function Rm_reich
  !
  !----------------------------------------------------------------------
  !
  subroutine stomata( gs , ci_dummy , lwp )

    ! determines stable ci for given gs !

    use gv_metab,              only: an, ci, metabolic_opt_temp, vcm, vjm, &
                                     vcmax_max_temp, jmax_max_temp, rn 
    use gv_meteo,              only: gbb, layer_LCA, nit, par, rad, temp, wdef, rad_pass
    use gv_scale_declarations, only: boltz_kW, user_opts, g_to_mol_water, pi, &
                                     time
    use gv_veg,                only: co2amb, emiss, latitude_radians
    use gv_soil_structure, only: weighted_SWP
    use gv_carbon_model,       only: twq, pars
    use math_tools,            only: zbrent,calculate_daylength_hours

    implicit none

    ! arguments..
    double precision,intent(in)  :: gs
    double precision,intent(out) :: ci_dummy, lwp  ! had to change ci to ci_dummy as ci is declared in metab

    ! local variables..
    double precision :: ad, netrad, vcmt, vjmt, x1, x2, xacc
    double precision :: max_dayl, max_declin, aoa, daylength, dayl_factor

    gs2 = gs            ! swap argument gs into common block gs2

    ! select method for calculating surface energy balance
    if (user_opts%solve_canopy_temperature == 1) then
       lt = leaf_temperature( gs )
    elseif (user_opts%solve_canopy_temperature == 2) then
       lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
    end if

    !  determine if radiation is isothermal or net radiation
    if (rad_pass > 1) then
        ! radiation already netrad
        netrad = rad
    else
        ! determine net radiation from isothermal, longwave correction
        ! (kW.m-2) based on temperature difference (Jones, p108)
        netrad = rad - 4d0 * emiss * boltz_kW * ( temp + 273.15d0 )**3 * ( lt - temp )
    end if 

    ! evaporation rate in g m-2 s-1, converted to clayer total mol m-2 s-1...
    et = max(0d0,evap( gs , lt , netrad , wdef , gbb ) * g_to_mol_water)
    !et = max(0d0,trans(gs , lt , netrad , wdef , gbb ) * g_to_mol_water)
    
    vcmt  = tempmet(vcmax_max_temp,metabolic_opt_temp,Vc_kurtosis,lt)    !  {temperature modifications begin below}
    vjmt  = tempmet(jmax_max_temp,metabolic_opt_temp,Vj_kurtosis,lt)
    vcmax = vcmt * vcm
    vjmax = vjmt * vjm

    ! Day length modification to Vcmax Ruehr et al., 2014
    !max_declin=0.409571 !23.44
    !aoa=(sin(latitude_radians)*sin(max_declin))/(cos(latitude_radians)*cos(max_declin))
    !max_dayl = 12d0*(1d0+2d0*asin(aoa)/pi)
    !daylength=calculate_daylength_hours(dble(time%day),latitude_radians)
    !dayl_factor=min(1d0,max(0.01,(daylength*daylength)/(max_dayl*max_dayl)))
    !vcmax = vcmax*dayl_factor
    !vjmax = vjmax*dayl_factor

    ! Temperature adjustments for Michaelis-Menten coefficients 
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    kc     = arrhenious(kc_saturation,kc_half_sat_conc,lt)
    ko     = arrhenious(ko_saturation,ko_half_sat_conc,lt)
    gamma1 = arrhenious(co2comp_saturation,co2comp_half_sat_conc,lt)

    ! leaf respiration (umolC.m-2.s-1)
!    resp = Rm_reich(lt,nit,twq,layer_LCA)
    resp = rn * nit * exp( log(2d0) * ( lt - 10d0 ) / 10d0 )

    x1   = co2amb
    xacc = 0.01d0
    
    if (c3path) then !c3
        x2 = 1d0
    else
        x2 = 0d0     
    end if

    ! determine internal CO2 concentration (ppm)
    ci = zbrent( 'stomata:cdiff' , cdiff , x1 , x2 , xacc )

    ! NOW UPDATE VARIABLES OUTSIDE OF THIS MODULE...
    ! all fluxes are umol.C.m-2 ground area s-1
    if ( c3path ) then
        an = farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , ci )   ! stable ci C uptake
    else
        an = collatz( vcmax , resp , par , ci )
    end if
    ad   = diffusion( gs , ci , et , gbb , lt)

    ! change in leaf water potential
    call deltapotential( lwp )

    ci_dummy = ci

  end subroutine stomata
  !
  !----------------------------------------------------------------------
  !
  double precision function stomdiff( gs )

    ! Determines whether the current increment in stomatal conductance 
    ! generates an increase in photosynthesis greater than the minimum
    ! efficiency parameter or results in leaf water potential breaching minimum

    use gv_scale_declarations, only: user_opts
    use gv_metab, only: an 
    use gv_meteo, only: psil, Rcon, temp
    use gv_clim,  only: atmos_press, vpd_top
    use gv_veg,   only: iWUE, WUE, minlwp 

    implicit none

    ! arguments..
    double precision,intent(in) :: gs

    ! local variables..
    double precision            :: an1, an2, ci1, ci2, et1, et2, delta, eff, lwp, minpsi

    ! estimates gs increment (m.s-1) equivalent to 0.001 mol.m-2.s-1 (i.e. 1 mmol.m-2.s-1)
    delta = 0.001d0 / ( atmos_press / ((temp+273.15d0)*Rcon) )
    delta = min(gs,delta)

    ! estimate net photosynthesis (An) and transpiration (et) for current gs and
    ! gs - 1 mmol.m-2.s-1 
    lwp   = psil                       ! hold initial psil value
    call stomata( gs-delta, ci2, lwp )
    et2 = et                           ! transpiration at lower gs (mol.m-2.s-1)
    an2 = an                           ! photosynthesis at lower gs
    lwp = psil                         ! reinitialise psil value
    call stomata( gs, ci1, lwp )
    et1 = et-et2                       ! transpiration at gs (molH2O.m-2.s-1)
    an1 = an-an2                       ! photosynthesis at current gs (umolC.m-2.s-1)
    ! use SPA original intrinsic WUE (dA/dgs) or WUE (dA/dE) efficiency (Williams et al 1996; Bonan et al 2014)
    if (user_opts%iWUE == 1) then
        eff = an1 - iWUE      ! is increase in photosynthesis due to gs increment
                              ! greater than efficiency parameter
                              ! (umolCO2.mmol-1H2O.m-2.s-1) in stomatal
                              ! conductance
    else

        ! checks to prevent 
        if (an1 == 0d0 .and. et1 == 0d0) then
            ! NaN error
            eff = 1d0 - WUE
        else if (an1 > 0d0 .and. et1 == 0d0) then
            ! infinity error
            eff = 9999d0 - WUE
        else
            eff = an1/et1 - WUE   ! is increase in photosynthesis due to gs increment 
                                  ! greater than efficiency parameter (umolCO2.mol-1H2O.m-2.s-1)
        endif

    endif ! end if for iWUE or WUE

    ! has leaf water potential falled further than the minimum leaf water
    ! potential. This is a check to prevent cavitation
    minpsi = lwp - minlwp   

    ! determine whether we are most limited by C efficiency or hydraulics
    ! (i.e. which of eff or minpsi is negative / most negative)
    stomdiff = min( eff, minpsi )

    ! explicit return
    return

  end function stomdiff
  !
  !----------------------------------------------------------------------
  !
  double precision function tempmet( max_temp , opt , q , x )

    ! > function summary? < !

    implicit none

    ! arguments..
    double precision,intent(in) :: max_temp, opt, q, x

    ! local variables..
    double precision :: dummy

    if ( x .ge. max_temp ) then
        tempmet = 0d0
    else
        dummy     = ( max_temp - x ) / ( max_temp - opt )
        dummy     = exp( log( dummy ) * q * ( max_temp - opt ) )
        tempmet = dummy * exp( q * ( x - opt ) )
    end if

  end function tempmet


  ! ----------------------------------------------------------------------
  double precision function CiFunc( ci_val ) 

    ! difference between metabolic assimilation !
    ! rate and diffusive assimilation rate      !
    ! using Ball-Berry or Medlyn                !

    use gv_clim, only: atmos_press, cp_air, vpd_top, temp_top, wind_spd
    use gv_meteo, only: lwp_pd, psil, gbb, par, wdef, gi, Rcon, rad_pass, temp, rad, air_density_kg 
    use gv_metab, only: ci
    use gv_veg, only: co2amb, emiss, cond_slope, cond_slope_b, dimen
    use gv_scale_declarations, only: boltz_kW,user_opts,g_to_mol_water,time,grid
    use math_tools, only: quadratic, zbrent
    use gv_soil_structure,      only: weighted_SWP, SWP

    implicit none

    ! ARGUMENTs..
    double precision,intent(in) :: ci_val

    ! local variables..
    double precision :: gleaf, ci_new 
    double precision :: adx, anx, g0, g1, cs, gbv, rh 
    double precision :: netrad, gs, gs1, convert 
    double precision :: aquad, bquad, cquad, r1, r2 
    double precision :: term 
    !double precision :: vpd_term
    double precision :: cinew, de, psych, rho 
    double precision :: cp, lambda, es, ea 
    double precision :: vpd_min, btran, lwp, a, molCO2_to_molC
    double precision :: gbc, dc, shc, shc_forced, shc_free, scc, gr, grav
    double precision :: visc, shc_lam, shc_turb, b1, re, dleaf

    ! --------- metabolic based photosynthesis umolCO2.m-2.s-1 --------
    if ( c3path ) then
        anx = farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , ci_val )
    else
        anx = collatz( vcmax , resp , par , ci_val )
    end if

    ! ---------------- Calculate boundary layer conductance to CO2 ----
    dleaf = 0.02 ! leaf dimension
    b1 = 1.5d0
    visc = 0.0000148 ! m2/s need to modify for temp and pressure
    dc = 0.000015 ! m2/s need to modify for temp and pressure
    re = wind_spd * dleaf / visc

    scc = visc / dc
    gr = 9.81d0*dleaf**3 * max(lt - temp_top,0d0)/(temp_top*visc*visc) 
    shc_free = 0.54d0 * scc**0.25d0 * gr**0.25d0
    
    shc_turb = b1 * 0.036d0*scc**0.33d0*re**0.8d0
    shc_lam = b1 * 0.66d0 * scc**0.33d0 * re**0.5d0
    shc_forced = max(shc_lam, shc_turb)

    shc = shc_forced + shc_free

    gbc = dc * shc / dleaf

    gbc = gbc*atmos_press/(Rcon*(lt+273.15d0)) ! convert m/s to mol/m2/s
    ! ---------------- co2 at leaf surface ----------------------------    

    ! convert m/s to molH2O/m2/s
    gbv = gbb*atmos_press/(Rcon*(lt+273.15d0))
    gbc = gbv/1.4d0
    gs1 = gs2*atmos_press/(Rcon*(lt+273.15d0))

    cs = co2amb
    cs = co2amb - anx / gbc

    ! ---------------- relative humidity from mixing ratio ------------
    rho = air_density_kg * 1d3  ! density of air g/m3 (t-dependent)
    cp = cp_air * 1d-6 ! convert to correct units
    psych = 0.1d0 * ( 0.646 * exp( 0.00097 * lt ) )   ! phsychometric constant is temperature dependent
    lambda = 0.001d0 * ( 2501d0 - 2.364d0 * lt )  ! latent heat of vapor

    de = wdef * lambda * psych / (rho * cp )  ! VPD (kPa)
    es = .61078d0 * exp( 17.269d0 * lt / ( 237.3d0 + lt) ) ! saturation es (kPa)
    ea = max(0d0,es - de)
    ea = min( max(ea,0.1d0*es),es) ! constrain ea to be <es and >0.05es
    rh = (gbv*ea + gs1*es)/((gbv+gs1)*es)
    
    ! calculate Medlyn VPD terms
    vpd_min = .1d0
    vpd_term = max((es-rh*es),vpd_min) ! kPa
    !print*, rh
    ! --------------- Stomatal constraint function --------------------
    if (user_opts%iWUE == 3) then
        ! Ball-Berry stomatal conductance
        ! Quadratic gs calculation given An. Valid for An>=0

        ! hydraulic limitaiton
        !lwp = psil
        !vcmax = vcmax * ((1+exp(-5.28d0*2.31d0))/(1+exp(-2.31d0-lwp)))

        g0 = 0.01d0
        g1 = cond_slope*exp(cond_slope_b*lwp_pd)
        !g1 = cond_slope

        if (anx .gt. 0d0) then
            aquad = cs
            bquad = cs*(gbv-g0) - (g1 * anx)
            cquad = -gbv * (cs*g0 + g1 * anx * rh)
            call quadratic (aquad,bquad,cquad,r1,r2)
            gsx = max(r1,r2)
        else
            gsx = g0
        end if
    elseif (user_opts%iWUE == 4) then
        ! Medlyn stomatal conductance
        ! Quadratic gs calculation given anx. Valid for anx>=0

        ! hydraulic limitation
        !lwp = psil
        !vcmax = vcmax * ((1+exp(-5.28d0*2.31d0))/(1+exp(-2.31d0-lwp)))

        g0 = 0.01d0
    !    g1 = cond_slope*exp(cond_slope_b*lwp_pd)
        g1 = cond_slope

    !    if (anx .gt. 0d0) then
    !        term = 1.6d0*anx / cs
    !        aquad = 1d0
    !        bquad = -(2d0 * (g0+term) + (g1*term)**2d0 / (gbv*vpd_term))
    !        cquad = g0 * g0 + (2d0*g0 + term * (1d0-g1*g1 / vpd_term))*term
    !        call quadratic (aquad,bquad,cquad,r1,r2)
    !        gsx = max(r1,r2)
    !    else
    !        gsx = g0
    !    end if

    end if
    

    ! medlyn equation
    if (anx .gt. 0d0) then
        gsx = 1.6d0 * (1 + g1/sqrt(vpd_term))*(anx/cs)
    else
        gsx=g0
    end if

    ! Convert gs (mol/m2/s) to (m/s)
    convert = ( Rcon * ( 273.15d0 + lt ) ) / atmos_press
    !g1_med = g1
    gs = gsx*convert
    gs2 = gs    ! swap argument gs into common block gs2

    ! diffusion rate based photosynthesis umolCO2.m-2.s-1 
    !adx = diffusion( gs , ci_val , et , gbb , lt )

    gleaf = 1d0 / (1.4d0/gbv + 1.6d0/gsx)
    !gleaf = gsx/1.6d0
    ci_new = co2amb - anx/gleaf

    ! difference in Ci from metabolic and ci from diffusion photosynthesis.
    CiFunc = ci_new - ci_val
    ci = min(co2amb,ci_new)

    if (anx .lt. 0d0) then
        CiFunc = 0d0
    end if

  end function CiFunc

  !
  !----------------------------------------------------------------------
  !
  subroutine empirical_stomata( tleaf )

    ! determines stable ci for given gs !
    use gv_clim,               only: vpd_top, atmos_press, cp_air
    use gv_metab,              only: an, ci, metabolic_opt_temp, vcm, vjm, &
                                     vcmax_max_temp, jmax_max_temp, rn 
    use gv_meteo,              only: gbb, layer_LCA, nit, par, psil, rad, temp, wdef, rad_pass, Rcon, air_density_kg, lwp_pd
    use gv_scale_declarations, only: boltz_kW, user_opts,g_to_mol_water,time,pi
    use gv_veg,                only: co2amb, emiss, latitude_radians, cond_slope, cond_slope_b, cond_slope_c
    use gv_carbon_model,       only: twq, pars
    use math_tools,            only: zbrent, calculate_daylength_hours

    implicit none

    ! arguments..
    double precision,intent(in)  :: tleaf

    ! local variables..
    double precision :: ci_diff,ad, netrad, vcmt, vjmt, x1, x2, xacc, vpd_min, &
                        ci0, ci1, tol, max_declin, aoa, max_dayl, daylength, dayl_factor, &
                        g0, g1, gsx, gbv, vpd_term, anx, cs, gleaf, convert, &
                        es, ea, rh, de, psych, rho, cp, lambda 


    !----------------------------------------------------------------------
    ! leaf respiration (umolC.m-2.s-1)
    lt = tleaf
    resp = rn * nit * exp( log(2d0) * ( lt - 10d0 ) / 10d0 )

    ! calculate daylength factor
    max_declin = 0.409571d0
    aoa = (sin(latitude_radians) * sin(max_declin)) / (cos(latitude_radians)*cos( max_declin ))
    max_dayl = 12.0d0 * ( 1d0 + 2d0*asin( aoa ) / pi )
    daylength = calculate_daylength_hours(dble(time%day),latitude_radians)
    dayl_factor = min(1d0,max(0.01d0,(daylength*daylength)/(max_dayl*max_dayl)))

    ! ----------------------------------------------------------------------
    !temperature modifications to parameters

    vcmt  = tempmet(vcmax_max_temp,metabolic_opt_temp,Vc_kurtosis,lt)    
    vjmt  = tempmet(jmax_max_temp,metabolic_opt_temp,Vj_kurtosis,lt)
    vcmax = vcmt * vcm
    vjmax = vjmt * vjm 
    !vcmax = vcmt * vcm * dayl_factor
    !vjmax = vjmt * vjm * dayl_factor
   
    ! hydraulic limitaiton on vcmax
    !vcmax = vcmax * ((1+exp(-5.28d0*2.31d0))/(1+exp(-2.31d0-lwp_pd)))


    ! Temperature adjustments for Michaelis-Menten coefficients 
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    kc     = arrhenious(kc_saturation,kc_half_sat_conc,lt)
    ko     = arrhenious(ko_saturation,ko_half_sat_conc,lt)
    gamma1 = arrhenious(co2comp_saturation,co2comp_half_sat_conc,lt)

    ! ---------------- relative humidity from mixing ratio ------------
    rho = air_density_kg * 1d3  ! density of air g/m3 (t-dependent)
    cp = cp_air * 1d-6 ! convert to correct units
    psych = 0.1d0 * ( 0.646 * exp( 0.00097 * lt ) )   ! phsychometric constant is temperature dependent
    lambda = 0.001d0 * ( 2501d0 - 2.364d0 * lt )  ! latent heat of vapor

    de = wdef * lambda * psych / (rho * cp )  ! VPD (kPa)
    es = .61078d0 * exp( 17.269d0 * lt / ( 237.3d0 + lt) ) ! saturation es (kPa)
    ea = max(0d0,es - de)
    ea = min( max(ea,0.2d0*es),es) ! constrain ea to be <es and >0.05es


    gbv = gbb * (atmos_press / (Rcon * (lt+273.15d0))) ! BLC to H2O mol/m2/s
    g0 = 0.01d0
    !g1 = cond_slope*(0.2d0 + 0.8d0/(1+exp(-(lwp_pd+cond_slope_c)/cond_slope_b)))
    g1 = cond_slope*(0.2d0 + 0.8d0*exp(-(-psil/cond_slope_b)**cond_slope_c))
    g1_med = g1
    ! 
    ! -------------------------
    ! Ci calculation with while loop
    ci0 = 0.7*co2amb
    ci1 = 0.99d0*ci0
    
    do while (abs(ci0-ci1)>0.00001)

        ! calculate rh and vpd with gs feedback
        rh = (gbv*ea + gsx*es)/((gbv+gsx)*es)
        vpd_min = .1d0 ! kPa
        vpd_term = max((es-rh*es),vpd_min) ! kPa

        ci0 = ci1
        ! ------------ metabolic based photosynthesis umol/m2/s
        if (c3path ) then
            anx = farquhar( vcmax, vjmax, kc, ko, gamma1, resp, par, ci1)
        else
            anx = collatz( vcmax, resp, par, ci1)
        end if

        if (anx .gt. 0d0 ) then
            ! ----------- Ball-berry/Medlyn ci calculation -----------------
            cs = co2amb - anx/gbv
           
            if ( user_opts%iWUE == 3 ) then ! ball-Berry
                gsx = g1 * rh * ( anx / cs )
            else if ( user_opts%iWUE == 4 ) then
                gsx = 1.6d0 * (1d0 + g1 / sqrt(vpd_term)) * ( anx / cs )! Medlyn
            end if
        else
            gsx = g0
        end if
        
        gleaf = 1d0 / (1.4d0/gbv + 1.6d0/gsx)
            
        ci1 = co2amb - anx/gleaf

        if (anx .lt. 0d0 ) then
            ci1=ci0
        end if
  
    end do

    ci = ci1


    ! update variables outside this module

    ! --------------- convert gsx from molco2/m2leaf/s to m/s 
    convert = (Rcon * (lt+273.15d0)) / atmos_press
    gs2 = gsx*convert ! common block gs in m/s used outside this function

    if ( c3path ) then
        an = farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , ci )   !stable gs 
    else
        an = collatz( vcmax , resp , par , ci )
    end if

    !  determine if radiation is isothermal or net radiation
    if (rad_pass > 1) then
        ! radiation already netrad
        netrad = rad
    else
        ! determine net radiation from isothermal, longwave correction
        ! (kW.m-2) based on temperature difference (Jones, p108)
        netrad = rad - 4d0 * emiss * boltz_kW * ( temp + 273.15d0 )**3 * ( lt - temp )
    end if

    ! evaporation rate in g m-2 s-1, converted to clayer total mol m-2 s-1...
    et = max(0d0,trans( gs2 , lt , netrad , wdef , gbb )*g_to_mol_water)

    ! explicit return
    return

  end subroutine empirical_stomata 

  !
  !----------------------------------------------------------------------
  !
  double precision function TleafFunc ( leaftemp_in )

    use gv_clim, only: atmos_press
    use gv_meteo, only: gbb, wdef,rad_pass, temp,rad
    use math_tools, only: zbrent
    use gv_veg, only: emiss
    use gv_scale_declarations, only: boltz_kW,user_opts,g_to_mol_water
    !
    ! Description: 
    ! Calculate leaf fluxes from input leaf temperature (tleaf_val) and compare
    ! the new temperature to the prior temperature. This function equals zero
    ! when tleaf does not change between iterations. 

    implicit none

    ! arguments..
    double precision, intent(in) :: leaftemp_in
    
    ! local variables..
    double precision :: tleaf2, tleaf1, netrad

    tleaf1 = leaftemp_in

    call empirical_stomata( leaftemp_in )

    ! select method for calculating surface energy balance
    if (user_opts%solve_canopy_temperature == 1) then
       lt = leaf_temperature( gs2 )
    elseif (user_opts%solve_canopy_temperature == 2) then
       lt = zbrent('assimilate:leaf_balance',leaf_balance,temp-20d0,temp+20d0,leaf_temp_precision)
    end if

    tleaf2 = lt

    TleafFunc = tleaf2 - tleaf1

    ! explicit return
    return

  end function TleafFunc

end module leaf
!
!------------------------------------------------------------------------
!
