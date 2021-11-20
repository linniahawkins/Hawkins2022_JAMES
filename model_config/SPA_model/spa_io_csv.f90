! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_io_csv

  !! contains routines specific to reading/writing csv files. !!

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: close_output_csv, open_output_csv, read_crops_csv, read_met_csv, read_phenology_csv, &
            read_carbon_csv, read_soil_csv, read_veg_csv,write_assimilate_output_csv, write_crops_output_csv,  &
            write_ecosystem_fluxes_output_csv,write_daily_ecosystem_fluxes_output_csv,write_daily_stocks_output_csv, & 
            write_daily_fluxes_output_csv,write_daily_water_output_csv,write_solar_output_csv, write_soil_output_csv

  ! to keep track of whether write-routine already called this timestep..
  double precision,allocatable :: laste_call(:),last_call_shade(:), last_call_sun(:)

  integer,parameter :: input_crops_unit       =  20, &
                       input_met_unit         =  21, &  ! Boundary condition
                       input_pheno_unit       =  43, &
                       input_soil_unit        =  22, &  ! Initial
                       input_veg_unit         =  23, &  !    conditions
                       input_carbon_unit      =  44, &  !
                       daily_unit             =  25, &
                       energy_unit            =  26, &
                       fluxes_unit            =  27, &
                       EcoFluxes_unit         =  29, &
                       ice_prop_unit          =  30, &
                       partition_unit         =  31, &
                       solar1_unit            =  32, &
                       solar2_unit            =  33, &
                       solar3_unit            =  34, &
                       solar4_unit            =  35, &
                       soil_status_unit       =  36, &
                       soil_temp_unit         =  37, &
                       soil_water_unit        =  38, &
                       statistics_unit        =  39, &
                       stocks_unit            =  40, &
                       up_frac_unit           =  41, &
                       water_fluxes_unit      =  42, &
                       snow_mass_unit         =  65, &
                       snow_temp_unit         =  66, &
                       layer_sun_start_unit   = 100, &  ! units 100-139 are reserved for 'layer_sun_[i].csv' files
                       layer_shade_start_unit = 140     ! units 140-179 are reserved for 'layer_shade_[i].csv' files
                       ! units 71 & 72 are used in spa_restart.

  save

contains
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! public procedures first..
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine close_output_csv

     use gv_scale_declarations, only: grid
     use gv_veg,                only: lafrac

     implicit none

     ! local variables..
     integer :: i

     close( unit=daily_unit        , status='keep' )
     close( unit=energy_unit       , status='keep' )
     close( unit=fluxes_unit       , status='keep' )
     close( unit=EcoFluxes_unit    , status='keep' )
     close( unit=ice_prop_unit     , status='keep' )
     do i = 1 , grid%canopy
       if (lafrac(i) .ne. 0.) close(unit=layer_sun_start_unit+i  ,status='keep')
       if (lafrac(i) .ne. 0.) close(unit=layer_shade_start_unit+i,status='keep')
     enddo
     close( unit=partition_unit    , status='keep' )
     close( unit=solar1_unit       , status='keep' )
     close( unit=solar2_unit       , status='keep' )
     close( unit=solar3_unit       , status='keep' )
     close( unit=solar4_unit       , status='keep' )
     close( unit=soil_status_unit  , status='keep' )
     close( unit=soil_temp_unit    , status='keep' )
     close( unit=soil_water_unit   , status='keep' )
     close( unit=statistics_unit   , status='keep' )
     close( unit=stocks_unit       , status='keep' )
     close( unit=up_frac_unit      , status='keep' )
     close( unit=water_fluxes_unit , status='keep' )
     close( unit=snow_mass_unit    , status='keep' )
     close( unit=snow_temp_unit    , status='keep' )
 
  end subroutine close_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine open_output_csv( outdir , plant_func_flag )

    ! This opens the standard output files !
    !  and writes a header to each.        !

    use gv_scale_declarations, only: fname_length, grid
    use gv_veg,                only: lafrac
    use gv_carbon_model,       only: fluxes_names, pools_names, nofluxes, nopools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: outdir
    integer,intent(in)          :: plant_func_flag

    ! local variables..
    character(fname_length) :: filename
    character(600)          :: full_header
    integer                 :: i

   ! Open all output files (alphabetical order!)..
    call open_file( trim(outdir)//'daily.csv', daily_unit , header = &
       'mod_ET(mm/d),mod_LE(MJ/m2/d),mod_ET(MJ.m-2.d-1),mod_SoilE(mm/d),mod_can_evap(mm/d),' &
     //'weighted_SWP(MPa),plant_resistance_layer1,canopy_soil_resistance_layer1,'          &
     //'LSC(mmol/m2/s/Mpa),sapflow(mm/d),mean-LAI' , daily = .true. )

    call open_file( trim(outdir)//'energy.csv', energy_unit , header = &
                     'Qh(W/m2),Qe(W/m2),Qn(W/m2),Qc(W/m2),Qm(W/m2),Qs(W/m2),air_temp(oC),surface_temp(oC),soilt2(oC),drythick(mm)' )

    call open_file( trim(outdir)//'ecosystem_fluxes.csv', EcoFluxes_unit , header = &
       'gpp(gC/m2/t),ra(gC/m2/t),rh1(gC/m2/t),rh2(gC/m2/t),tsk(K),canopy_hfx(W/m2),' &
     //'hfx(W/m2),le(W/m2),trans(W/m2),ess(W/m2),wetle(W/m2),gpp(umol/m2/s),nee(umol/m2/s)' )

    call open_file( trim(outdir)//'iceprop.csv', ice_prop_unit , header=&
            'i1(frac),i2(frac),i3(frac),i4(frac),i5(frac),i6(frac),i7(frac),i8(frac),i9(frac),i10(frac),i11(frac),i12(frac),i13(frac),i14(frac),i15(frac)' )

    do i = 1 , grid%canopy
      if ( lafrac(i) .ne. 0. ) then 
        write(filename,"(A,i0.2,A)")"layer_sun_",i,".csv"
        call open_file( trim(outdir)//trim(filename), layer_sun_start_unit+i , header = &
                        'gs(mmol.m-2.s-1),agr(umol.m-2.s-1),res(umol.m-2.s-1),psil(MPa),' &
                      //'ci(ppm),etr(W.m-2),tempdf(oC),wdef(g.m-3),par(umol.m-2.s-1),rad(kW.m-2),' &
                      //'coa(ppm),la(m2/m2),an(umol.m-2.s-1),cistar(mmol.m-2.s-1),N(gN.m-2),LCA(gC.m-2),D13C,gsx(mmol.m-2.s-1),vpd_term(kPa),g1,lwp_pd(MPa)' )
        write(filename,"(A,i0.2,A)")"layer_shade_",i,".csv"
        call open_file( trim(outdir)//trim(filename), layer_shade_start_unit+i , header = &
                        'gs(mmol.m-2.s-1),agr(umol.m-2.s-1),res(umol.m-2.s-1),psil(MPa),' &
                      //'ci(ppm),etr(W.m-2),tempdf(oC),wdef(g.m-3),par(umol.m-2.s-1),rad(kW.m-2),' &
                      //'coa(ppm),la(m2/m2),an(umol.m-2.s-1),cistar(mmol.m-2.s-1),N(gN.m-2),LCA(gC.m-2),D13C,gsx(mmol.m-2.s-1),vpd_term(kPa),g1,lwp_pd(MPa)' )
      end if
    enddo

    if ( plant_func_flag == 2 .or. plant_func_flag == 3) then ! just crops..
      call open_file( trim(outdir)//'partitioning.csv', partition_unit , header=&
                      'shoot,root,leaves,stem,DS,DR,lai,fT,fP,fV,VD,raso,max_raso' )
      call open_file( trim(outdir)//'statistics.csv' , statistics_unit , header=&
                      'NEE_harvest,NEE_endofyear,yield,exp_residue,'//&
                      'total_exp,NBP_harvest,NBP_endofyear' )
    end if

    call open_file( trim(outdir)//'solar_part1.csv', solar1_unit , header=&
                   'soilnet,soilabsn,soilabls,soil_absorbed_ppfd,par_top' )

    call open_file( trim(outdir)//'solar_part2.csv', solar2_unit , header=&
       'sum(abspar_sun),(1.-fdiff)*par_top,(parabs_sun(i),i=1,10)' )

    call open_file( trim(outdir)//'solar_part3.csv', solar3_unit , header=&
     'sum(abspar_shade),fdiff*par_top,(parabs_shade(i),i=1,10),skyabsp' )

    call open_file( trim(outdir)//'solar_part4.csv', solar4_unit , header=&
                           'fdiff,(leaf_sun(i),i=1,10),check' )

    call open_file( trim(outdir)//'soil_status.csv', soil_status_unit , header=&
                   'runoff(mm),totet(mm),drainage(mm),surfwater(mm),total_water(mm),' &
                 //'within_soil_flux(mm),outwith_soil_flux(mm),flux_check(mm),watergain(1;mm),' &
                 //'waterloss(1;mm),del_soilwater_stock(mm),stock_flux_check(mm),total_precip_gain(mm),' &
                 //'total_water_gains(mm),total_water_losses(mm),water->core(mm),core->water(mm),snow(kg/m2),snow_height(m)' )

    call open_file( trim(outdir)//'soil_temp.csv', soil_temp_unit , header=&
       'soil-temp(layer 1),soil-temp(layer 2),soil-temp(layer 3),'     &
     //'soil-temp(layer 4),soil-temp(layer 5),soil-temp(layer 6),'     &
     //'soil-temp(layer 7),soil-temp(layer 8),soil-temp(layer 9),'     &
     //'soil-temp(layer 10),soil-temp(layer 11),soil-temp(layer 12),'  &
     //'soil-temp(layer 13),soil-temp(layer 14),soil-temp(layer 15),soil-temp(layer 21)' )

    call open_file( trim(outdir)//'snow_mass.csv', snow_mass_unit , header=&
       'ice-mass(layer_1),ice-mass(layer_2),ice-mass(layer_3),ice-mass(layer_4),ice-mass(layer_5),'  &
     //'liquid-mass(layer_1),liquid-mass(layer_2),liquid-mass(layer_3),liquid-mass(layer_4),liquid-mass(layer_5)')
    call open_file( trim(outdir)//'snow_temp.csv', snow_temp_unit , header=&
       'snow-temp(layer_1),snow-temp(layer_2),snow-temp(layer_3),snow-temp(layer_4),snow-temp(layer_5),'  &
     //'layer-height(layer_1),layer-height(layer_2),layer-height(layer_3),layer-height(layer_4),layer-height(layer_5)')

    call open_file( trim(outdir)//'soilwater.csv',  soil_water_unit , header=&
     'w1(m3/m3),w2(m3/m3),w3(m3/m3),w4(m3/m3),w5(m3/m3),w6(m3/m3),w7(m3/m3),w8(m3/m3),w9(m3/m3),w10(m3/m3),w11(m3/m3),w12(m3/m3),w13(m3/m3),w14(m3/m3),w15(m3/m3),w_swp(MPa),swp1(MPa),swp2(MPa),swp3(MPa),swp4(MPa),swp5(MPa),swp6(MPa),swp7(MPa),swp8(MPa),swp9(MPa),swp10(MPa)' )

    call open_file( trim(outdir)//'upfrac.csv',     up_frac_unit , header=&
     'up1,up2,up3,up4,up5,up6,up7,up8,up9,up10,up11,up12,up13,up14,up15' )

    call open_file( trim(outdir)//'waterfluxes.csv', water_fluxes_unit , &
                    header='modess,runoff,ppt,trans,delw,disc,cwtr,check,diff,modwet,unint,canst,snow_watermm' , &
                    daily=.true. )

    full_header=trim(fluxes_names(1))
    do i = 2, nofluxes
       full_header=trim(full_header)//','//trim(fluxes_names(i))
    enddo
    call open_file( trim(outdir)//'fluxes.csv', fluxes_unit , header=full_header, daily=.true. )

    full_header=trim(pools_names(1))
    do i = 2, nopools
       full_header=trim(full_header)//','//trim(pools_names(i))
    enddo
    call open_file( trim(outdir)//'stocks.csv', stocks_unit , header=full_header, daily=.true. )

  end subroutine open_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_carbon_csv( filename )

    ! open and read the contents of the carbon.csv file. !
    use gv_carbon_model,       only: pars, nopars, pars_names
    use gv_veg,                only: avN, stock_roots
    use gv_scale_declarations, only: fname_length 
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename

    ! local variables..
    character(len=200) :: header
    integer            :: i

    ! Open the plant physiological parameters
    call open_file( filename , input_carbon_unit , readonly=.true. )

    ! read in carbon model parameters
    do i = 1, nopars
       read(unit=input_carbon_unit,fmt=*)header,pars(i)
       if (trim(header) /= pars_names(i)) then
           print*,trim(header)," found in carbon parameter file when expecting ",trim(pars_names(i)) ; stop
       endif
    end do

    ! overwrite the avgN value from the *_veg.csv when using DALEC
    avN = 10d0**(pars(11))
    ! pass roots initial condition 
    stock_roots = pars(20)

    ! appropriate sanity checks have yet to be developmed for this....

  end subroutine read_carbon_csv
  !
  !-----------------------------------------------------------------------------------------------------------------------------
  !
  subroutine read_crops_csv( filename )

    ! Read the contents of the development.csv file !

    use carbon_model_crop_mod, only: DR_P, DR_T_PRA, DR_T_POA, DRAO_P, DRAO_T_POA, DRAO_T_PRA, &
                                     DS_LRLV, DS_LRRT, DS_root, DS_shoot, fol_frac, LCA_DS,    &
                                     LCA_ratio, LRLV, LRRT, root_frac, stem_frac, &
                                     stock_seed_labile
    use gv_scale_declarations, only: fname_length

    implicit none

    ! arguments..
    character(len=*) :: filename

    ! local variables..
    integer                 :: columns, i, rows
    character(fname_length) :: variables

    ! crop development file
    call open_file( filename , input_crops_unit , readonly=.true. )

    ! read in the amount of carbon available (as labile) in each seed..
    read(unit=input_crops_unit,fmt=*) variables,stock_seed_labile,variables,variables

    ! read in C partitioning/fraction data and corresponding developmental stages (DS)
    ! shoot
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_shoot(rows) , fol_frac(rows) , stem_frac(rows)  )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_shoot(i), fol_frac(i), stem_frac(i)
    enddo

    ! root
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_root(rows) , root_frac(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_root(i), root_frac(i)
    enddo

    ! loss rates of leaves and roots
    ! leaves
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRLV(rows) , LRLV(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRLV(i), LRLV(i)
    enddo

    ! roots
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DS_LRRT(rows) , LRRT(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DS_LRRT(i), LRRT(i)
    enddo

    ! developmental rate as a function of temperature (not needed when calculated through modified Wang&Engel model,
    ! which is the current set-up)
    !preanthesis (before flowering)
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DR_T_PRA(rows) , DRAO_T_PRA(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DR_T_PRA(i), DRAO_T_PRA(i)
    enddo

    ! postanthesis
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DR_T_POA(rows) , DRAO_T_POA(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DR_T_POA(i), DRAO_T_POA(i)
    enddo
 
    ! photoperiod (daylength effect on development)
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( DR_P(rows) , DRAO_P(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) DR_P(i), DRAO_P(i)
    enddo

    ! LCA ratios (if LCA to be dynamic, but currently not used)
    read(unit=input_crops_unit,fmt=*) variables
    read(unit=input_crops_unit,fmt=*) rows , columns
    allocate( LCA_DS(rows)  , LCA_ratio(rows) )
    do i = 1 , rows
      read(unit=input_crops_unit,fmt=*) LCA_DS(i), LCA_ratio(i)
    enddo

  end subroutine read_crops_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_met_csv( filename , met )

    ! open the met driver (csv) file, and load the data from it !

    use gv_scale_declarations, only: met_drivers, time, user_opts
    use gv_carbon_model, only: twq
    use gv_hourscale, only: freeze
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    type(met_drivers)           :: met

    ! local variables..
    character(len=250) :: line
    character(len=250), dimension(25) ::header
    integer            :: i, j, num_lines, status, &
                          airt_col,lwrad_col,ppfd_col, &
                          precip_col,sfc_pressure_col, &
                          snowfall_col,swrad_col, &
                          vpd_col,sp_moist_col, &
                          wind_spd_col,co2_col, &
                          clearfrac_col,clearman_col, &
                          firefrac_col, &
                          sw_diffuse_col, &
                          swc10_col, swc20_col, swc30_col, swc40_col, &
                          swc50_col, swc60_col, swc70_col, swc80_col, &
                          swc90_col, swc100_col, swc110_col, swc120_col, &
                          swc130_col, swc140_col, swc150_col, &                          
                          num_drivers
    double precision  :: temp_var(25)

    ! open file
    call open_file( filename , input_met_unit , readonly=.true. )

    ! get header out of the way..
    read(unit=input_met_unit,fmt=*) line

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 0
    do
      read(input_met_unit,fmt=*,iostat=status) line
      if ( status .ne. 0. ) exit
      num_lines = num_lines + 1
    enddo

    ! make the variables the correct size..
    allocate( met%ambient_co2( num_lines ) , &
                   met%lw_rad( num_lines ) , &
                     met%ppfd( num_lines ) , &
                   met%precip( num_lines ) , &
                   met%sw_rad( num_lines ) , &
             met%sfc_pressure( num_lines ) , &
                     met%temp( num_lines ) , &
                      met%vpd( num_lines ) , &
                 met%wind_spd( num_lines ) , &
       met%clearance_fraction( num_lines ) , &
     met%clearance_management( num_lines ) , &
                    met%swc10( num_lines ) , &
                    met%swc20( num_lines ) , &
                    met%swc30( num_lines ) , &
                    met%swc40( num_lines ) , &
                    met%swc50( num_lines ) , &
                    met%swc60( num_lines ) , &
                    met%swc70( num_lines ) , &
                    met%swc80( num_lines ) , &
                    met%swc90( num_lines ) , &
                   met%swc100( num_lines ) , &
                   met%swc110( num_lines ) , &
                   met%swc120( num_lines ) , &
                   met%swc130( num_lines ) , &
                   met%swc140( num_lines ) , &
                   met%swc150( num_lines ) , &
            met%fire_fraction( num_lines ) )

    if (user_opts%met_file_has_sw_diffuse) then
        allocate(met%sw_diffuse(num_lines))
    end if
    if (user_opts%met_file_has_snowfall) then
        allocate(met%snowfall(num_lines))
    endif


    ! Go back to start of file..
    rewind( input_met_unit )

    ! get header out of the way again but at the same time record the header
    ! names (MATCH WITH THE DO LOOP BELOW)
    read(unit=input_met_unit,fmt=*) header(1),header(2),header(3),header(4) &
                                   ,header(5),header(6),header(7),header(8) &
                                   ,header(9),header(10),header(11),header(12) &
                                   ,header(13),header(14),header(15) &
                                   ,header(16), header(17), header(18) &
                                   ,header(19), header(20), header(21) &
                                   ,header(22), header(23), header(24) &
                                   ,header(25)

    ! set initial values
    airt_col = 0 ; lwrad_col = 0 ; ppfd_col = 0 ; precip_col = 0
    sfc_pressure_col = 0 ; snowfall_col = 0 ; swrad_col = 0
    vpd_col = 0 ; sp_moist_col = 0 ; wind_spd_col = 0
    co2_col = 0 ; sw_diffuse_col = 0 ; firefrac_col = 0
    clearfrac_col = 0 ; clearman_col = 0 ; num_drivers = 0
    swc10_col = 0 ; swc20_col = 0 ; swc30_col = 0 
    swc40_col = 0 ; swc50_col = 0 ; swc60_col = 0
    swc70_col = 0 ; swc80_col = 0 ; swc90_col = 0
    swc100_col = 0 ; swc110_col = 0 ; swc120_col = 0
    swc130_col = 0 ; swc140_col = 0 ; swc150_col = 0

    ! use the headers in the file to determine which columns contain which
    ! drivers. NOTE also that this system assumes that the last identified
    ! drivers must be the widest column we want to read
    do i = 1, 25
       if (trim(header(i)) == "airt") then
          airt_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "lw_rad") then
          lwrad_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "ppfd" .or. trim(header(i)) == "PAR" .or. trim(header(i)) == "par") then
          ppfd_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "precip") then
          precip_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "sfc_pressure") then
          sfc_pressure_col=i ;  num_drivers = i 
       endif
       if (trim(header(i)) == "snowfall") then
           snowfall_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "sw_rad") then
           swrad_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "sw_diffuse") then
           sw_diffuse_col = i ; num_drivers = i
       endif
       if (trim(header(i)) == "vpd") then
          vpd_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "sp_moist") then
          sp_moist_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "wind_spd") then
          wind_spd_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "co2") then
          co2_col=i ; num_drivers = i 
       endif
       if (trim(header(i)) == "clearance_fraction") then
          clearfrac_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "clearance_management") then
          clearman_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "fire_fraction") then
          firefrac_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc10") then
          swc10_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc20") then
          swc20_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc30") then
          swc30_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc40") then
          swc40_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc50") then
          swc50_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc60") then
          swc60_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc70") then
          swc70_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc80") then
          swc80_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc90") then
          swc90_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc100") then
          swc100_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc110") then
          swc110_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc120") then
          swc120_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc130") then
          swc130_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc140") then
          swc140_col=i ; num_drivers = i
       endif
       if (trim(header(i)) == "swc150") then
          swc150_col=i ; num_drivers = i
       endif 
  
    enddo ! checking driver columns

    ! rewind the file and read over the first line
    ! Go back to start of file..
    rewind( input_met_unit )
    ! get header out of the way again
    read(unit=input_met_unit,fmt=*) line

    ! check to ensure that we have located all the variables, or at least all
    ! the mandantory drivers
    if (airt_col == 0) then
        print*,"have not located 'airt' driver"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    end if
    if (swrad_col == 0) then
        print*,"have not located 'sw_rad' driver"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (precip_col == 0) then
        print*,"have not located 'precip' driver"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (snowfall_col == 0 .and. user_opts%met_file_has_snowfall) then
        print*,"'snowfall' could not be found but config indicates that we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (wind_spd_col == 0) then
        print*,"have not located 'wind_spd' driver"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (ppfd_col == 0 .and. user_opts%use_ppfd_from_met_file) then
        print*,"'PAR' could not be found in met file but config indicates we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
    endif
    if (sw_diffuse_col == 0 .and. user_opts%met_file_has_sw_diffuse) then
        print*,"'sw_diffuse' could not be found in met file but config indicates we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
    endif
    if (lwrad_col == 0 .and. user_opts%met_file_has_lw_rad) then
        print*,"'lw_rad' could not be found in met file but config indicates we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (vpd_col == 0 .and. user_opts%vpd_or_specific_humidity == "vpd") then
        print*,"'vpd' was not found met file but config indicates we are expecting VPD"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (sp_moist_col == 0 .and. user_opts%vpd_or_specific_humidity == "specific_humidity") then
        print*,"'sp_moist' was not found met file but config indicates we are expecting specific_humidity"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (sfc_pressure_col == 0 .and. user_opts%vpd_or_specific_humidity == "specific_humidity") then
        print*,"you had indicated that you want to use specific_humidity as a driver, &
              &however you have not provided sfc_pressure. To use specific_humidity &
              &I am afraid you need the correct sfc_pressure or things will go very wrong!"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (clearfrac_col == 0 .and. user_opts%met_file_has_disturbance) then
        print*,"'clearance_fraction' could not be found in met file but config indicates we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (clearman_col == 0 .and. user_opts%met_file_has_disturbance) then
        print*,"'clearance_management' could not be found in met file but config indicates we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif
    if (firefrac_col == 0 .and. user_opts%met_file_has_disturbance) then
        print*,"'fire_fraction' could not be found in met file but config indicates we are expecting it"
        print*,"accepted variable names are: airt, sw_rad, sw_diffuse, precip, snowfall, &
               &PAR, sfc_pressure, wind_spd, lw_rad, co2, clearance_fraction, &
               &clearance_management, fire_fraction and vpd or sp_moist"
        stop
    endif

    !  now read back through the met file to read in the correct drivers
    do i = 1 , num_lines

       ! first read in the data from each column blindly to a dummy variable
       read(unit=input_met_unit,fmt=*) (temp_var(j), j = 1, num_drivers)

       ! then allocate these based on the above column markers to the correct
       ! output
       if (vpd_col /= 0 .and. user_opts%vpd_or_specific_humidity == "vpd") then
           met%vpd(i)=temp_var(vpd_col)
       endif
       if (sp_moist_col /= 0 .and. user_opts%vpd_or_specific_humidity == "specific_humidity") then
           met%vpd(i)=temp_var(sp_moist_col)
       endif
       if (lwrad_col /= 0 .and. user_opts%met_file_has_lw_rad) then
           met%lw_rad(i)=temp_var(lwrad_col)
       endif
       if (sw_diffuse_col /= 0 .and. allocated(met%sw_diffuse)) then
           met%sw_diffuse(i)=temp_var(sw_diffuse_col)
       endif
       if (snowfall_col /= 0 .and. allocated(met%snowfall)) then
           met%snowfall(i)=temp_var(snowfall_col)
       endif
       if (clearfrac_col /= 0 .and. allocated(met%clearance_fraction)) then
           met%clearance_fraction(i)=temp_var(clearfrac_col)
           met%clearance_management(i)=int(temp_var(clearman_col))
           met%fire_fraction(i)=temp_var(firefrac_col)
       endif

       if (airt_col /= 0) met%temp(i)=temp_var(airt_col)
       if (co2_col /= 0) met%ambient_co2(i)=temp_var(co2_col)
       if (wind_spd_col /= 0) met%wind_spd(i)=temp_var(wind_spd_col)
       if (swrad_col /= 0) met%sw_rad(i)=temp_var(swrad_col)
       if (ppfd_col /= 0) met%ppfd(i)=temp_var(ppfd_col)
       if (precip_col /= 0) met%precip(i)=temp_var(precip_col)
       if (sfc_pressure_col /= 0) met%sfc_pressure(i)=temp_var(sfc_pressure_col)
       if (swc10_col /= 0) met%swc10(i) = temp_var(swc10_col)
       if (swc20_col /= 0) met%swc20(i) = temp_var(swc20_col)
       if (swc30_col /= 0) met%swc30(i) = temp_var(swc30_col)
       if (swc40_col /= 0) met%swc40(i) = temp_var(swc40_col)
       if (swc50_col /= 0) met%swc50(i) = temp_var(swc50_col)
       if (swc60_col /= 0) met%swc60(i) = temp_var(swc60_col)
       if (swc70_col /= 0) met%swc70(i) = temp_var(swc70_col)
       if (swc80_col /= 0) met%swc80(i) = temp_var(swc80_col)
       if (swc90_col /= 0) met%swc90(i) = temp_var(swc90_col)
       if (swc100_col /= 0) met%swc100(i) = temp_var(swc100_col)
       if (swc110_col /= 0) met%swc110(i) = temp_var(swc110_col)
       if (swc120_col /= 0) met%swc120(i) = temp_var(swc120_col)
       if (swc130_col /= 0) met%swc130(i) = temp_var(swc130_col)
       if (swc140_col /= 0) met%swc140(i) = temp_var(swc140_col)
       if (swc150_col /= 0) met%swc150(i) = temp_var(swc150_col)

    end do ! for num_lines

    ! Check sanity of temperature..
    if (user_opts%temp_in_kelvin .and. minval(met%temp) .lt. 150d0) then
       write(message,*) "Config file says we are expecting kelvin but met drive appear not..."
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! Check sanity of temperature..
    if (.not.user_opts%temp_in_kelvin .and. maxval(met%temp) .gt. 100d0) then
       write(message,*) "Config file says we are expecting oC but met drive appear not..."
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! change air temperature to oC if needed
    if (user_opts%temp_in_kelvin) met%temp = met%temp - freeze

    ! Is precip in volume-per-timestep (mm), or rate (mm/sec)?
    if ( user_opts%precip_is_rate ) then
       call write_log("Converting precipitation from rate to volume by "&
                   //"multiplying by nos of seconds per timestep.")
       met%precip = met%precip * time%seconds_per_step
    end if

    ! Is precip in volume-per-timestep (mm), or rate (mm/sec)?
    ! assume same value for snowfall
    if ( user_opts%met_file_has_snowfall .and. .not.user_opts%precip_is_rate ) then
       call write_log("Converting snowfall to rate from volume by "&
                   //"dividing by nos of seconds per timestep.")
       met%snowfall = met%snowfall / time%seconds_per_step
    end if

    ! If Surface pressure not provided..
    if (.not. user_opts%met_file_has_sfc_press ) then
        met%sfc_pressure = met%const_sfc_pressure
        if (sfc_pressure_col /= 0) then
            print*,'WARNING: config file indicates that no surface pressure is provided,' &
                   //' however the met file contains "sfc_pressure" column'
        endif
    endif

    ! In case user doesn't want to use co2 in the file..
    if ( .not. user_opts%use_co2_from_met_file ) then
        met%ambient_co2 = met%const_ambient_co2
        if (co2_col /= 0) then
            print*,"WARNING: config file indicates that no co2 is provided, however the met file contains 'co2' column"
        endif
    endif

    ! In case user doesn't want to use the par in the file..
    if ( .not. user_opts%use_ppfd_from_met_file )  then
        met%ppfd = 2.3d0 * met%sw_rad ! W.m-2 SW -> PAR (umol.m-2s-1)
        if (ppfd_col /= 0) then
            print*,'WARNING: config file indicates that no PAR / ppfd is provided, ' &
                   //' however the met file contains "ppfd" or "PAR" column'
        endif
    endif

    ! Adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2d0 ) met%wind_spd = 0.2d0
    call write_log('Please note that any wind speeds below 0.2 m.s-1 will be' &
                 //'assumed to be = 0.2')

  end subroutine read_met_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_phenology_csv( filename , pheno )

    ! open the met driver (csv) file, and load the data from it !

    use gv_scale_declarations, only: phenology_drivers
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    type(phenology_drivers)     :: pheno

    ! local variables..
    character(len=200) :: line
    integer            :: i, num_lines, status
    double precision   :: dummy

    ! open file
    call open_file( filename , input_pheno_unit , readonly=.true. )

    ! get header out of the way..
    read(unit=input_pheno_unit,fmt=*) line

    ! count the number of remaining lines in the file..
    status = 0 ; num_lines = 0
    do
      read(input_pheno_unit,fmt=*,iostat=status) line
      if ( status .ne. 0. ) exit
      num_lines = num_lines + 1
    enddo

    ! make the variables the correct size..
    allocate( pheno%lai( num_lines ) , &
              pheno%rootC( num_lines ) )

    ! Go back to start of file..
    rewind( input_pheno_unit )

    ! get header out of the way again..
    read(unit=input_pheno_unit,fmt=*) line

    ! and now load all data..
    do i = 1 , num_lines
      read(unit=input_pheno_unit,fmt=*)  dummy, pheno%lai(i), pheno%rootC(i)
    enddo
!pheno%rootC = pheno%rootC*5
    ! Check sanity of lai
    if ( ( maxval(pheno%lai) > 20d0 ) .or. ( minval(pheno%lai) < 0d0 ) ) then
        write(message,*) "Input lai greater than 20 m2/m2 or negative"
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! Check sanity of root carbon
    if ( minval(pheno%rootC) < 0d0 ) then
        write(message,*) "Input root carbon negative"
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if

  end subroutine read_phenology_csv
  !
  !------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_soil_csv( filename , plant_func_flag )

    ! open and read the contents of the soil.csv file. !
 
    use gv_hydrol,            only: soil_frac_clay, soil_frac_sand
    use gv_scale_declarations,only: grid
    use gv_snow_info,         only: Dsfix, Dsnow, Nsnow, Sice, Sliq, snowheight, &
                                    snowalb_nir, snowalb_par, snowweight, Tsnow
    use gv_soil_structure,    only: iceprop, layer_depth, mineralfrac, organicfrac, &
                                    thickness, waterfrac, soil_temp

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    integer,intent(in)          :: plant_func_flag

    ! local variables..
    character(len=200) :: header
    integer            :: i
    double precision               :: dl

    ! Open the soil parameters file..
    call open_file( filename , input_soil_unit , readonly=.true. )

    ! read soil properties..
    read(unit=input_soil_unit,fmt=*)header
    read(unit=input_soil_unit,fmt=*)header,(thickness(i),  i=1,grid%core)
    if (trim(header) /= "layer_thickness") then
        print*,"layer_thickness: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(layer_depth(i),i=1,grid%core)
    if (trim(header) /= "layer_depth") then
        print*,"layer_depth: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(organicfrac(i),i=1,grid%core)
    if (trim(header) /= "organic_fraction") then 
        print*,"organic_fraction: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(mineralfrac(i),i=1,grid%core)
    if (trim(header) /= "mineral_fraction") then
        print*,"mineral_fraction: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(waterfrac(i),  i=1,grid%core)
    if (trim(header) /= "initial_water_fraction") then
        print*,"initial_water_fraction: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(soil_temp(i),  i=1,grid%core)
    if (trim(header) /= "initial_soil_temp") then 
        print*,"initial_soil_temp: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(iceprop(i),    i=1,grid%core)
    if (trim(header) /= "initial_ice_proportion") then
        print*,"initial_ice_proportion: parameter miss read" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(soil_frac_sand(i),i=1,grid%core)
    if (trim(header) /= "%_of_sand_in_soil_layer") then
        print*,"%_of_sand_in_soil_layer: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,(soil_frac_clay(i),i=1,grid%core)
    if (trim(header) /= "%_of_clay_in_soil_layer") then 
        print*,"%_of_clay_in_soil_layer: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,snowweight
    if (trim(header) /= "snow_water_equivalent") then
        print*,"snow_water_equivalent: parameter miss read in soil parameter file" ; stop
    end if
    read(unit=input_soil_unit,fmt=*)header,snowheight
    if (trim(header) /= "depth_of_snow") then 
        print*,"depth_of_snow: parameter miss read in soil parameter file" ; stop
    end if

    ! initialize multi-layer snow model. These are parameters but are currently
    ! hardcoded for simplicity
    Dsfix(:) = thickness(1)
    Dsnow(:) = 0d0
    Sice(:) = 0d0
    Sliq(:) = 0d0
    Tsnow(:) = 270d0 ! Temperature (K)
    snowalb_nir = 0.73d0 ! Near infer-red reflectance
    snowalb_par = 0.95d0 ! PAR reflectance
    Nsnow = 0
    if (snowheight > 0d0) then
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

  end subroutine read_soil_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine read_veg_csv( filename )

    ! open and read the contents of the veg.csv file. !
    use gv_scale_declarations, only: fname_length, grid, pi
    use gv_soil_structure,     only: max_depth,max_storage,root_k,through_fall, swp_params
    use gv_veg
    use gv_carbon_model,       only: pars
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename

    ! local variables..
    character(len=212) :: header
    integer            :: i, dummy(grid%canopy)

    ! Open the plant physiological parameters
    call open_file( filename , input_veg_unit , readonly=.true. )

    ! read vegetation parameters
    read(unit=input_veg_unit,fmt=*)header,( lafrac(i)       , i=1,grid%canopy ) ! LA fraction
    if (trim(header) /= "lafrac") then
        print*,"lafrac: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header, ( lafrac_dead(i), i=1,grid%canopy ) ! dead lai fraction
    if (trim(header) /= "lafrac_dead") then
        print*,"lafrac_dead: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,( nfrac(i)        , i=1,grid%canopy ) ! nitrogen fraction in each layer
    if (trim(header) /= "nfrac") then 
        print*,"nfrac: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,( layer_height(i) , i=1,grid%canopy ) ! layer heights
    if (trim(header) /= "layer_height") then
        print*,"layer_height: parameter miss read in veg parameter file" ; stop
    end if
    dummy = 1 ! default is c3..
    read(unit=input_veg_unit,fmt=*)header,( dummy(i) , i=1,grid%canopy )      ! PS pathway; 1=>c3, anything else=>c4 
    if (trim(header) /= "c3") then 
        print*,"c3: parameter miss read in veg parameter file" ; stop
    end if 
    do i = 1 , grid%canopy
      if (dummy(i) .eq. 0 ) then
        c3(i) = .false.
      else
        c3(i) = .true.
      end if
    enddo
    read(unit=input_veg_unit,fmt=*)header,avN                         ! average foliar N, gN m-2 leaf area
    if (trim(header) /= "avN") then
        print*,"avN: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,pars(17)                    ! average foliar N, gN m-2 leaf area
    if (trim(header) /= "LCA") then
        print*,"LCA: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,gplant                      ! conductivity
    if (trim(header) /= "gplant") then
        print*,"gplant: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,minlwp                      ! critical LWP
    if (trim(header) /= "minlwp") then
        print*,"minlwp: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,iWUE                        ! stomatal efficiency
    if (trim(header) /= "iWUE") then 
        print*,"iWUE: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,WUE                        ! stomatal efficiency
    if (trim(header) /= "WUE") then
        print*,"WUE: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,cond_slope                 ! Medlyn/Ball-Berry slope param
    if (trim(header) /= "cond_slope") then
        print*,"cond_slope: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,cond_slope_b
    if (trim(header) /= "cond_slope_b") then
        print*,"cond_slope_b: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,cond_slope_c
    if (trim(header) /= "cond_slope_c") then
        print*,"cond_slope_c: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,capac                       ! leaf capacitance    
    if (trim(header) /= "capac") then
        print*,"capac: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,dimen(1)                    ! leaf length
    if (trim(header) /= "leaf_length") then
        print*,"leaf_length: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,dimen(2)                    ! leaf / cone width
    if (trim(header) /= "cone_width") then
        print*,"cone_width: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,rootresist                  ! root resistivity
    if (trim(header) /= "rootresist") then
        print*,"rootresist: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,tower_height                ! height of measurement tower
    if (trim(header) /= "tower_height") then
        print*,"tower_height: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,conductivity                ! Does conductance vary with stem length? 0=NO, 1=YES
    if (trim(header) /= "conductivity") then
        print*,"conductivity: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,kappac                      ! coefficient relating max Vcmax with foliar N
    if (trim(header) /= "kappac") then
        print*,"kappac: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,kappaj                      ! coefficient relating max Jmax with foliar N
    if (trim(header) /= "kappaj") then
        print*,"kappaj: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,max_depth                   ! max rooting depth
    if (trim(header) /= "max_depth") then
        print*,"max_depth: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,root_k                      ! root mass for reaching 50% max depth
    if (trim(header) /= "root_k") then
        print*,"root_k: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,( swp_params(i),i=1,4) ! parameters for SWP
    if (trim(header) /= "swp_params") then
        print*,"swp_params: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,through_fall
    if (trim(header) /= "through_fall") then 
        print*,"through_fall: parameter miss read in veg parameter file" ; stop
    end if
    read(unit=input_veg_unit,fmt=*)header,max_storage
    if (trim(header) /= "max_storage") then
        print*,"max_storage: parameter miss read in veg parameter file" ; stop
    end if
   
  end subroutine read_veg_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_assimilate_output_csv( time, output , in_sun_flag )

    ! moved writes from assimilate (leaf.f90) to here !

    use gv_clim,               only: coa
    use gv_metab,              only: an, ci
    use gv_meteo,              only: gbb, la, par, psil, rad, temp, wdef
    use gv_scale_declarations, only: grid, time_holder
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    double precision,             intent(in) :: output(:)
    logical,          intent(in) :: in_sun_flag

    ! local variables..
    integer :: fileunit, output_layer
    double precision    :: timeflag

    ! If first time through then we need to allocate..
    if ( .not. allocated(last_call_shade) ) then
        allocate(last_call_shade(grid%canopy))
        last_call_shade = 0d0
    end if
    if ( .not. allocated(last_call_sun) ) then
        allocate(last_call_sun(grid%canopy))
        last_call_sun   = 0d0
    end if

    output_layer = int( output(1) )        ! convert back from real to integer
    output_layer = max( output_layer , 1 ) ! make sure index is at least 1

    if ( in_sun_flag ) then
        ! canopy layer is in sunlight..
        timeflag = last_call_sun( output_layer )
        fileunit = layer_sun_start_unit + output_layer
    else
        ! canopy layer is in shade..
        timeflag = last_call_shade( output_layer )
        fileunit = layer_shade_start_unit + output_layer
    end if

    ! Check we have not been called already this time step!
    if ( timeflag .lt. time%daytime ) then

        call write_to_file( fileunit , (/ output(2:4), psil, ci, output(5), output(6)-temp, &
                                          wdef, par, rad, coa, la, an, gbb*output(7),output(8),output(9), output(10), output(11), output(12), output(13), output(14) /) )

        if ( in_sun_flag ) then        
            last_call_sun(output_layer) = time%daytime
        else
            last_call_shade(output_layer) = time%daytime
        end if

    else

!        write(message,*)'write_assimilate_output_csv has already been called this timestep: ignoring additional call'
!        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )

    end if

  end subroutine write_assimilate_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_ecosystem_fluxes_output_csv( time )

    ! write to output files. HARD-CODED. !

    use gv_carbon_model,        only: GPP, resp_auto, resp_h_litter, resp_h_soilOrgMatter
    use gv_hourscale,           only: Qh
    use gv_scale_declarations,  only: time_holder
    use gv_veg,                 only: tsk, can_sensible_layer, lemod, transt, wetle &
                                     ,gppt, neemod, ess, sensible_heat

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    call write_to_file( EcoFluxes_unit , (/ GPP, resp_auto, resp_h_litter, resp_h_soilOrgMatter, tsk, &
                        sum(can_sensible_layer), sensible_heat(time%step), lemod(time%step), transt(time%step), &
                        ess(time%step), wetle(time%step), gppt(time%step), neemod(time%step) /) )

  end subroutine write_ecosystem_fluxes_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_daily_stocks_output_csv( time )

    ! write to output files. HARD-CODED. !

    use gv_carbon_model, only: POOLS
    use gv_scale_declarations,  only: time_holder
    use log_tools

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! output stocks at the end of the day only...
    if (time%step .eq. time%steps_per_day) then
       ! stocks..
       call write_to_file(stocks_unit,POOLS(1,:),daily = .true. )
    else ! this has been called at the wrong time
       write(message,*)"Subroutine 'write_daily_stocks_output_csv' has been called not at the end of the day"
       call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )   
    end if ! if at end of the day

  end subroutine write_daily_stocks_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_daily_fluxes_output_csv( time , plant_func_flag )

    ! write to output files. HARD-CODED. !

    use carbon_model_crop_mod, only: shoot_frac_intpol, root_frac_intpol &
                                      ,fol_frac_intpol, stem_frac_intpol   &
                                      ,DS, DR, fT, fP, fV, VD, raso, max_raso
    use gv_carbon_model
    use gv_scale_declarations,   only: time_holder
    use gv_veg,                  only: lai

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    integer,intent(in)           :: plant_func_flag

    ! local variables..
    integer :: i

    ! output cumulative flux information at the end of the day
    if (time%step .eq. time%steps_per_day) then

        ! if this is the first time of outputting then we need to declare this
        if (.not.allocated(daily_fluxes)) allocate(daily_fluxes(nofluxes))

        ! convert fluxes (gC.m-2.day-1 -> gC.m-2.t-1)
        do i = 1, nofluxes
           daily_fluxes(i) = sum(FLUXES(:,i))
        enddo
        daily_fluxes = daily_fluxes / time%steps_per_day

        ! fluxes..
        call write_to_file( fluxes_unit,daily_fluxes,daily=.true.)

      if ( plant_func_flag == 2 .or. plant_func_flag == 3 ) then ! crops..

        call write_to_file( partition_unit ,                                         &
                            (/ shoot_frac_intpol, root_frac_intpol, fol_frac_intpol, &
                               stem_frac_intpol, DS, DR, sum(lai), fT, fP, fV, VD,   &
                               raso, max_raso /) ,                                   &
                            daily = .true.                                            )
      end if ! for crop

    end if ! if at end of the day

  end subroutine write_daily_fluxes_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_daily_water_output_csv( time )

    ! write to output files. HARD-CODED. !

    use gv_clim,                only: dayppt, wetev
    use gv_hourscale,           only: canopy_store, discharge, runoff, unintercepted
    use gv_irradiance_sunshade, only: check
    use gv_scale_declarations,  only: time_holder
    use gv_snow_info,           only: snow_watermm
    use gv_soil_structure,      only: watericemm, prevwater
    use gv_veg,                 only: soiletmm, totevap

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    double precision    :: check2, currentwater, delwater, modess, modwet 

    if (time%step .eq. time%steps_per_day) then

        ! modelled soil evaporation (mm d-1)
        modess    = sum(soiletmm)
        ! modelled evap from wetted canopy (mm d-1)
        modwet    = sum(wetev)             

        ! sanity check of water balance
        currentwater = sum( watericemm )
        delwater     = currentwater - prevwater    ! change in soil water storage (mm)
        ! total water fluxes..
        check2       = -1d0 * ( modess + runoff * 1d3 + totevap + discharge - unintercepted )
        ! check and delwater should be the same magnitude..

        call write_to_file( water_fluxes_unit ,                                       &
                          (/ modess, runoff*1d3, dayppt, totevap, delwater,         &
                             discharge, currentwater, check2, check-delwater ,      &
                             modwet, unintercepted, canopy_store, snow_watermm /) , &
                          daily = .true.                                             )
       ! save for next day
       prevwater = currentwater

    end if ! if at end of the day

  end subroutine write_daily_water_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_daily_ecosystem_fluxes_output_csv( time )

    ! write to output files. HARD-CODED. !

    use gv_clim,                only: wetev
    use gv_scale_declarations,  only: time_holder, g_to_mol_carbon, mol_to_g_water
    use gv_soil_structure,      only: weighted_SWP
    use gv_veg

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    double precision    :: daymm, daytrans, lsc, modess, modle, modwet, rplant, totalflux
    double precision,dimension(time%steps_per_day) :: posess

    ! reset sums
    if ( time%step .eq. 1 ) then
        timeflux = 0d0
    end if

    if ( time%step .eq. time%steps_per_day ) then
         modle     = sum(lemod) * time%seconds_per_step * 1d-6      ! modelled ecosystem water loss to atmosphere (MJ m-2 d-1)
         daytrans  = sum(transt) * time%seconds_per_step * 1d-6     ! modelled canopy transpiration (MJ m-2 d-1)
         modess    = sum(soiletmm)                                  ! modelled soil evaporation (mm d-1)
         posess    = max(0d0,ess)                                   ! remove dew component
         modwet    = sum(wetev)                                     ! modelled evap from wetted canopy (mm d-1)
         daymm     = sum(mmmod)                                     ! total ET (mm d-1)
         totalflux = sum(timeflux) * time%seconds_per_step * 0.001 * mol_to_g_water * 0.001  ! total flux  (mm d-1 or kg m-2 ga d-1)

         ! determine leaf specific conductance for 1st canopy layer
         if ( lai(1) .gt. 0d0 ) then
             if ( conductivity .eq. 1 ) then
                 rplant = layer_height(1) / ( gplant * lai(1) ) ! ! MPa s m2 mmol-1
             else
                 rplant = 1d0 / ( gplant * lai(1) )         ! conductance is constant with height
             end if
             lsc = ( 1d0 / ( rplant + canopy_soil_resistance(1) ) ) / lai(1)
         else
             rplant = -999d0
             lsc    = -999d0
         end if ! 

         call write_to_file( daily_unit,                                            &
                             (/ daymm, modle, daytrans, modess, modwet, weighted_SWP, rplant, &
                                canopy_soil_resistance(1), lsc, totalflux, sum(lai) /) , &
                             daily = .true.                                          )

    end if ! if at end of the day

  end subroutine write_daily_ecosystem_fluxes_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_solar_output_csv( output )

    ! moved writes from solar (light.f90) to here !

    use gv_scale_declarations, only: grid

    implicit none

    ! arguments..
    double precision, intent(in) :: output(:)

    ! local variables..
    double precision,dimension(grid%canopy) :: abspar_shade_out, abspar_sun_out, &
                                   leaf_sun_out, parabs_shade_out, parabs_sun_out

    ! output =  (/ soilnet, soilabsn, soilabsl, soilabsp, par_top, fdiff, check, skyabsp,  &
    !                      abspar_sun(grid%canopy),  &
    !                      parabs_sun(grid%canopy),  &
    !                    abspar_shade(grid%canopy),  &
    !                    parabs_shade(grid%canopy),  &
    !                        leaf_sun(grid%canopy)   /)

    abspar_sun_out   = output( 9               : 9+  grid%canopy-1 )
    parabs_sun_out   = output( 9+  grid%canopy : 9+2*grid%canopy-1 )
    abspar_shade_out = output( 9+2*grid%canopy : 9+3*grid%canopy-1 )
    parabs_shade_out = output( 9+3*grid%canopy : 9+4*grid%canopy-1 )
    leaf_sun_out     = output( 9+4*grid%canopy : 9+5*grid%canopy-1 )

    call write_to_file( solar1_unit , output(1:5) )

    call write_to_file( solar2_unit , (/ sum(abspar_sun_out), (1.-output(6))*output(5), parabs_sun_out /) )

    call write_to_file( solar3_unit , (/ sum(abspar_shade_out), output(6)*output(5) , parabs_shade_out , output(8) /) )

    call write_to_file( solar4_unit , (/ output(6), leaf_sun_out, output(7) /) )

  end subroutine write_solar_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_crops_output_csv
!
!    ! moved crops (main) output here !
!
!    use carbon_model_crop_mod, only: BM_EX, yield
!    use gv_scale_declarations, only: time_holder
!
    implicit none
!
!    call write_to_file( statistics_unit , (/ nee_cum_harvest, nee_cum, yield, BM_EX,      &
!                                             BM_EX + yield, NBP_harvest, -nee_cum - BM_EX - yield /) )
!
  end subroutine write_crops_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_soil_output_csv( daily , output )

    ! moved writes from soilday (soil_functions.f90) to here !

    use gv_hourscale,          only: freeze, hourtemp, overflow, Qm, Qs, Qc, Qe, Qh, Qn, surface_watermm, &
                                     totet, underflow
    use gv_scale_declarations, only: grid
    use gv_snow_info,          only: Dsnow, Nsnow, Sice, Sliq, snowheight, snowweight, Tsnow
    use gv_soil_structure,     only: drythick, fraction_uptake, iceprop, pptgain, soil_temp, waterfrac, &
                                     watergain, watericemm, waterloss, weighted_SWP, SWP
    use log_tools

    implicit none

    ! arguments..
    logical,intent(in) :: daily
    double precision,   intent(in) :: output(:)

    ! locatl variables
    integer :: i    
    double precision, dimension(grid%snow) :: snow_layer_height, snow_temp

    if (.not.daily) then
        call write_to_file( energy_unit , (/ Qh, Qe, Qn, Qc, Qm, Qs, hourtemp-freeze, &
                                            output(1)-freeze, soil_temp(2)-freeze, drythick*1d3/) )

        call write_to_file( ice_prop_unit , iceprop(1:15) )

        call write_to_file( soil_status_unit, (/ 1d3*overflow, 1d3*totet, 1d3*underflow,             &
                                              surface_watermm, sum(watericemm), output(2:4),         &
                                              1d3*watergain(1), 1d3*waterloss(1), output(5:6),       &
                                              sum(pptgain), sum(watergain), sum(waterloss) ,         &
                                              watergain(grid%core), waterloss(grid%core), snowweight, snowheight /) )

        call write_to_file( soil_temp_unit , (/soil_temp(1:15)-273.15d0, soil_temp(grid%core)-273.15d0/) )

        call write_to_file( soil_water_unit , (/ waterfrac(1:15), weighted_swp, SWP(1:10) /) )
 
        snow_temp(:) = -9999d0
        do i = 1, Nsnow
           snow_temp(i) = Tsnow(i) - 273.15d0
        enddo
        snow_layer_height(:) = 0d0
        if (Nsnow > 0) then
            snow_layer_height(Nsnow) = 0.5d0*Dsnow(Nsnow)
            do i = Nsnow-1, 1, -1
               snow_layer_height(i) = snow_layer_height(i+1) + 0.5d0*(Dsnow(i) + Dsnow(i+1))
            enddo
        endif
        call write_to_file( snow_temp_unit , (/snow_temp(1:5), snow_layer_height(1:5)/) )
        call write_to_file( snow_mass_unit , (/Sice(1:5), Sliq(1:5)/) )

    elseif (daily) then

        call write_to_file( up_frac_unit , fraction_uptake(1:15) )

    else

        write(message,*)"flag supplied to handle_output:",daily," was not recognised"
        call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )

    endif ! daily output or not

  end subroutine write_soil_output_csv
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine open_file( filename , unitnos , readonly , header , daily )
 
    ! This s/r is a wrapper to open(). It performs !
    ! sanity checks when opening files.  Control   !
    ! of write-permission is available through the !
    ! optional readonly argument, and if header is !
    ! provided it will be written to the file.     !

    use log_tools

    implicit none

    ! arguments..
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: unitnos
    logical,optional, intent(in) :: readonly, daily
    character(len=*),optional, &
                      intent(in) :: header

    ! local variables..
    integer            :: ios
    logical            :: print_day_msg, read_mode
    character(len=6)   :: header_fmt
    character(len=20)  :: file_status, write_status
    character(len=800) :: full_header

    ios = 0

    ! determine whether reading or writing..
    ! (assume that it is okay to write)
    read_mode    = .false.
    file_status  = 'replace'
    write_status = 'readwrite'
    if ( present( readonly ) ) then
      if ( readonly ) then
        read_mode    = .true.
        file_status  = 'old'
        write_status = 'read'
      end if
    end if

    ! determine if output data will be daily-avg'd..
    ! (assume that msg is not wanted)
    print_day_msg = .false.
    if ( present( daily ) ) print_day_msg = daily

    ! open file..
    open( unit=unitnos , file=trim(filename) , iostat=ios ,    &
          status=trim(file_status) , action=trim(write_status) )

    ! check opened file okay..
    if ( ios .ne. 0 ) then
      write(message,*)"Problem opening file : ",trim(filename)
      call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    end if

    ! write header if provided..
    if ( present( header ) ) then
      if ( read_mode ) then
        write(message,*)"You cannot write a header to a file that "&
                      //"you have opened in read-only mode!"
        call write_log( message , msg_warning , __FILE__ , __LINE__ )
      else
        if ( print_day_msg ) then
          write(unitnos,*)'!NB. This file contains data that is averaged or summed over daily periods.!'
        end if
        full_header='Time (days),'//trim(header)
        write(header_fmt,'(A2,I3.3,A1)')'(A',len(trim(full_header)),')'
        write(unitnos,header_fmt)trim(full_header)
      end if
    end if

  end subroutine open_file
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
  subroutine write_to_file( unit_nos , out_data , daily )

    use gv_scale_declarations, only: time
    use log_tools

    implicit none

    ! arguments..
    integer,intent(in)          :: unit_nos
    double precision,   intent(in)          :: out_data(:)
    logical,optional,intent(in) :: daily

    ! local variables..
    logical           :: daily_data, file_open
    integer           :: n
    character(len=7)  :: write_permit
    character(len=20) :: write_fmt

    if (present(daily)) then
        daily_data = daily
    else 
        daily_data = .false.
    end if

    ! check file is open and can be written to..
    inquire( unit=unit_nos, opened=file_open, write=write_permit )

    if ( .not. file_open ) then

        ! Unit-nos is not for a file that is open!!
        write(message,*)"Unit_nos:",unit_nos," is not connected to an open file!"
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    else if ( trim(write_permit) .ne. 'YES' ) then

        ! Unit-nos is for a file in read-only setting!!
        write(message,*)"Unit number:",unit_nos," is not open in write-mode!"
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    else

        ! Okay to do the write..
        n = size( out_data )
        if ( daily_data ) then
            write(write_fmt,'(a,i0,a)') "(i0.3," , n ,'(",",f0.7))'
            write( unit=unit_nos , fmt=trim(write_fmt) )  ((time%year-1)*time%days_in_year(time%year)+time%day) , out_data
        else
            write(write_fmt,'(a,i0,a)') "(f0.2," , n ,'(",",f0.7))'
            write( unit=unit_nos , fmt=trim(write_fmt) )  time%daytime , out_data
        end if

    end if

  end subroutine write_to_file
  !
  !--------------------------------------------------------------------------------------------------------------------------------!
  !
end module spa_io_csv
!
!----------------------------------------------------------------------------------------------------------------------------------!
!
