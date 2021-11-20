

subroutine rspatrunk(output_dim,met,pheno,pars,out_var,lat,long,nopheno &
                    ,nopars,nomet,nofluxes,nopools,pft,nodays,deltat &
                    ,nos_iter,nosteps,out_steps,year_in)

  use main_spa, only: spa

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  implicit none
  ! declare input variables
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,year_in        & !
                        ,out_steps      & !
                        ,output_dim     & !
                        ,pft            & ! plant functional type
                        ,nos_iter       & !
                        ,nosteps        & ! model steps for each simulation
                        ,nopheno        & ! number of phenology fields
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(in) :: met(nomet,nosteps) & ! met drivers, note reverse of needed
                        ,pheno(nodays,nopheno)         &
                        ,deltat                        & ! time step in decimal days
                        ,pars(nopars)                  & ! number of parameters
                        ,long                          &
                        ,lat                 ! site latitude (degrees)

  ! output declaration
  double precision, intent(out), dimension(nos_iter,out_steps,output_dim) :: out_var

  ! local variables
  ! vector of ecosystem pools
  double precision, dimension(out_steps,nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(out_steps,nofluxes) :: FLUXES
  integer i, steps_per_day
  double precision, dimension(out_steps) :: lai &        ! leaf area index
                                         ,daily_weighted_SWP   &
                                         ,daily_weighted_soilR &
                                         ,sun_canopy_temperature &
                                         ,shade_canopy_temperature &
                                         ,canopy_temperature &
                                         ,soil_conductance &
                                         ,ground_heat &
                                         ,wetcanopyevap &
                                         ,potevap    &
                                         ,netrad     &
                                         ,porosity   &
                                         ,soilevap   &
                                         ,evapotrans & ! evapotranspiration (W.m-2/mm.day-1)
                                         ,sensible   & ! sensible heat (W.m-2 / MJ.m-2.day-1)
                                         ,soilwater  & ! soil water content (mm)
                                         ,GPP &        ! Gross primary productivity (umolC.m-2.s-1/gC.m-2.day-1)
                                         ,NEE          ! net ecosystem exchange of CO2 (umolC.m-2.s-1/gC.m-2.day-1)
! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0,total_time = 0
!real :: Rtot_track_time=0, aero_time=0 , soilwater_time=0 , acm_et_time = 0 , Rm_time = 0
  ! zero initial conditions
  lai = 0d0 ; GPP = 0d0 ; NEE = 0d0 ; POOLS = 0d0 ; FLUXES = 0d0 ; out_var = 0d0

  ! steps_per_day
  steps_per_day = int(deltat**(-1d0))
!call cpu_time(begin)
  ! begin iterations
  do i = 1, nos_iter

     ! call the models
     call spa(met,pheno(1:nodays,1:nopheno),pars(1:nopars),deltat,nosteps &
             ,out_steps,nodays,lat,long,year_in,lai,GPP,FLUXES,POOLS,nopheno &
             ,nopars,nomet,nopools,nofluxes,soilwater,evapotrans,sensible &
             ,daily_weighted_SWP,daily_weighted_soilR,soilevap,wetcanopyevap &
             ,potevap,netrad,porosity,sun_canopy_temperature,shade_canopy_temperature &
             ,canopy_temperature,ground_heat,soil_conductance)

     ! now allocate the output the our 'output' variable
     out_var(i,1:out_steps,1)  = lai(1:out_steps)
     out_var(i,1:out_steps,2)  = GPP(1:out_steps)
     out_var(i,1:out_steps,3)  = pars(1)
     out_var(i,1:out_steps,4)  = sensible(1:out_steps)
     out_var(i,1:out_steps,5)  = FLUXES(1:out_steps,1)
     out_var(i,1:out_steps,6)  = soilwater(1:out_steps)
     out_var(i,1:out_steps,7)  = evapotrans(1:out_steps)
     out_var(i,1:out_steps,8)  = daily_weighted_SWP(1:out_steps)
     out_var(i,1:out_steps,9)  = daily_weighted_soilR(1:out_steps)
     out_var(i,1:out_steps,10) = soilevap(1:out_steps)
     out_var(i,1:out_steps,11) = wetcanopyevap(1:out_steps)
     out_var(i,1:out_steps,12) = potevap(1:out_steps)
     out_var(i,1:out_steps,13) = netrad(1:out_steps)
     out_var(i,1:out_steps,14) = porosity(1:out_steps)
     out_var(i,1:out_steps,15) = sun_canopy_temperature(1:out_steps)
     out_var(i,1:out_steps,16) = shade_canopy_temperature(1:out_steps)
     out_var(i,1:out_steps,17) = canopy_temperature(1:out_steps)
     out_var(i,1:out_steps,18) = ground_heat(1:out_steps)
     out_var(i,1:out_steps,19) = soil_conductance(1:out_steps)

  end do ! nos_iter loop
!  call cpu_time(done)
!  print*,"time taken per step",(done-begin) / real(nodays), nodays

  ! return back to the subroutine then
  return

end subroutine rspatrunk
