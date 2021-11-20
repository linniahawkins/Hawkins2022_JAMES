! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2014 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module math_tools

  !! This module contains general maths routines necessary for !!
  !! calculating derivatives and bisections.  See individual   !!
  !! routines for details.                                     !!

  use gv_scale_declarations, only: max_nos_iterations, pi, deg_to_rad

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: tridiag, interpolate, ode_int, zbrent, quadratic, &
            calculate_daylength_hours
  ! Variables..
  public :: dxsav, kmax

  ! variables shared across SPA..
  integer :: kmax   !
  double precision    :: dxsav  !

  ! variables private to this module..
  integer, parameter :: kmaxx  = 200,  & ! descriptions
                       maxstp = 10000   !   would
  double precision,parameter    :: tiny   = 1.d-30  !        nice!

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  double precision function calculate_daylength_hours(doy,latitude_in_radians)

    implicit none

    ! arguments
    double precision, intent(in) :: doy, latitude_in_radians

    ! local arguments
    double precision :: declination,sinld,cosld,aob

    ! calculate daylength in hours..(warning declination is a module variable)
    declination  = - asin ( sin ( 23.45d0 * deg_to_rad ) * cos ( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
    sinld     = sin ( latitude_in_radians ) * sin ( declination )
    cosld     = cos ( latitude_in_radians ) * cos ( declination )
    aob = max(-1d0,min(1d0,sinld / cosld))
    calculate_daylength_hours = 12d0 * ( 1d0 + 2d0 * asin ( aob ) / pi )

    return

  end function calculate_daylength_hours
  !
  !----------------------------------------------------------------------
  !
  double precision function interpolate( x , reference_x , reference_y , row )

    ! Interpolation function.                    !
    ! x is input value, interpol is output value !
    ! reference_x/y are reference input data.    !

    implicit none

    ! arguments..
    integer,intent(in)             :: row
    double precision,intent(in)                :: x
    double precision,dimension(row),intent(in) :: reference_x , reference_y

    ! local variables..
    integer::i

    do i = 1 , row

       if ( x .le. reference_x(1) ) then
           interpolate = reference_y(1)
           exit
       end if

       ! cycling means growth rate remains constant between DS levels
       if ( ( x .gt. reference_x(i) ) .and. ( i .lt. row ) ) cycle

       if ( x .eq. reference_x(i) ) then 
           interpolate = reference_y(i) 
           exit
       end if

       if ( x .lt. reference_x(i) ) then
           interpolate = reference_y(i-1) + ( x - reference_x(i-1) ) &
                        * ( reference_y(i) - reference_y(i-1) )      &
                        / ( reference_x(i) - reference_x(i-1) )
           exit
       else
           interpolate = reference_y(row)
       end if

    enddo

  end function interpolate
  !
  !----------------------------------------------------------------------
  !
  subroutine ode_int( called_from , ystart , nvar , x1 , x2 , eps , h1 , hmin , nok , nbad , derivs )

    ! This is an integrator for ordinary differential equations. We use !
    ! the RUNGE_KUTTA to track the dynamic behaviour of various model   !
    ! state variables.. uch a leaf water potential, and water stored on !
    ! leaf surfaces. The RUNGE_KUTTA finds the time step that ensures   !
    ! dynamics are smooth and that the relevant feedbacks are properly  !
    ! incorporated. For a full description see Press et al. (1986).     !

    use gv_scale_declarations, only: max_nos_iterations
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    integer,intent(in)          :: nvar 
    double precision,intent(in)             :: h1, hmin, x1, x2
    double precision,intent(inout)          :: ystart(nvar), eps
    integer,intent(out)         :: nbad, nok

    ! Interfaces are the correct way to pass procedures as arguments.
    ! (note 'derivs' is actually one of canopy_water_store,soil_water_store,lwp_diff_eqn)
!    external :: derivs
    interface
       subroutine derivs(time,y,dydt)
         use gv_scale_declarations, only: max_nos_iterations
         double precision,intent(in)    :: time
         double precision,intent(in)    :: y(max_nos_iterations)
         double precision,intent(out)   :: dydt(max_nos_iterations)
       end subroutine derivs
    end interface

    ! local variables..
    integer :: i, kount, nstp
    double precision    :: h, hdid, hnext, x, xsav ,         &
               dydx(max_nos_iterations), y(max_nos_iterations), yscal(max_nos_iterations), &
               xp(kmaxx), yp(max_nos_iterations,kmaxx) 

    ! calculations..
    x = x1
    h = sign( h1 , x2-x1 )
    nok = 0
    nbad = 0
    kount = 0
    do i = 1 , nvar
       y(i) = ystart(i)
    enddo
    if ( kmax .gt. 0 ) xsav = x - 2d0 * dxsav
    do nstp = 1 , MAXSTP
       call derivs( x , y , dydx )
       do i = 1 , nvar
          yscal(i) = abs(y(i))+abs(h*dydx(i))+TINY
       enddo
       if ( kmax .gt. 0 ) then
          if ( abs( x - xsav ) .gt. abs( dxsav ) ) then
             if ( kount .lt. kmax - 1 ) then
                kount = kount + 1
                xp(kount) = x
                do i = 1 , nvar
                   yp(i,kount) = y(i)
                enddo
                xsav = x
             end if
          end if
       end if
       if ( (x+h-x2) * (x+h-x1) .gt. 0d0 ) h = x2 - x

       call runge_kutta_q_step( trim(called_from)//":ode_int" , y , dydx , nvar , x , h , eps , yscal , hdid , hnext , derivs )
       if (hdid.eq.h) then
          nok = nok+1
       else
          nbad = nbad+1
       end if
       if ( ( x-x2 ) * ( x2-x1 ) .ge. 0d0 ) then
          do i = 1 , nvar
            ystart(i) = y(i)
          enddo
          if ( kmax .ne. 0 ) then
            kount = kount + 1
            xp(kount) = x
            do i = 1 , nvar
              yp(i,kount) = y(i)
            enddo
          end if
          return
       end if
       if ( abs(hnext) .lt. hmin ) then
         write(message,*) "stepsize smaller than permitted minimum in ode_int",new_line('x'),&
                          "ode_int was called by: ",trim(called_from)
         call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
       end if
       h  =  hnext
    enddo

  end subroutine ode_int
  !
  !----------------------------------------------------------------------
  !
  subroutine tridiag( Nvec, Nmax, a, b, c, r, x )
    ! Tridiagonal matrix solver !

    implicit none

    ! arguments..
    integer, intent(in) :: Nvec, Nmax
    double precision, intent(in)    :: a(Nmax), b(Nmax), c(Nmax), r(Nmax)
    double precision, intent(out)   :: x(Nmax)

    ! local variables..
    integer :: n
    double precision    :: beta, gamma(Nvec)

    beta = b(1)
    x(1) = r(1) / beta

    do n = 2, Nvec
      gamma(n) = c(n-1) / beta
      beta = b(n) - a(n)*gamma(n)
      x(n) = (r(n) - a(n)*x(n-1)) / beta
    end do

    do n = Nvec - 1, 1, -1
      x(n) = x(n) - gamma(n+1)*x(n+1)
    end do

  end subroutine tridiag
  !
  !----------------------------------------------------------------------
  !
  double precision function zbrent( called_from , func , x1 , x2 , tol )

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    !  reference to a particular function and also two values which bound !
    !  the arguments for the function of interest. ZBRENT finds a root of !
    !  the function (i.e. the point where the function equals zero), that !
    !  lies between the two bounds.                                       !
    ! For a full description see Press et al. (1986).                     !

    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    double precision,intent(in)             :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision,intent(in) :: xval
       end function func
    end interface

!   external :: derivs
!    interface
!       subroutine derivs(time,y,dydt)
!         use gv_scale_declarations, only: max_nos_iterations
!         double precision,intent(in)    :: time
!         double precision,intent(in)    :: y(max_nos_iterations)
!         double precision,intent(out)   :: dydt(max_nos_iterations)
!       end subroutine derivs
!    end interface


    ! local variables..
    integer            :: iter
    integer,parameter  :: ITMAX = 30
    double precision               :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision,parameter     :: EPS = 3.d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0d0 ) then
        zbrent = a
        return
    elseif ( fb .eq. 0d0 ) then
        zbrent = b
        return
    end if

    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
    if ( sign(1d0,fa) .eq. sign(1d0,fb) ) then
       fa = func( a )
       fb = func( b )
       if (trim(called_from) == "leaf_temperature_hfx_sun:leaf_balance" .or. &
           trim(called_from) == "leaf_temperature_hfx_shade:leaf_balance" .or. &
           trim(called_from) == "leaf_temperature_hfx_dead:leaf_balance" .or. &
           trim(called_from) == "assimilate:leaf_balance") then
           ! do nothing if one of these   
           else
           ! tell me otherwise what is going on
           write(message,*)"Supplied values must bracket the root of the function.",new_line('x'),  &
             "     ","You supplied x1:",x1,new_line('x'),                     &
             "     "," and x2:",x2,new_line('x'),                             &
             "     "," which give function values of fa :",fa,new_line('x'),  &
             "     "," and fb:",fb," .",new_line('x'),                        &
             " zbrent was called by: ",trim(called_from)
           call write_log( trim(message) , msg_error , __FILE__ , __LINE__ )
       end if
       fa = func( a )
       fb = func( b )
    end if
    c = b
    fc = fb
    do iter = 1 , ITMAX

       ! If the new value (f(c)) doesn't bracket
       ! the root with f(b) then adjust it.. 
       if ( sign(1d0,fb) .eq. sign(1d0,fc) ) then
           c  = a
           fc = fa
           d  = b - a
           e  = d
       end if
       if ( abs(fc) .lt. abs(fb) ) then
           a  = b
           b  = c
           c  = a
           fa = fb
           fb = fc
           fc = fa
       end if
       tol1 = 2d0 * EPS * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
           zbrent = b
           return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
           s = fb / fa
           if ( a .eq. c ) then
              p = 2d0 * xm * s
              q = 1d0 - s
           else
              q = fa / fc
              r = fb / fc
              p = s * ( 2d0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1d0 ) )
              q = ( q - 1d0 ) * ( r - 1d0 ) * ( s - 1d0 )
           end if
           if ( p .gt. 0d0 ) q = -q
           p = abs( p )
           if ( (2d0*p) .lt. min( 3d0*xm*q-abs(tol1*q) , abs(e*q) ) ) then
              e = d
              d = p / q
           else
              d = xm
               e = d
           end if
       else
           d = xm
           e = d
       end if
       a  = b
       fa = fb
       if ( abs(d) .gt. tol1 ) then
          b = b + d
       else
          b = b + sign( tol1 , xm )
       end if
       fb = func(b)
    enddo
    write(message,*) "zbrent has exceeded maximum iterations",new_line('x'),&
                     "zbrent was called by: ",trim(called_from)
    call write_log( message , msg_warning , __FILE__ , __LINE__ )
    zbrent = b

  end function zbrent

  !
  !----------------------------------------------------------------------
  !
  subroutine quadratic (a, b, c, r1, r2)
    !
    ! DESCRIPTION:
    ! Solve a quadratic equation for its two roots
    !

    implicit none

    double precision, intent(in)        :: a,b,c        ! Terms for quadratic equation
    double precision, intent(out)       :: r1,r2        ! Roots of quadratic

    ! Local variables..
    double precision :: q ! Temporary term

    if (b >= 0d0) then
        q = -0.5d0 * (b + sqrt(b*b - 4d0*a*c))
    else
        q = -0.5d0 * (b - sqrt(b*b - 4d0*a*c))
    end if

    r1 = q/a
    if (q .eq. 0d0) then
        r2 = 1d0**36d0
    else
        r2 = c / q
    end if

  end subroutine quadratic

  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  double precision function bilinear_interp( data , corners , location )

    ! Bilinearly interpolates to find value at location based  !
    ! on values in data occuring at corners.  Assumes that     !
    ! data  starting at bottom left and !
    ! moving clockwise,ie 2nd is top left, 3rd is top right.   !

    implicit none

    ! arguments..
    double precision,dimension(4),  intent(in) :: data     ! the values of the 4 nearest grid pts
    double precision,dimension(2,2),intent(in) :: corners  ! (1,*) = x1,x2, (2,*) = y1,y2
    double precision,dimension(2),  intent(in) :: location ! (x,y)

    ! local variables..,
    double precision :: xfrac, yfrac

    ! check that there is not almost-zero x distance...
    if ( ( corners(1,2) - corners(1,1) ) .lt. tiny ) then
       xfrac = 1
    else
       xfrac = ( location(1) - corners(1,1) ) / ( corners(1,2) - corners(1,1) )
    end if

    ! check that there is not almost-zero y distance...
    if ( ( corners(2,2) - corners(2,1) ) .lt. tiny ) then
       yfrac = 1
    else
       yfrac = ( location(2) - corners(2,1) ) / ( corners(2,2) - corners(2,1) )
    end if

    ! calculate the solution...
    bilinear_interp =     yfrac     *     xfrac     * data(3)  &
                    +     yfrac     * ( 1 - xfrac ) * data(2)  &
                    + ( 1 - yfrac ) *     xfrac     * data(4)  &
                    + ( 1 - yfrac ) * ( 1 - xfrac ) * data(1)

  end function bilinear_interp
  !
  !----------------------------------------------------------------------
  !
  subroutine runge_kutta_check( y , dydx , n , x , h , yout , yerr , derivs )

    ! > subroutine summary? < !

    use gv_scale_declarations, only: max_nos_iterations

    implicit none

    ! arguments..
    integer,intent(in) :: n
    double precision,intent(in)    :: h,x,dydx(n),y(n)
    double precision,intent(out)   :: yerr(n),yout(n)

    ! Interfaces are the correct way to pass procedures as arguments.
    ! (note 'derivs' is actually one of canopy_water_store,soil_water_store,lwp_diff_eqn)
    external :: derivs
    interface
       subroutine derivs( time , y , dydt )
         use gv_scale_declarations, only: max_nos_iterations
         double precision,intent(in)   :: y(max_nos_iterations)
         double precision,intent(in)   :: time
         double precision,intent(out)  :: dydt(max_nos_iterations)
       end subroutine derivs
    end interface

    ! local variables..
    integer        :: i
    double precision           :: ak2(max_nos_iterations),ak3(max_nos_iterations) &
                     ,ak4(max_nos_iterations),ak5(max_nos_iterations) &
                     ,ak6(max_nos_iterations),ytemp(max_nos_iterations)
    double precision, parameter :: A2 = 0.2d0, A3 = 0.3d0, A4 = 0.6d0, A5 = 1d0,    &
         A6 = 0.875d0, B21 = 0.2d0, B31 = 3d0/40d0, B32 = 9d0/40d0, B41 = 0.3d0,    &
         B42 = -0.9d0, B43 = 1.2d0, B51 = -11d0/54d0, B52 = 2.5d0, B53 = -70d0/27d0,&
         B54 = 35d0/27d0, B61 = 1631d0/55296d0, B62 = 175d0/512d0,                  &
         B63 = 575d0/13824d0, B64 = 44275d0/110592d0, B65 = 253d0/4096d0,           &
         C1 = 37d0/378d0, C3 = 250d0/621d0, C4 = 125d0/594d0, C6 = 512d0/1771d0,    &
         DC1 = C1-2825d0/27648d0, DC3 = C3-18575d0/48384d0,                         &
         DC4 = C4-13525d0/55296d0, DC5 = -277d0/14336d0, DC6 = C6-0.25d0

    ! calculations...
    do i = 1 , n
       ytemp(i) = y(i) + B21 * h * dydx(i)
    enddo
    call derivs( x + A2 * h , ytemp , ak2 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B31 * dydx(i) + B32 * ak2(i) )
    enddo
    call derivs( x + A3 * h , ytemp , ak3 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i) )
    enddo
    call derivs( x + A4 * h , ytemp , ak4 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + B54 * ak4(i) )
    enddo
    call derivs( x + A5 * h , ytemp , ak5 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + B64 * ak4(i) + B65 * ak5(i) )
    enddo
    call derivs( x + A6 * h , ytemp , ak6 )
    do i = 1 , n
       yout(i) = y(i) + h * ( C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + C6 * ak6(i) )
       yerr(i) = h * ( DC1 * dydx(i) + DC3 * ak3(i) + DC4 * ak4(i) + DC5 * ak5(i) + DC6 * ak6(i) )
    enddo

  end subroutine runge_kutta_check
  !
  !----------------------------------------------------------------------
  !
  subroutine runge_kutta_q_step( called_from , y , dydx , n , x , htry , eps , yscal , hdid , hnext , derivs )

    ! > subroutine summary? < !

    use gv_scale_declarations, only: max_nos_iterations
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    integer,intent(in)          :: n
    double precision,intent(in)             :: eps, htry, dydx(n), yscal(n)
    double precision,intent(inout)          :: x
    double precision,intent(out)            :: hdid, hnext, y(n)

    ! Interfaces are the correct way to pass procedures as arguments.
    ! (note 'derivs' is actually one of canopy_water_store,soil_water_store,lwp_diff_eqn)
    interface
       subroutine derivs( time , y , dydt )
         use gv_scale_declarations, only: max_nos_iterations
         double precision,intent(in)  :: y(max_nos_iterations)
         double precision,intent(in)  :: time
         double precision,intent(out) :: dydt(max_nos_iterations)
       end subroutine derivs
    end interface

    ! local variables..
    integer        :: i
    double precision            :: errmax, h, htemp, xnew, yerr(max_nos_iterations), ytemp(max_nos_iterations)
    double precision, parameter :: ERRCON = 1.89d-4, & !
                                   PGROW  = -0.2d0,  & !
                                   PSHRNK = -0.25d0, & !
                                   SAFETY = 0.9d0      !

    ! calculations...
    h = htry
1   call runge_kutta_check( y , dydx , n , x , h , ytemp , yerr , derivs )
    errmax = 0d0
    do  i = 1 , n
       errmax = max( errmax , abs( yerr(i) / yscal(i) ) )
    enddo
    errmax = errmax / eps
    if ( errmax .gt. 1d0 ) then
       htemp = SAFETY * h * ( errmax**PSHRNK )
       h = sign( max( abs(htemp) , 0.1d0*abs(h) ) , h )
       xnew = x + h
       if ( xnew .eq. x ) then
         write(message,*) "stepsize underflow in runge_kutta_q_step",new_line('x'),&
                          "runge_kutta_q_step called from: ",trim(called_from)
         call write_log( message , msg_warning , __FILE__ , __LINE__ )
       end if
       goto 1
    else
       if ( errmax .gt. ERRCON ) then
          hnext = SAFETY * h * ( errmax**PGROW )
       else
          hnext = 5d0 * h
       end if
       hdid = h
       x = x + h
       do i = 1 , n
          y(i) = ytemp(i)
       enddo

    end if

  end subroutine runge_kutta_q_step
  !
  !----------------------------------------------------------------------
  !
end module math_tools
!
!------------------------------------------------------------------------
!
