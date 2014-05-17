PROGRAM DISCRETE
   use real16
    use global_variables
    use discrete_helpers
    use DVODE_F90_M
    implicit none

! Population balace differential equations are solved in the form
! dN_i / dt = F_i

! general syntax: using caps only for global variables!
! declaring the basic variables

    real(real_16)          ::     t,          &  ! time
                                  dt_r16,    &    ! time step in r16
                                  tmax,        &     ! end of execution time
                                  ntot,        &    ! total cluster number in system
                                  vtot,        &     ! total molecule number in system
                                  temp_s             ! temperature in system

    real(real_16), dimension(:),allocatable    :: N_main   ! the

    real(real_16)            :: big_clusters,    &    ! number of clusters grown out of the system
                                     big_flux                 !  rate of clusters moving out of the system (per dt)

    integer                     :: m1, m2  ! molecule indexes




! here start the declarations of less important variables

    INTEGER                    :: alloc_error,dummy  ! error message for array allocation

! dvode specific variables, 
! THESE ARE ARRAYS!
    double precision, dimension(:),allocatable :: atol_main, &
                                                 rtol_main, &
                                                 y0_main  ! dimension = MAXNEQ

    double precision :: t0,  &          ! start time
                       dt,  &          ! time step
                       dtd, &          ! time step change factor 
                       atol_ini,  &    ! absolute tolerance
                       rtol_ini        ! relative tolerance


call set_globals




!! allocating memory...
write(*,'(a)',advance='no'),'Allocating memory for concentration arrays...'
allocate(N_main(0:IMAX),y0_main(CURNEQ),atol_main(CURNEQ),rtol_main(CURNEQ), COAGSINK(CURNEQ),&
             STAT = alloc_error)
print *,'done. Status=',alloc_error

! setting the collision (and evaporation) coefficients
! also allocates COE, ECOE and IN_USE
call init_coeff_table

print *, DRYDIAM(50)
print *, WETDIAM(50)

print *, 'SA mol vol:', molecvol_h2so4
print *, 'SA mol radius:', (molecvol_h2so4*3./(4.*pi))**(1./3.)

call make_coag_sink

! make the initial distribution into y0
N_main(:) = 0.01
N_main(NUCSIZE) = 0.0e12
N_main(1) = CVAP_0
y0_main = N_main(1:IMAX)

atol_ini = 10.0
rtol_ini = 0.1

atol_main(:)=atol_ini
rtol_main(:)=rtol_ini

print *,'MAXNEQ: ', MAXNEQ, ' CURNEQ:', CURNEQ

! still to do: take care of stuff going out of the system....
! y0_main(MAXNEQ) = big_clusters
! atol_main(MAXNEQ)= atol_ini
! rtol_main(MAXNEQ)= rtol_ini


!! setting the start time
t0 = 0.0
dt = 10.0
dtd= 1.0
tmax = 4.0*3600.
print *, 'Starting dvode...'

call int_dvode_f90(t0,dt,dtd,tmax,y0_main,rtol_ini,atol_main)
print *,'Getting out of int_dvode...deallocating...'

deallocate(COE,ECOE)
! deallocate(J)

if (allocated(y0_main)) then
    print*,'Deallocating y0_main...'
     deallocate(y0_main,STAT = alloc_error)
    print *,'Done. STATUS = ', alloc_error
else
    print*,'y0_main already deallocated!'
end if

if (allocated(rtol_main)) then
 print*,'Deallocating rtol_main...'
 deallocate(rtol_main,STAT = alloc_error)
    print *,'Done. STATUS = ', alloc_error
else
    print*,'rtol_main already deallocated!'
end if

if (allocated(atol_main)) then
 print*,'Deallocating atol_main...'
 deallocate(atol_main,STAT = alloc_error)
    print *,'Done. STATUS = ', alloc_error
else
    print*,'atol_main already deallocated!'
end if

if (allocated(N_main)) then
print*,'Deallocating N_main...'
 deallocate(N_main,STAT = alloc_error)
    print *,'Done. STATUS = ', alloc_error
else
    print*,'N_main already deallocated!'
end if

END program DISCRETE


!! *** THE MAIN INTEGRATION LOOP   ************************
subroutine int_dvode_f90(t0,dt_in,dtd_in,tmax,y0,rtol,atol)
  use global_variables
  use DVODE_F90_M  ! the solver module
  use DISCRETE_DVODE_F90 ! the diff equation declarations
  implicit none


    double precision, intent(in) :: t0, dt_in, dtd_in,tmax
    double precision             :: dt


    double precision, dimension(CURNEQ) :: y0, atol ! dimension = NEQ
    double precision                    :: rtol
    intent(in)                          :: y0
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.

    double precision, dimension(CURNEQ)  :: y  ! dimension = NEQ

    double precision :: tim         !  The initial value of the independent variable.
    double precision :: tout        !  First point where output is desired (.ne. T)
    double precision :: rpar        !  RPAR,IPAR = user-defined real and integer arrays passed to F and JAC.
    integer :: ipar                 !  ...
    integer :: itask                !  1 for normal computation of output values of Y at t = TOUT.
    integer :: istate               !  Integer flag (input and output).  Set ISTATE = 1.
    integer :: itol                 !  1 or 2 according as ATOL (above) is a scalar or array.
    integer :: iopt                 !  IOPT   = 0 to indicate no optional input used.
    integer :: iout                 !
    integer :: i,ww,aa,ww2,aa2,ix   ! for a loop

    TYPE (VODE_OPTS) :: OPTIONS

    double precision :: test
	character(25)             :: form ! for saving in dynamic format

! set up variables

  dt = dt_in
  y = y0
  tim = t0
  itol = 2
  itask = 1
  istate = 1
  iopt = 0
  tout = dt

OPTIONS = SET_OPTS(RELERR=rtol, ABSERR_VECTOR=atol)


!pause
! start calculating....
    do iout = 1,STEPMAX

        call dvode_f90(f_kine_vode_f90,CURNEQ,y,tim,tout,itask,istate,OPTIONS)

        print *,'Got out of Dvode now'

!! Stopping calculation because of error: 
        if (istate .lt. 0) then
            write(6,'(a,i4)') 'Error halt: istate = ',istate
!            print *, 'Distribution: ', y

            return
        end if

        print *,'t =',tim
        print *,y(NUCSIZE:NUCSIZE+10)

! test saving 
        open(unit=10,file='DISC_TEST.TXT',status='unknown',position='append')
  		form = '(d14.6,1x,XXXXXD14.6)'
		write(form(11:15),"(I5.5)"),IMAX
		print *,"Save format string = "
		print *,form     
        write(10,form),tim,y(1:IMAX)
        close(10)
    !pause
        if (tout .gt. tmax) then
            print *,'Calculation time over, stopping...'
            return
        end if
        tout = tout + dt


    end do

end subroutine int_dvode_f90

subroutine set_globals
use global_variables
implicit none
INTEGER :: alloc_err, ios
character(79) :: dummytxt

print *, 'Defaulting globals...'

! setting the defaults...
  IMAX = 1000
  TEMP = 273.15  ! K
  RH = 0.75      ! 0-1
  PRES = PSTAND  ! 
  NUCSIZE = 8    ! the size (in No. of H2SO4 molecules in "nucleating" particle
  NUCRATE = 1000.*1.0d6 ! The nucleation rate   
  PULSE_LENGTH = 900.0  ! length of nucleation pulse in seconds
  COND_ON = 1           ! condensation turned on
  EVAP_ON = 0           ! evaporation turned on
  SINK_ON = 1                ! external coagulation sink turned on
  COAG_ON = 0                 ! intermodal coagulation turned on
  NUC_MECH =1                ! the nucleation mechanism used (1 = const, 2 = act 3 = kine, 4 = free 5 = org)
  NUC_COEFF=1.0d-14     ! the nucleation coefficient (if NUC_MECH>1)
  NUC_EXP=2.0                 ! the nucleation exponent (if NUC_MECH>3)
  NUC_COEFF_ORG = 1.0d-14 ! the nucleation coefficient for organics (if NUC_MECH=5)
  NUC_EXP_ORG = 2.0       ! the nucleation exponent for organics (if NUC_MECH=5)
  CVAP_0 = 1.00d7*1e6     ! the (initial) sulphuric acid concentration
  QVAP_0 = 1.00d5*1e6     ! the vapour source rate
  CONST_CVAP = 1;
  

print *, 'Reading globals...'
! reading the globals from a file
  open(unit=10, file = 'MODEL_SETUP.TXT' ,status = 'old', action='read', position = 'rewind',iostat=ios)


  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),IMAX
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),TEMP
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),RH
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),PRES
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),NUCSIZE
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),NUCRATE
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),PULSE_LENGTH
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),COND_ON
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),EVAP_ON
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),SINK_ON
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),COAG_ON
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),NUC_MECH
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),NUC_COEFF
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),NUC_EXP
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),NUC_COEFF_ORG
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),NUC_EXP_ORG
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),CVAP_0
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(d14.5)'),QVAP_0
  read (unit=10, iostat = ios,fmt = '(a79)'),dummytxt
  read (unit=10, iostat = ios,fmt = '(i6)'),CONST_CVAP

close(10)





print *, 'IMAX=         ',IMAX
print *, 'TEMP=         ',TEMP
print *, 'RH=           ',RH
print *, 'PRES=         ',PRES
print *, 'NUCSIZE=      ',NUCSIZE
print *, 'NUCRATE=      ',NUCRATE
print *, 'PULSE_LENGTH= ',PULSE_LENGTH
print *, 'COND_ON=      ',COND_ON
print *, 'EVAP_ON=      ',EVAP_ON
print *, 'SINK_ON=      ',SINK_ON
print *, 'COAG_ON=      ',COAG_ON
print *, 'NUC_MECH=     ',NUC_MECH
print *, 'NUC_COEFF=    ',NUC_COEFF
print *, 'NUC_EXP=      ',NUC_EXP
print *, 'NUC_COEFF_ORG=',NUC_COEFF_ORG
print *, 'NUC_EXP_ORG=  ',NUC_EXP_ORG
print *, 'CVAP_0=       ',CVAP_0
print *, 'CONST_CVAP=   ',CONST_CVAP

!pause

MAXNEQ = IMAX
CURNEQ = MAXNEQ

    write(*,'(a)',advance='no'), 'Allocating usage matrix space...'
    allocate(IN_USE(0:IMAX),STAT=alloc_err)
    write(*,'(a,i3)'), 'done. Error:',alloc_err

! mark the bins in the distribution that are used by the code
    IN_USE(0:(NUCSIZE-1)) = 0
    IN_USE(NUCSIZE:IMAX) = 1

print *,'Globals done.'
end subroutine set_globals



subroutine init_coeff_table
use discrete_helpers
integer                   :: mm, m1, m2, alloc_err
character(15)             :: form


write(*,'(a)',advance='no'), 'Allocating coefficient space...'
allocate(COE(0:IMAX,0:IMAX),ECOE(0:IMAX,0:IMAX),STAT=alloc_err)
write(*,'(a,i3)'), 'done. Error:',alloc_err


write(*,'(a)',advance='no'), 'Allocating diameter matrix space...'
allocate(DRYDIAM(0:IMAX),WETDIAM(0:IMAX),DRYMASS(0:IMAX), & 
         WETMASS(0:IMAX),STAT=alloc_err)
write(*,'(a,i3)'), 'done. Error:',alloc_err

print *, 'Setting dry diameter vector...'
! set global arrays for dry diamater and mass
do mm = 0,IMAX
    DRYDIAM(mm) = (mm*molecvol_h2so4*3./(4.*pi))**(1./3.) ! m^3
    DRYMASS(mm) = mm*molecmass_h2so4 ! kg
end do
print *, 'done.'


print *, 'Setting wet diameters etc...'
! set global arrays for wet diameter and mass 
! only does it from acid content 3->
call binsize
print *, 'done.'

!  this is probably wrong; should use hydrated acid mass and diameter

WETDIAM(1) = 3.1762d-10  ! weighted mean diameter from Noppel distribution 
                         ! from the vehkam�ki model at 273.15 and rh = 100 and 
                                 !    ra = 100

WETMASS(1) = 2.2448d-25 ! from same

print *, 'Volume of H2SO4 molecule: ', molecvol_h2so4
print *, 'Diameter of WET H2SO4 molecule: ', WETDIAM(1)
print *, 'Max DRY mass: ', DRYMASS(IMAX)
print *, 'Max WET mess: ', WETMASS(IMAX)
print *, 'WET masses 1-20: ', WETMASS(1:20)
print *, 'DRY masses 1-20: ', DRYMASS(1:20)
print *, 'Max DRY diameter: ', DRYDIAM(IMAX)
print *, 'Max WET diameter: ', WETDIAM(IMAX)

! save the dry diameters, wet diameters
open(unit=20,file='DISC_DIAMETERS.TXT',status='replace',position='rewind')
! some trickery to get a dynamic format... looks a bit stoopid i guess there's a better way but oh well
form = '(XXXXXD14.6)'
write(form(2:6),"(I5.5)"),IMAX
print *, form
write(20,form),DRYDIAM(1:IMAX)
write(20,form),WETDIAM(1:IMAX)
close(20)
!pause


! set the coagulation and evaporation matrices
do m1 = 1,IMAX
 do m2 = 1,IMAX
    COE(m1,m2) = coag_coeff(m1,m2)
 end do
!print *, 'Coagulation coef... m1 =  ',m1,'COEF(m1, 1) = ', COE(m1, 1)
end do

! save the coagulation coefficient of vapour molecule with particles
open(unit=30,file='COND_RATE.TXT',status='replace',position='rewind')
! some trickery to get a dynamic format... looks a bit stoopid i guess there's a better way but oh well
form = '(XXXXXD14.6)'
write(form(2:6),"(I5.5)"),IMAX
print *, form
write(30,form),COE(1,1:IMAX)
close(30)
!pause



! no evaporation
print *, 'Zeroing evaporation...'
forall (m1 = 0:IMAX, m2 = 0:IMAX) ECOE(m1,m2) = 0.
print *, 'done.'

end subroutine init_coeff_table
