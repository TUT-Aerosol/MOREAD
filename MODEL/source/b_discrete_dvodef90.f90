MODULE DISCRETE_DVODE_F90
    CONTAINS

!! *** FUNCTION DESCRIBING THE DIFF EQUATIONS FOR DVODE
subroutine f_kine_vode_f90(neq_f, t, y, ydot)
     use real16
     use global_variables

    implicit none

    real(real_16)            ::     big_clusters,    &    ! number of clusters grown out of the system
                                    big_flux                 !  rate of clusters moving out of the system (per dt)
  integer :: neq_f
  double precision :: t
  double precision, dimension(CURNEQ) :: y
  double precision, dimension(CURNEQ) :: ydot

  INTENT(IN)  :: neq_f, t, y
  INTENT(OUT) :: ydot

  real(real_16)            :: N(0:IMAX)  ! to eliminate subzero elements

  integer                     :: mm, m1, m2

 print *,'Fstart, t= ',t
print *,'COAGS 1 = ', COAGSINK(1)

do mm=0,IMAX
     N(mm) = max(0.0,y(mm))
    ydot(mm) = 0.0
    N(2)=0.0 ! growth out of range
    
end do

big_flux = 0.0




if (COAG_ON .EQ. 1) THEN
    print *,'Coag...'
    do m1=NUCSIZE,IMAX        ! calculate source and loss
        do m2=NUCSIZE,IMAX
            if ((m1 .gt. 0.0) .and. (m2 .gt. 0.0)) then ! eliminating nonzero elements

                ! internal coagulation loss of class m1
                ydot(m1)=ydot(m1)-0.5*coe(m1,m2)*N(m1)*N(m2)
                ydot(m2)=ydot(m2)-0.5*coe(m1,m2)*N(m1)*N(m2)

                if (((m1+m2) .LE. IMAX)) THEN
                    ydot(m1+m2)=ydot(m1+m2) + 0.5*coe(m1,m2)*N(m1)*N(m2)
                else
                    big_flux = big_flux + 0.5*coe(m1,m2)*N(m1)*N(m2)
                end if

            end if ! nonzero

        end do         ! do m1
    end do            ! do m2
end if


if (EVAP_ON .EQ. 1) THEN
    print *,' Evap...'
    do m1=NUCSIZE,IMAX        ! calculate source and loss
        do m2=NUCSIZE,IMAX
                ! evaporation fluxes
                if ((m1+m2) .LE. IMAX)  THEN
                    ydot(m1+m2)=ydot(m1+m2) - 0.5*N(m1+m2) * ecoe(m1,m2)
                    ydot(m1) = ydot(m1)+0.5*N(m1+m2) * ecoe(m1,m2)
                    ydot(m2) = ydot(m2)+0.5*N(m1+m2) * ecoe(m1,m2)
                end if
        end do         ! do m1
    end do            ! do m2
end if

! condensation
if (COND_ON .eq. 1) then
    print *,'Cond...'
    ! monomer collisions
    do mm = NUCSIZE,IMAX

        ydot(mm) = ydot(mm)-coe(1,mm)*max(N(mm),0.0)*max(N(1),0.0)
        ydot(mm+1) = ydot(mm+1)+coe(1,mm)*N(mm)*max(N(1),0.0)
        
        ! make condensing vapour variable
        if (CONST_CVAP .eq. 1) then 
       		ydot(1) = 0.0
        else
        	ydot(1) = ydot(1)-coe(1,mm)*max(N(mm),0.0)*max(N(1),0.0)        
        end if 
    end do
    
    ! get the number of particles growing out of the matrix
    big_flux = big_flux + ydot(IMAX)+coe(1,IMAX)*N(IMAX)*max(N(1),0.0)
    
    
end if

! external coagulation and external condensation
if (SINK_ON .eq. 1) then
    print *,'Sink...'
    do mm = NUCSIZE,IMAX
        ydot(mm) = ydot(mm)-COAGSINK(mm)*N(mm)
    end do
       ! condensing vapour is lost
        if (CONST_CVAP .eq. 1) then 
        	ydot(1) = 0.0
        else
        	ydot(1) = ydot(1) - COAGSINK(1)*max(N(1),0.0)
        end if 
        print *,'dCVAP = ',ydot(1)
 
end if

! nucleation
if (t .le. PULSE_LENGTH) then
!    ydot(NUCSIZE) = ydot(NUCSIZE)+NUCRATE
!    print *, 'Nucleating... nucrate/ydot(nucsize)', NUCRATE, ydot(NUCSIZE)


select case (NUC_MECH)

!   1. Activation nucleation
	case (1)
   		ydot(NUCSIZE) = ydot(NUCSIZE)+NUC_COEFF*N(1)
   	case (2) ! kinetic nucleation
        ydot(NUCSIZE) = ydot(NUCSIZE)+(NUC_COEFF*N(1)**2)
    case (3) ! constant J
        ydot(NUCSIZE) = ydot(NUCSIZE)+NUCRATE
    case (4) !  3. free exponent nucleation
    	ydot(NUCSIZE) = ydot(NUCSIZE)+(NUC_COEFF*N(1)**nuc_exp)
	case default
		ydot(NUCSIZE) = 0.0
end select

if (CONST_CVAP .eq. 1) then 
   	ydot(1) = 0.0
else
   	ydot(1) = ydot(1)+QVAP_0; 
end if 

end if

if (t .ge. 1.0e20) then
 !print *,'Fend: ydot', ydot
!pause
else
  !  print *,'Fend: ydot', ydot
   ! print *,'Fend: N(0)', N(0)
    !print *,'Fend: N', N
end if

ydot(2) = big_flux;

! make absolutely sure that if we want a constant CVAPOR then it really stays constant
if (CONST_CVAP .eq. 1) then
	ydot(1) = 0.0
end if

end subroutine f_kine_vode_f90


END MODULE DISCRETE_DVODE_F90
