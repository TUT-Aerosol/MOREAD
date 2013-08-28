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

do mm=0,IMAX
     N(mm) = max(0.0,y(mm))
    ydot(mm) = 0.0
    N(0)=0.0
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

        ydot(mm) = ydot(mm)-coe(1,mm)*max(N(mm),0.0)*N(1)
        ydot(mm+1) = ydot(mm+1)+coe(1,mm)*N(mm)*N(1)
        
        ! make condensing vapour variable
        ydot(1) = ydot(1)-coe(1,mm)*max(N(mm),0.0)*N(1)
        
    end do
end if

! external coagulation and external condensation
if (SINK_ON .eq. 1) then
    print *,'Sink...'
    do mm = NUCSIZE,IMAX
        ydot(mm) = ydot(mm)-COAGSINK(mm)*N(mm)
        ! condensing vapour is lost 
        ydot(1) = ydot(1) - COAGSINK(1)*N(1)
    end do
end if

! nucleation
if (t .le. PULSE_LENGTH) then
!    ydot(NUCSIZE) = ydot(NUCSIZE)+NUCRATE
!    print *, 'Nucleating... nucrate/ydot(nucsize)', NUCRATE, ydot(NUCSIZE)

!   1. Activation nucleation
   ydot(NUCSIZE) = ydot(NUCSIZE)+2.0e-7*N(1)

!    2. Kinetic nucleation
!  ydot(NUCSIZE) = ydot(NUCSIZE)+(1.0e-20*N(1)**2)
 !  ydot(NUCSIZE) = ydot(NUCSIZE)+(1.0e-21*N(1)**2)
!  3. free exponent nucleation
!   ydot(NUCSIZE) = ydot(NUCSIZE)+(Xe-14*N(1)**nuc_exp)*1e6
!  4. free exponent nucleation, acid and organic
!   ydot(NUCSIZE) = ydot(NUCSIZE)+(Xe-14*N(1)**nuc_exp)*((Xe-14*N(2)**nuc_exp_org))*1e6


end if


! ydot(0) = big_flux

if (t .ge. 1.0e20) then
 !print *,'Fend: ydot', ydot
!pause
else
  !  print *,'Fend: ydot', ydot
   ! print *,'Fend: N(0)', N(0)
    !print *,'Fend: N', N
end if
end subroutine f_kine_vode_f90


END MODULE DISCRETE_DVODE_F90
