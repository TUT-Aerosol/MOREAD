MODULE real16
  INTEGER,PARAMETER :: REAL_16 = SELECTED_REAL_KIND (p=14,r=30)
END MODULE real16

MODULE global_variables
use real16
implicit none


INTEGER   ::IMAX    ! max. number of molecules in clusters

INTEGER ::MAXNEQ,&! max number of diff. equations
CURNEQ      ! number of equations in use

integer,allocatable,save            :: IN_USE(:) ! cluster classes being active

real(real_16), allocatable,save  ::   COE(:,:),& ! collision coefficient matrix
     ECOE(:,:) ! evaporation coefficient matrix

real(real_16), allocatable, save    ::DRYDIAM(:),  &    ! the diameters of clusters
WETDIAM(:), &
DRYMASS(:), &! the masses of the clusters
WETMASS(:), &
COAGSINK(:)! the coagulation sink value


INTEGER, PARAMETER  :: STEPMAX = 10000 ! max number of steps


real(real_16), parameter :: &
pi = 3.141592653589793238, &
pstand_atm = 1.0, & ! atm
boltz=1.3807d-23, & !J/K
rconst=8.3143, &  !J/(molK)
pstand= 1.01325e5, &  !Pa  1 atm pressure
atomunit=1.661e-27  ! kg atomic mass unit




! ******************************* CONSTANTS *********************************************

 real(real_16), parameter :: &
avogadro = 6.0221d23, &   ! Avogadro number (#/mol)
boltzmann = 1.3807d-23, &! Boltzmann constant (J/K)
gravitation = 9.81,&! gravitational acceleration (m/s^2)
zero = 0.0,  &
one  = 1.0,  &
unity = 1.0
! ******************************** Vapour properties***********************************
  real(real_16), parameter ::&
density_h2so4 = 1830., &  !  density (kg/m^3)
molar_h2so4 = 98.08,&  ! molar mass (g/mol)
c_sat_h2so4 = 0.d0,  &! saturation concentration above flat surface (#/m^3)
surf_ten_h2so4= 55.d-3,&! surface tension (N/m)
diff_vol_h2so4 = 51.96,&! diffusion volume (????)
alpha_h2so4 = 1., &! mass accomodation coefficient
molecmass_h2so4   = molar_h2so4*1.d-3/avogadro, & ! mass of one molecule (kg)
molecvol_h2so4   = molecmass_h2so4/density_h2so4 ! volume of one molecule (kg)

! ***************************  OTHER STUFF **********************************************
  real(real_16), parameter :: &
molar_h2o = 18.016,&! molar mass of water (g/mol)
molar_air = 28.965,&!- " -air
molar_nh3 = 17.03!- " -ammonia

  real(real_16), parameter :: &! mass of one molecule (kg)
molec_h2o = molar_h2o*1.d-3/avogadro, &
molec_nh3 = molar_nh3*1.d-3/avogadro

  real(real_16), parameter :: &
scale = 1.d6! scaling factor (m -> um)

  real(real_16), parameter :: &
BIG_DENS = 1.0d3! density of background particles kg/m^3
 
!**********************************************************************************




! physical variables

real(real_16) :: RH, TEMP, PRES

integer  :: NUCSIZE  !  size of nucleating particles
real(real_16) :: NUCRATE ! rate of particle appearance per m^3 per s
real(real_16) :: PULSE_LENGTH ! LENGTH of time that particles appear
real(real_16) :: CVAP_0 ! H2SO4 concentration
real(real_16) :: QVAP_0 ! H2SO4 source rate


integer  :: COND_ON, EVAP_ON, SINK_ON, COAG_ON, NUC_MECH
real(real_16) ::  NUC_COEFF,NUC_EXP,NUC_COEFF_ORG,NUC_EXP_ORG

end MODULE global_variables


MODULE unixenv
 IMPLICIT NONE
 INTERFACE
  FUNCTION iargc()
   ! argumenttien lukum��r�
   INTEGER :: iargc
  END FUNCTION iargc
 END INTERFACE
 INTERFACE
  SUBROUTINE getarg(k, arg)
   ! palauttaa k:nnen argumentin
   INTEGER, INTENT(IN) :: k
   CHARACTER(LEN=*), INTENT(OUT) :: arg
  END SUBROUTINE getarg
 END INTERFACE
END MODULE unixenv

module discrete_helpers
use global_variables
 contains

! coag_coeff
! binsize
! binvalidity

function coag_coeff(ix1,ix2)
! the coagulation coefficient of particles 
! with the diameters DIAM(ix1) and DIAM(ix2)

! adopted from UHMA by H. Korhonen
implicit none
   integer, intent(in) :: ix1, ix2  ! the indexes into the distribution arrays
   real(real_16) :: coag_coeff
   real(real_16) :: free_path, viscosity, vel_12, rad_12, dif_12
   real(real_16) :: prov1, prov2, dist, continuum, free_molec
   real(real_16), dimension(2) :: radii, m_part, knudsen, corr, diff, veloc, omega

 ! air mean free path (m) and viscosity (kg/ms) (where from?)
  free_path = (6.73d-8*TEMP*(1.+110.4/TEMP))/(296.*PRES/PSTAND*1.373)
  viscosity = (1.832d-5*406.4*TEMP**1.5)/(5093*(TEMP+110.4))

  ! for both sections
  radii(1)  = WETDIAM(ix1)! radii of colliding particles (m)
  radii(2)  = WETDIAM(ix2)
  m_part(1) = WETMASS(ix1)! masses of colliding particles (kg)
  m_part(2) = WETMASS(ix2)

  knudsen = free_path/radii! particle Knudsen number
  corr = 1. + knudsen*(1.142+0.558*exp(-0.999/knudsen))! Cunninghan correction factor (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
  diff = boltzmann*TEMP*corr/(6*pi*viscosity*radii)! particle diffusion coefficient (m^2/s)
  veloc = sqrt((8.*boltzmann*TEMP)/(pi*m_part))! mean thermal velocity of a particle (m/s)
  omega = 8.*diff/(pi*veloc)! ????? (m)

  vel_12 = sqrt(veloc(1)**2 + veloc(2)**2)! mean relative thermal velocity
  rad_12 = sum(radii)
  dif_12 = diff(1)+diff(2)! relative diffusion coefficient
  continuum = 4.*pi*rad_12*dif_12! flux in continuum regime
  free_molec = pi*vel_12*rad_12**2! flux in free molecular regime

 ! flux matching according to Fuchs (1964) (e.g. Seinfeld & Pandis p. 661)
 prov1 = (rad_12+omega(1))**3 - (rad_12**2 + omega(1)**2)**1.5 
 prov1 = prov1/(3.*rad_12*omega(1)) - rad_12
 prov2 = (rad_12+omega(2))**3 - (rad_12**2 + omega(2)**2)**1.5 
 prov2 = prov2/(3.*rad_12*omega(2)) - rad_12
 dist = sqrt(prov1**2 + prov2**2)! distance at which fluxes are matched

 ! coagulation coefficient between particles [m^3/s]
 coag_coeff = continuum / (rad_12/(rad_12+dist) + continuum/free_molec)

end function coag_coeff


function coag_coeff_sink(ix1,diam_big)
! the coagulation coefficient of particles 
! with the diameters DIAM(ix1) and diam_big
! this one is used for the extramodal (coagulation sink)
! calculation

! adopted from UHMA ( H. Korhonen) by Miikka Dal Maso
implicit none
   integer, intent(in) :: ix1  ! the indexes into the distribution arrays
   real(real_16) :: coag_coeff_sink, &
  diam_big
   real(real_16) :: free_path, viscosity, vel_12, rad_12, dif_12
   real(real_16) :: prov1, prov2, dist, continuum, free_molec
   real(real_16), dimension(2) :: radii, m_part, knudsen, corr, diff, veloc, omega

 ! air mean free path (m) and viscosity (kg/ms) (where from?)
  free_path = (6.73d-8*TEMP*(1.+110.4/TEMP))/(296.*PRES/PSTAND*1.373)
  viscosity = (1.832d-5*406.4*TEMP**1.5)/(5093*(TEMP+110.4))

  ! for both sections
  radii(1)  = WETDIAM(ix1)! radii of colliding particles (m)
  radii(2)  = diam_big
  m_part(1) = WETMASS(ix1)! masses of colliding particles (kg)
  m_part(2) = 4./3.*pi*(diam_big**3.)*BIG_DENS 

  knudsen = free_path/radii! particle Knudsen number
  corr = 1. + knudsen*(1.142+0.558*exp(-0.999/knudsen))! Cunninghan correction factor (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
  diff = boltzmann*TEMP*corr/(6*pi*viscosity*radii)! particle diffusion coefficient (m^2/s)
  veloc = sqrt((8.*boltzmann*TEMP)/(pi*m_part))! mean thermal velocity of a particle (m/s)
  omega = 8.*diff/(pi*veloc)! ????? (m)

  vel_12 = sqrt(veloc(1)**2 + veloc(2)**2)! mean relative thermal velocity
  rad_12 = sum(radii)
  dif_12 = diff(1)+diff(2)! relative diffusion coefficient
  continuum = 4.*pi*rad_12*dif_12! flux in continuum regime
  free_molec = pi*vel_12*rad_12**2! flux in free molecular regime

 ! flux matching according to Fuchs (1964) (e.g. Seinfeld & Pandis p. 661)
 prov1 = (rad_12+omega(1))**3 - (rad_12**2 + omega(1)**2)**1.5 
 prov1 = prov1/(3.*rad_12*omega(1)) - rad_12
 prov2 = (rad_12+omega(2))**3 - (rad_12**2 + omega(2)**2)**1.5 
 prov2 = prov2/(3.*rad_12*omega(2)) - rad_12
 dist = sqrt(prov1**2 + prov2**2)! distance at which fluxes are matched

 ! coagulation coefficient between particles [m^3/s]
 coag_coeff_sink = continuum / (rad_12/(rad_12+dist) + continuum/free_molec)

end function coag_coeff_sink

subroutine binsize()
! routine to get the wet size and mass of sulphuric acid particles
! adopted from UHMA by H. Korhonen (et al.)
implicit none

integer :: ii
real(real_16) :: vol_sa, vol_tot, nacid, nwater, wet_rad, wet_mass, rd, z, &
vol_org, norg, r_dorg, r_worg

call binvalidity()! temperature and rh limits of the parameterizations
! this is broken atm;
do ii = 3,IMAX

nacid = real(ii)! number of acid molecules -"-

if (nacid < 3.or.nacid > 2e11) stop &! validity of parameterizations
'number of sulphuric acid molecules in particles out of bounds'

! parameterization for particle 'wet' radius (accuracy: < +/- 10%)
rd = DRYDIAM(ii)*1d9/2.! radius of sulphuric acid content (nm)

rd = log(rd)
z = 0.3571062410312164 + log(RH)*0.1005557226829263 + rd*(1.072418249000919 - rd*0.007225150816512805)
wet_rad = exp(z)*1.d-9! radius of particle containing sulphuric acid - water content (m)

! parameterization for number of water molecules in the particle
! (accuracy: - 15%-+22%; for rh 0.20-0.90: -8%-+18%)
nwater = exp( 0.7175349331751585 + 3.516681581495114*RH - 7.799949839852581*RH**2 &
+ 4.399706826114728*RH**3 - 0.003003960532620574*TEMP &
+ 0.003600567492582062*RH*TEMP + 1.168794974365219*log(nacid) &
               + 0.07550157096542515*RH*log(nacid) &
+ 0.1230059777461879*RH**2*Log(nacid) - 0.00004695130053660258*TEMP*Log(nacid) &
- 0.01003371715122718*Log(nacid)**2 - 0.005132252262920538*RH*Log(nacid)**2 &
+ 0.0002242756945069555*Log(nacid)**3)

! mass of sulphuric acid-water content and the whole particle (kg)
wet_mass = nacid*molecmass_h2so4 + nwater*molec_h2o


! updating global arrays 
WETMASS(ii) = wet_mass  ! (kg)
WETDIAM(ii) = wet_rad*2.! m


end do
 
  end subroutine binsize
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  subroutine binvalidity()

if (temp < 190.) then
print *,'Warning: temperature < 190 K, using 190 K'
temp = 190.
stop
else if(temp > 330.) then
print *,'Warning: temperature > 330 K, using 330 K'
temp = 330.
stop
end if
if (rh < 0.15) then
print *,'Warning: rh < 15%, using 15%'
rh = 0.15
stop
else if(rh > 0.94) then
print *,'Warning: rh > 94%, using 94%'
rh = 0.94
stop
end if

  end subroutine binvalidity


subroutine make_coag_sink

implicit none
real(real_16), allocatable,save  ::   sinkdistN(:), &
  sinkdistD(:)! pre-existing distribution

integer ::   sinksize, &! size of pre-existing distribution
  i1,i2, ios, alloc_error
character(90) :: distfile

distfile = 'SINKDIST.TXT'

open(unit=10, file = distfile ,status = 'old', action='read', position = 'rewind',iostat=ios)
read (unit=10, iostat = ios,fmt = '(i3)'),sinksize

allocate(sinkdistN(sinksize), sinkdistD(sinksize), &
 STAT = alloc_error)

print *, 'Sink size: ', sinksize


read(10,*),(sinkdistD(i1),i1=1,sinksize)
print *, 'In dist: ', sinkdistD
read(10,*),(sinkdistN(i1),i1=1,sinksize)
print *, 'In dist N: ', sinkdistN
close(10)




do i1 = 1,IMAX
do i2 = 1,sinksize
COAGSINK(i1) = COAGSINK(i1) + coag_coeff_sink(i1,sinkdistD(i2))*sinkdistN(i2)*1.0e6
end do
!print *,'i = ', i1, ' sink = ', COAGSINK(i1)
end do


end subroutine make_coag_sink



end module discrete_helpers
