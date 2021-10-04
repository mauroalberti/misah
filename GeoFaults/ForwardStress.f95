module kinds_mod

	integer, parameter :: i1b = selected_int_kind(1)                    
	integer, parameter :: i2b = selected_int_kind(4)                    
	integer, parameter :: i4b = selected_int_kind(9)  
	integer, parameter :: r4b = selected_real_kind( 6, 37)              
	integer, parameter :: r8b = selected_real_kind(15,307) 

end module kinds_mod

!------------------------------------------------ 


 
!------------------------------------------------  


module data_types

use kinds_mod

! variable declarations

implicit none

type :: range_1D
	real (kind=r8b) :: min,max
end type range_1D

type :: range_3D
	type (range_1D) :: r1D(3)
end type range_3D

type :: vector  
	sequence             
	real (kind=r8b)	:: x,y,z    
end type vector

type :: axis
	real (kind=r8b)	:: trend, plunge
end type axis

type :: orthonormal_triad
	type(vector) :: X,Y,Z
end type orthonormal_triad

type :: plane
	real (kind=r8b)	:: strike_rhr,dipdirection,dipangle
end type plane
    
type :: fault_datum
	integer (kind=i4b) :: id
	type (plane) :: fltplane
	type (axis) :: slickenline    
	real (kind=r8b)	:: rake_aki
	type(vector) :: faultnorm_vect,slickenl_vect
    real(kind=r8b) :: x,y,z,t
end type fault_datum

type :: stress_princcompon
	type (axis) :: S1_pole, S2_pole, S3_pole
	type (vector) :: S1_vector, S2_vector, S3_vector
    real (kind=r8b) :: sigma1, sigma2, sigma3, phi
end type stress_princcompon

type :: stress_solution
	logical :: valid_solution
	type(vector) :: traction_vect, normalstress_vect, shearstress_vect
	real (kind=r8b)	:: theor_rake,tractionvect_magn,normalstress_magn,shearstress_magn
    real (kind=r8b)	:: modif_sliptendency, deformation_index 
    type (axis) :: theor_slickenline
end type stress_solution

type :: stressanalysis_param
    type(stress_princcompon) :: stressprinccomp1
	real (kind=r8b) :: tens(3,3)
	real (kind=r8b) :: rot_matrix(3,3)
end type stressanalysis_param

type :: record_variables
	logical :: spatial_var,time_var 
    integer (kind=i1b) :: spat_coord_system, readcase    
end type record_variables

type :: faultsimulation_param    
	integer (kind=i4b) :: faults_totnumb
    type (range_3D) :: spat_bound
	type (range_1D) :: temp_bound
end type faultsimulation_param

type :: program_parameters
	integer (kind=i4b) :: time_initial(8)
	type (record_variables) :: fault_loctimevar 
	integer (kind=i1b) :: input_case 
 	character (len=50) :: inputfile_name 
    integer (kind=i2b) :: inputfile_headerrownumb   
	type (faultsimulation_param) :: faultsimulpar
	integer (kind=i4b) :: faults_totnumb
end type program_parameters


integer :: ios !status of input/output connection
real (kind=r8b) :: pi  	! pi radians
real (kind=r8b) :: r2d, d2r	! for conversion from radians to degrees

type(orthonormal_triad), parameter :: frame0 = orthonormal_triad(vector(1.0,0.0,0.0),vector(0.0,1.0,0.0),vector(0.0,0.0,1.0))
real (kind=r8b),parameter :: unitary_vect_tolerance = 1.0e-6 !maximum accepted value for unitary vector magnitude difference with unit
real (kind=r8b),parameter :: vect_normaliz_tolerance = 1.0e-5 !minimum accepted value for vector magnitude to apply normalization
real (kind=r8b), parameter :: shearmagnitude_minthresh = 1.0e-5 !minimum accepted value for shear stress to be considered meaningful

end module data_types


!------------------------------------------------  


!------------------------------------------------    

module vector_processing

use data_types


contains

!------------------------
! calculates the sum of two vectors

type(vector) function vector_sum(vector1,vector2)

type(vector), intent(in) :: vector1, vector2 

vector_sum%x = vector1%x+vector2%x
vector_sum%y = vector1%y+vector2%y
vector_sum%z = vector1%z+vector2%z

end function vector_sum

!------------------------


!------------------------
! calculates the difference between two vectors

type(vector) function vector_diff(vector1,vector2)

type(vector), intent(in) :: vector1, vector2 

vector_diff%x = vector1%x-vector2%x
vector_diff%y = vector1%y-vector2%y
vector_diff%z = vector1%z-vector2%z

end function vector_diff

!------------------------


!------------------------
! angle (in radians) between two vectors (radiants, 0-pi)

real (kind=r8b) function vector_angle_rad(vector1,vector2)

type(vector), intent(in) :: vector1, vector2 
    
	vector1_magn = vector_magn(vector1)
	vector2_magn = vector_magn(vector2)
	if ((vector1_magn < vect_normaliz_tolerance).or.(vector2_magn < vect_normaliz_tolerance)) then
		write(*,*) 'Error in vector magnitude (function vector_angle_rad). Hit any key to stop'
		read(*,*)
    	stop	! STOP PROGRAM FOR MAJOR ERROR IN DATA INPUT
	end if
	! scalar product between two vectors
	scaledscalarproduct = vector_scalprod(vector1, vector2)/(vector1_magn*vector2_magn)  
     
	! angle between vectors (in radians)
    
	if (scaledscalarproduct < -1.) then 
		vector_angle_rad = pi
    else if (scaledscalarproduct > 1.) then 
		vector_angle_rad = 0.0
    else
    	vector_angle_rad = acos(scaledscalarproduct)
    end if

end  function vector_angle_rad 
  
      
!------------------------


!------------------------
! angle (in radians) between two axes (radiants, 0-pi/2)

real (kind=r8b) function axes_angle_rad(vector1,vector2)

type(vector), intent(in) :: vector1, vector2 

	! angle between vectors (in radians)
	axes_angle_rad = vector_angle_rad(vector1,vector2)
    axes_angle_rad = min(axes_angle_rad, pi-axes_angle_rad)

end function axes_angle_rad


!------------------------


!------------------------
! vector normalization 

type(vector) function vector_normalization(vector1) result(vector2)

type(vector), intent(in) :: vector1 

	vector1_magn = vector_magn(vector1)

    if (vector1_magn < vect_normaliz_tolerance) then
      write(*,*) 'Error in vector magnitude processing. Hit any key to stop'
      read (*,*)
      stop ! STOP PROGRAM FOR MAJOR ERROR IN DATA INPUT
    end if

	vector2%x = vector1%x/vector1_magn
	vector2%y = vector1%y/vector1_magn
    vector2%z = vector1%z/vector1_magn  
    

end function vector_normalization

!------------------------


!------------------------
! vector magnitude 

real (kind=r8b) function vector_magn(vector1)

type(vector), intent(in) :: vector1 

   	vector_magn = sqrt((vector1%x)**2+(vector1%y)**2+(vector1%z)**2)

end function vector_magn

!------------------------


!------------------------
! calculates the product of a vector by a scalar
type (vector) function vectorbyscalar(vector1, scalar1)

type(vector), intent(in) :: vector1
real (kind=r8b), intent(in) :: scalar1

vectorbyscalar%x = scalar1*vector1%x
vectorbyscalar%y = scalar1*vector1%y
vectorbyscalar%z = scalar1*vector1%z

end function vectorbyscalar

!------------------------
! scalar product of two vectors (given as their cartesian coordinates)
      
real (kind=r8b) function vector_scalprod(vector1, vector2)

type(vector), intent(in) :: vector1, vector2 

	vector_scalprod = vector1%x*vector2%x+vector1%y*vector2%y+vector1%z*vector2%z  

end function vector_scalprod

!------------------------


!------------------------
! vectorial product of two vectors (given as their cartesian coordinates) 
      
type(vector) function vector_vectprod(vector1,vector2)

type(vector), intent(in) :: vector1, vector2 

	vector_vectprod%x=(vector1%y*vector2%z)-(vector1%z*vector2%y)
	vector_vectprod%y=(vector1%z*vector2%x)-(vector1%x*vector2%z)
	vector_vectprod%z=(vector1%x*vector2%y)-(vector1%y*vector2%x)

end function vector_vectprod


!------------------------


!------------------------
! vetctor1 projection on vector2
type (vector) function vector_projection(vector1,vector2)

type(vector), intent(in) :: vector1, vector2 
real (kind=r8b) :: scalprod

scalprod = vector_scalprod(vector1, vector_normalization(vector2))
vector_projection = vectorbyscalar(vector2, scalprod)

end function vector_projection

!------------------------


!------------------------
! test if a vector has magnitude = 1
logical function vect_normaliztest(vector1) result(vect_is_normalized)

!  QUATv_normalizedvectortest
!  created 2005-02-15

type(vector), intent(in) :: vector1
real (kind=r8b) :: vector1_magn


vector1_magn = vector_magn(vector1)

abs_diff = dabs(1-vector1_magn)

if (abs_diff > unitary_vect_tolerance) then
    vect_is_normalized = .false.
else
    vect_is_normalized = .true.
endif

end function vect_normaliztest


!------------------------


!------------------------
! converts a 3D vector to a 3x1 array

function vector2array(vector1) result(array1)

type(vector), intent(in) :: vector1
real(kind=r8b) :: array1(3)

array1 = (/vector1%x, vector1%y, vector1%z/)

end function vector2array

!------------------------



!------------------------
! converts a 3x1 array to a 3D vector 

function array2vector(array1) result(vector1)

real(kind=r8b), intent(in) :: array1(3)
type(vector):: vector1

vector1%x = array1(1)
vector1%y = array1(2)
vector1%z = array1(3)

end function array2vector

!------------------------



end module vector_processing


!------------------------------------------------


!------------------------------------------------  

module geometric_processing

use vector_processing

contains


!------------------------
! calculates the polar components from the cartesian ones
type(vector) function pole2cartesian(axis1) result(vector1)

!  QUATg_pole2cartesian
!  created 2002-12-30

type(axis), intent(in) :: axis1


vector1%x = cos(d2r*axis1%plunge)*cos(d2r*axis1%trend)
vector1%y  = cos(d2r*axis1%plunge)*sin(d2r*axis1%trend)
vector1%z  =  sin(d2r*axis1%plunge)


end function pole2cartesian

!------------------------

        
!------------------------
! calculates polar components from cartesian ones

type(axis) function cartesian2pole(vector1) result(axis1)

!  QUATg_cartesian2pole
!  created 2003-01-01
!  modified: 2005-02-15

type(vector):: vector1 
logical :: vect_is_normalized

vect_is_normalized = vect_normaliztest(vector1)

if (.not.vect_is_normalized) then                            
  vector1 = vector_normalization(vector1)
endif

! polar coordinates calculation

if (vector1%z > 1.0_r8b) then
	vector1%z = 1.0_r8b
elseif (vector1%z < -1.0_r8b) then	
	vector1%z = -1.0_r8b
endif

axis1%plunge = r2d*dasin(vector1%z)


if (dabs(axis1%plunge)> 89.5) then
		axis1%trend = 0.0
else
	axis1%trend = r2d*atan2(vector1%y,vector1%x)
	if (axis1%trend < 0.0) then
		axis1%trend = 360.0 + axis1%trend
    endif
endif

end function cartesian2pole

!------------------------


!------------------------
! calculates the down axis from the axis
type(axis) function axis2downaxis(axis1) result(axis2)

!  QUATg_axis2downaxis
!  created 2003-02-01

type(axis), intent(in) :: axis1

if (axis1%plunge < 0.0) then
  axis2%plunge = -axis1%plunge
  axis2%trend = axis1%trend + 180.0
  if (axis2%trend >= 360.0) then
    axis2%trend = axis2%trend - 360.0
  endif 
else
  axis2%trend = axis1%trend
  axis2%plunge = axis1%plunge
endif 

end function axis2downaxis

  
!--------------------------



end module geometric_processing

!----------------------------------------------



!----------------------------------------------

module fault_processing

use geometric_processing

implicit none

contains

!--------------------------
! calculates dip direction of a fault 
subroutine dipdir_calc(faultplane1) 

type(plane), intent(inout) :: faultplane1
faultplane1%dipdirection = faultplane1%strike_rhr + 90.0
if (faultplane1%dipdirection >= 360.0) then
  faultplane1%dipdirection = faultplane1%dipdirection - 360.0
endif

end subroutine dipdir_calc

!--------------------------


!--------------------------
! calculates fault normal

type(vector) function faultplanenormal_calc(faultplane1) result(faultnorm)

type(plane), intent(in) :: faultplane1

! Fault Normal cartesian coordinates
! formulas from Aki and Richards, 1980
faultnorm%x = -sin(d2r*faultplane1%dipangle)*sin(d2r*faultplane1%strike_rhr)
faultnorm%y = sin(d2r*faultplane1%dipangle)*cos(d2r*faultplane1%strike_rhr)
faultnorm%z = -cos(d2r*faultplane1%dipangle)

end function faultplanenormal_calc

!--------------------------



!--------------------------
! calculates vector (cartesian) components of fault record
subroutine fault_cartcmp(faultplane1, rake1, faultnorm_vect, slickenl_vect) 

! QUATf_faultpole2faultvector

type(plane), intent(in) :: faultplane1
real (kind=r8b)	:: rake1 
type (vector) , intent(out) :: faultnorm_vect, slickenl_vect

! Fault Normal and Slickenline cartesian coordinates
! formulas from Aki and Richards, 1980
faultnorm_vect%x = -sin(d2r*faultplane1%dipangle)*sin(d2r*faultplane1%strike_rhr)
faultnorm_vect%y = sin(d2r*faultplane1%dipangle)*cos(d2r*faultplane1%strike_rhr)
faultnorm_vect%z = -cos(d2r*faultplane1%dipangle)

slickenl_vect%x = cos(d2r*rake1)*cos(d2r*faultplane1%strike_rhr)   &
 +sin(d2r*rake1)*cos(d2r*faultplane1%dipangle)*sin(d2r*faultplane1%strike_rhr)
slickenl_vect%y = cos(d2r*rake1)*sin(d2r*faultplane1%strike_rhr)  &
 -sin(d2r*rake1)*cos(d2r*faultplane1%dipangle)*cos(d2r*faultplane1%strike_rhr)
slickenl_vect%z = -sin(d2r*rake1)*sin(d2r*faultplane1%dipangle)


end subroutine fault_cartcmp

!--------------------------


!--------------------------
! calculate the slickenline (trend and plunge) from the rake angle
type(axis) function rake2slickenline(strike, dip, rake) result(slicken_pole) 

!  Structural_Med_CalculateSlickenline
!  created: 2005-09-11
!  modified: 2008-01-16

implicit none

real (kind=r8b), intent(in) :: strike, dip, rake
real (kind=r8b) :: strike_rd, dip_rd, rake_rd
type (vector) :: slick_vect

strike_rd = d2r*strike
dip_rd = d2r*dip
rake_rd = d2r*rake

slick_vect%x = cos(rake_rd)*cos(strike_rd)+sin(rake_rd)*sin(strike_rd)*cos(dip_rd)
slick_vect%y = cos(rake_rd)*sin(strike_rd)-sin(rake_rd)*cos(strike_rd)*cos(dip_rd) 
slick_vect%z = -sin(rake_rd)*sin(dip_rd)
 
! determination of trend and plunge of Slickenline
slicken_pole = cartesian2pole(slick_vect)


end function rake2slickenline
!--------------------------



end module fault_processing


!----------------------------------------------



!----------------------------------------------


module stress_processing

use fault_processing


implicit none


contains


!--------------------------
! calculation of S2_pole as vector product of S3_pole and S1_pole
subroutine S2_calc(stressprinccomp1)

type(stress_princcompon), intent(inout) :: stressprinccomp1

stressprinccomp1%S2_vector = vector_vectprod(stressprinccomp1%S3_vector,stressprinccomp1%S1_vector)
stressprinccomp1%S2_pole = cartesian2pole(stressprinccomp1%S2_vector)
stressprinccomp1%S2_pole = axis2downaxis(stressprinccomp1%S2_pole)

end subroutine S2_calc

!--------------------------



!--------------------------
! ! calculation of sigma2 from S1_pole, S3_pole and PHI values
subroutine sigma2_calc(stressprinccomp1)

type(stress_princcompon), intent(inout) :: stressprinccomp1

stressprinccomp1%sigma2 = stressprinccomp1%phi*stressprinccomp1%sigma1 				&
							+ (1-stressprinccomp1%phi)*stressprinccomp1%sigma3 


end subroutine sigma2_calc

!--------------------------



!--------------------------
! define rotation matrix based on stress principal axes
! based on Kuipers, 2002, p.161, eqs. 7.8
function rotmatr(stressprinccomp1) result(rot_matrix)

type(stress_princcompon), intent(in) :: stressprinccomp1
real (kind=r8b) :: rot_matrix(3,3)

rot_matrix(1,1) = vector_scalprod(stressprinccomp1%S1_vector, frame0%X)
rot_matrix(1,2) = vector_scalprod(stressprinccomp1%S2_vector, frame0%X) 
rot_matrix(1,3) = vector_scalprod(stressprinccomp1%S3_vector, frame0%X) 

rot_matrix(2,1) = vector_scalprod(stressprinccomp1%S1_vector, frame0%Y)
rot_matrix(2,2) = vector_scalprod(stressprinccomp1%S2_vector, frame0%Y) 
rot_matrix(2,3) = vector_scalprod(stressprinccomp1%S3_vector, frame0%Y) 

rot_matrix(3,1) = vector_scalprod(stressprinccomp1%S1_vector, frame0%Z)
rot_matrix(3,2) = vector_scalprod(stressprinccomp1%S2_vector, frame0%Z) 
rot_matrix(3,3) = vector_scalprod(stressprinccomp1%S3_vector, frame0%Z)
                    

end function rotmatr

!--------------------------


!--------------------------
! calculate stress tensor expressed in frame components
! from Kagan and Knopoff, 1985a, p. 433
function stresstensorcalc(stressprinccomp1,rot_matrix) result(tens)

type(stress_princcompon), intent(in) :: stressprinccomp1
real (kind=r8b), intent(in)  :: rot_matrix(3,3)

real (kind=r8b) :: stresstens0(3,3), tens(3,3)

! stress eigentensor
stresstens0 = reshape((/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/),shape=(/3,3/),order=(/2,1/))
stresstens0(1,1) = stressprinccomp1%sigma1
stresstens0(2,2) = stressprinccomp1%sigma2
stresstens0(3,3) = stressprinccomp1%sigma3

tens = matmul(rot_matrix,matmul(stresstens0,transpose(rot_matrix)))

end function stresstensorcalc

!--------------------------


!--------------------------
! calculation of stress solution for the given plane and stress tensor
! based on Xu, 2004 (Geoph. J. Int., 157,1316-1330) and references within.
type(stress_solution) function stresssolution_calc(tens1,faultrec1) result(stresssolution1)

! QUATstress_stress&fault2rake
!last modified: 2008/01/05


! INPUT
!    SELF with 2 objects:
!      (0): fault plane orientation: strike rhr, dip direction, dip angle
!      (1): stress field tensor components: x,y,z (sigma 1,2,3)

!  OUTPUT
!    normal and tangential stresses magnitudes
!    slip tendency
!    theoretical rake angle of tangential stress

!  STEPS
!    input of data:
!      fault plane orientation: strike rhr, dip direction, dip angle
!      stress field tensor components: x,y,z (sigma 1,2,3)
!    fault normal calculation: n (1,2,3)
!    traction vector c.: sigma (1,2,3)
!    normal stress c.: sigma n (1,2,3)
!    shear stress c: tau s (1,2,3)
!    slip tendency c.: Ts
!    theoretical rake angle c: lambda
!    output of results: 
!      magnitudes of normal and shear stresses, 
!      slip tendency,
!      theoretical rake angle


!    input of data:
!      fault plane orientation: strike rhr, dip direction, dip angle
!      stress field tensor components: x,y,z (sigma 1,2,3)

implicit none

real (kind=r8b), intent(in)  :: tens1(3,3)
type(fault_datum), intent(in)  :: faultrec1

type(vector) :: shearstress_unitvect
type(axis) :: strike_pole, dipdir_pole
type(vector) :: strike_vect, dipdir_vect, faultnorm_vect
real (kind=r8b) :: faultnormal_array(3), scalprod_shearstress_strike, scalprod_shearstress_dipdir


! calculation of fault normal calculation: n (1,2,3)
faultnorm_vect = faultplanenormal_calc(faultrec1%fltplane)
faultnormal_array = vector2array(faultnorm_vect)

!    traction vector c.: sigma (1,2,3)
stresssolution1%traction_vect = array2vector(- matmul(tens1,faultnormal_array))
stresssolution1%tractionvect_magn = vector_magn(stresssolution1%traction_vect)
if (vector_scalprod(stresssolution1%traction_vect, faultnorm_vect)<0.0_r8b) then
	stresssolution1%tractionvect_magn = -stresssolution1%tractionvect_magn
endif

!    normal stress c.: SigmaN
stresssolution1%normalstress_vect = vector_projection(stresssolution1%traction_vect, faultnorm_vect)
stresssolution1%normalstress_magn = vector_magn(stresssolution1%normalstress_vect)
if (vector_scalprod(stresssolution1%normalstress_vect, faultnorm_vect)<0.0_r8b) then
	stresssolution1%normalstress_magn = -stresssolution1%normalstress_magn
endif

!    shear stress c: tau s (1,2,3)
stresssolution1%shearstress_vect = vector_diff(stresssolution1%traction_vect,stresssolution1%normalstress_vect)
stresssolution1%shearstress_magn = vector_magn(stresssolution1%shearstress_vect)


if (stresssolution1%shearstress_magn <= shearmagnitude_minthresh) then 
  stresssolution1%valid_solution = .false.
else
  stresssolution1%valid_solution = .true.
end if

if (stresssolution1%valid_solution) then

	shearstress_unitvect = vector_normalization(stresssolution1%shearstress_vect)

	! theoretical rake angle c: lambda
	strike_pole = axis(faultrec1%fltplane%strike_rhr,0)
	dipdir_pole = axis(faultrec1%fltplane%dipdirection,faultrec1%fltplane%dipangle)

	strike_vect = vector_normalization(pole2cartesian(strike_pole))
	dipdir_vect = vector_normalization(pole2cartesian(dipdir_pole))

	scalprod_shearstress_strike = vector_scalprod(shearstress_unitvect, strike_vect)
                                             
	scalprod_shearstress_dipdir = vector_scalprod(shearstress_unitvect, dipdir_vect)    
                                      
	if (scalprod_shearstress_strike>1) then
		scalprod_shearstress_strike=1
	elseif (scalprod_shearstress_strike<-1) then
		scalprod_shearstress_strike=-1
	endif

	stresssolution1%theor_rake = r2d*acos(scalprod_shearstress_strike)

	if (scalprod_shearstress_dipdir > 0.0_r8b) then
  		stresssolution1%theor_rake = - stresssolution1%theor_rake
	endif

	stresssolution1%theor_slickenline = rake2slickenline(   &
    	faultrec1%fltplane%strike_rhr,faultrec1%fltplane%dipangle,stresssolution1%theor_rake)
    
	stresssolution1%modif_sliptendency = stresssolution1%shearstress_magn/dabs(stresssolution1%tractionvect_magn)

	stresssolution1%deformation_index = &
	(dabs(stresssolution1%tractionvect_magn) - stresssolution1%shearstress_magn)/stresssolution1%tractionvect_magn

endif 

end function stresssolution_calc

!--------------------------

end module stress_processing



!----------------------------------------------





!----------------------------------------------

module inputoutput_defs

use stress_processing

implicit none

contains


!------------------------
! prints initial message to screen
subroutine initial_message()

    write (*,*)
	write (*,"(10x,A)")'----------------------------------------------------------'
	write (*,"(25x,A)") 'ForwardStress'
	write (*,"(25x,A)") 'vers. 1.0, 2008-07'
	write (*,"(10x,A)")'----------------------------------------------------------'
    write (*,*)
    write (*,"(5x,A)") ' Program for forward simulation of homogeneous stress field'
    write (*,"(5x,A)") ' on faults/focal mechanism'      
    write (*,*)
	write (*,"(10x,A)")'----------------------------------------------------------'
    write (*,*) 
    
   
end subroutine initial_message

!------------------------



!-------------------------
! inputs the stress tensor as principal axes orientations and magnitudes
type(stress_princcompon) function define_princcompstresstensor() result(stressprinccomp1)

	do
	!  input of stress tensor values:
	!  S1 (Magnitude, trend and plunge), S3 (Magnitude, trend and plunge) e PHI (scalare, 0-1)

	! stress field input: S1
		write (*,"(A)",advance='no') 'Enter S1 trend and plunge: '
		read (*,*) stressprinccomp1%S1_pole%trend,stressprinccomp1%S1_pole%plunge     
			       
		stressprinccomp1%S1_vector = pole2cartesian(stressprinccomp1%S1_pole) 
		
	! stress field input: S3
		write (*,"(A)",advance='no') 'Enter S3 trend and plunge: '
		read (*,*) stressprinccomp1%S3_pole%trend,stressprinccomp1%S3_pole%plunge

		stressprinccomp1%S3_vector = pole2cartesian(stressprinccomp1%S3_pole) 

	!  check of orthogonality between S1 and S3: if angle lower than 89.5 degrees then 
    !   the user has to enter again the orientations
		if (.not.(r2d*axes_angle_rad(stressprinccomp1%S1_vector,stressprinccomp1%S3_vector)>=89.5)) then
			write (*,*) 'Non-orthogonality between S1 and S3. Enter again values'
      		cycle
        else
          	exit
		endif 

	end do

	!  calculation of S2 as vector product of S3 and S1
		call S2_calc(stressprinccomp1)

	do
		! input of stress field magnitudes
		write (*,"(A)",advance='no') 'Enter Sigma1 (magnitude): ' 
		read (*,*)stressprinccomp1%sigma1

		write (*,"(A)",advance='no') 'Enter Sigma3 (magnitude): ' 
		read (*,*) stressprinccomp1%sigma3

		if (stressprinccomp1%sigma1<=stressprinccomp1%sigma3) then
        	write (*,*) 'Invalid stress magnitudes. Define again'
            cycle
        else
          	exit
        endif
        
	end do

    do
        ! definition of phi value
		write (*,"(A)",advance='no') 'Enter Phi value (0-1): ' 
		read (*,*) stressprinccomp1%phi
		if ((stressprinccomp1%phi<=1).and.(stressprinccomp1%phi>=0)) then
        	exit
        endif		
	end do     
   
	!  calculation of sigma2 (magnitude) from sigma1, sigma3 and PHI values
	call sigma2_calc(stressprinccomp1)


end function define_princcompstresstensor

!-------------------------


!------------------------
! define fault characteristics in input/output
type(record_variables) function faultparam_define() result(fault_loctimevar1)

! spatial location
write(*,"(A)", ADVANCE="no") 'Location information present? T(rue)/F(alse): '
read (*,"(L1)") fault_loctimevar1%spatial_var
! in positive case, ask if planar or spherical coordinates, and boundaries
if (fault_loctimevar1%spatial_var) then
	do
		write(*,"(A)", ADVANCE="no") 'Choose coordinate system: 1 - planar (x-y-z) &
         & , 2 - spherical (lat-lon-z): '
        read (*,*) fault_loctimevar1%spat_coord_system
		if (fault_loctimevar1%spat_coord_system==1 &
          .or. fault_loctimevar1%spat_coord_system ==2) exit
	end do
endif
  
write(*,*)

! time information
write(*,"(A)", ADVANCE="no") 'Time information present? T(rue)/F(alse): '
read (*,"(L1)") fault_loctimevar1%time_var

! definition of input type 
fault_loctimevar1%readcase = 0 !case without spatial and/or time information
if (fault_loctimevar1%spatial_var .and. (.not.fault_loctimevar1%time_var)) fault_loctimevar1%readcase = 1 
if (fault_loctimevar1%spatial_var .and. fault_loctimevar1%time_var) fault_loctimevar1%readcase = 2
if ((.not.fault_loctimevar1%spatial_var) .and. fault_loctimevar1%time_var) fault_loctimevar1%readcase = 3

write(*,*)

end function faultparam_define

!------------------------


!------------------------
! define if input from file or from internal simulation

integer (kind=i1b) function define_input_case() result(input_case)

	do
		write (*,"(A)",advance='no') 'Read data from file (1) or create using simulation (2)?: ' 
		read (*,*) input_case
        if (input_case==1 .or.input_case==2) exit        
    end do
	
end function define_input_case

!------------------------


!------------------------
! define input file

subroutine inputfiledef(programparams1) 

 implicit none

 type(program_parameters), intent(inout) :: programparams1
 character (len=50) :: inputfilename 
 integer :: ios
 logical :: exists

 do	
	write(*,"(A)", ADVANCE="no") 'Enter the name of input file: '
	read(*,*) inputfilename

    programparams1%inputfile_name = trim(inputfilename)
    
	inquire(file=programparams1%inputfile_name,exist=exists)
	if (.NOT.exists) then
		write(*,*) 'INPUT FILE: not found. Program will stop'
		cycle        
	end if

	! open the input file with sequential access
	open(unit=17,file=programparams1%inputfile_name,status='old',access='sequential'  &
        ,form='formatted',iostat=ios)
	if (ios /= 0) then
		write(*,*) 'INPUT FILE: not opened'
		cycle
	end if

	exit
    
 end do


write(*,"(2x,A)", ADVANCE="no") 'Enter number of header rows in input file (default: 1): '
read(*,*) programparams1%inputfile_headerrownumb
      

end subroutine inputfiledef

!------------------------



!------------------------
! definition of fault simulation parameters
subroutine simulationparamdefine(programparams1) 

type(program_parameters), intent(inout) :: programparams1

! define fault number
do
	write(*,"(A)", ADVANCE="no") 'Enter number of faults to simulate: '
	read(*,*) programparams1%faultsimulpar%faults_totnumb
    if (programparams1%faultsimulpar%faults_totnumb >= 1) exit
end do

write(*,*)

! spatial boundaries definition
if (programparams1%fault_loctimevar%spatial_var) then
	do
		write(*,"(A)", ADVANCE="no") 'Enter spatial ranges &
        & (min and max for x,y,z or Lat, Lon, z; e.g. 0 100 0 100 0 100): '
        read (*,*) programparams1%faultsimulpar%spat_bound%r1D
		if (programparams1%faultsimulpar%spat_bound%r1D(1)%max>=programparams1%faultsimulpar%spat_bound%r1D(1)%min .and. &
        	programparams1%faultsimulpar%spat_bound%r1D(2)%max>=programparams1%faultsimulpar%spat_bound%r1D(2)%min .and. &
            programparams1%faultsimulpar%spat_bound%r1D(3)%max>=programparams1%faultsimulpar%spat_bound%r1D(3)%min) exit
	end do 
endif

write(*,*)

! time interval definition
if (programparams1%fault_loctimevar%time_var) then
	do
		write(*,"(A)", ADVANCE="no") 'Enter time ranges (min and max values, e.g. 0 100): '
        read (*,*) programparams1%faultsimulpar%temp_bound
		if (programparams1%faultsimulpar%temp_bound%max>=programparams1%faultsimulpar%temp_bound%min) exit
	end do  

endif

write(*,*)
  
end subroutine simulationparamdefine

!------------------------


!------------------------
! calculates simulations for faults

subroutine generatesimulatedfaults(programparams1)

	type(program_parameters), intent(inout) :: programparams1
 	integer :: ios    
 	integer (kind=i2b) :: num_randomvar, arraycounter_time
    real (kind=r8b), allocatable :: randnumb_array(:)
	integer (kind=i4b) :: i
    real (kind=r8b) :: x_range, y_range, z_range, t_range
    real (kind=r8b) :: x_0, y_0, z_0, t_0    
    character (len=100) :: faultsim_values
	character, parameter  :: SeparatorString = ","
    
	open(unit=17,file='temp_fault_xx.dat',status='replace',access='sequential'  &
        ,form='formatted',iostat=ios)
        
	if (ios /= 0) then
		write(*,*) 'temp fault file cannot be created/opened. Press any key to stop'
        read(*,*)
		stop
	endif
 
	! set to 0 the header number in temp file (for compatibility with case of external input file)
	programparams1%inputfile_headerrownumb = 0

    num_randomvar = 2  ! strike and dip_angle
    ! spatial locations
    if (programparams1%fault_loctimevar%spatial_var) then
      num_randomvar = num_randomvar + 3 ! added 3 loc comp.
      x_range = programparams1%faultsimulpar%spat_bound%r1D(1)%max - programparams1%faultsimulpar%spat_bound%r1D(1)%min
      y_range = programparams1%faultsimulpar%spat_bound%r1D(2)%max - programparams1%faultsimulpar%spat_bound%r1D(2)%min
      z_range = programparams1%faultsimulpar%spat_bound%r1D(3)%max - programparams1%faultsimulpar%spat_bound%r1D(3)%min      

      x_0 = programparams1%faultsimulpar%spat_bound%r1D(1)%min
      y_0 = programparams1%faultsimulpar%spat_bound%r1D(2)%min
      z_0 = programparams1%faultsimulpar%spat_bound%r1D(3)%min
      
	endif
    ! time 
    if (programparams1%fault_loctimevar%time_var) then
       num_randomvar = num_randomvar + 1 ! added 1 time comp.
       t_range = programparams1%faultsimulpar%temp_bound%max - programparams1%faultsimulpar%temp_bound%min
       t_0 = programparams1%faultsimulpar%temp_bound%min 
	endif
     
	! allocate array containing random values
	allocate(randnumb_array(num_randomvar))
  
    do i=1,programparams1%faultsimulpar%faults_totnumb

        call random_seed()
    	call random_number(randnumb_array)
        
        
      	! fault id
		write(faultsim_values,"(i10,A)") i,SeparatorString                


        if (programparams1%fault_loctimevar%spatial_var) then
			! X/Lat
	    	write(faultsim_values,"(A,f12.5,A)") trim(faultsim_values),x_0+(randnumb_array(3)*x_range),SeparatorString
			! Y/Lon
	    	write(faultsim_values,"(A,f12.5,A)") trim(faultsim_values),y_0+(randnumb_array(4)*y_range),SeparatorString
			! Z
	    	write(faultsim_values,"(A,f12.5,A)") trim(faultsim_values),z_0+(randnumb_array(5)*z_range),SeparatorString
        endif
        
    	if (programparams1%fault_loctimevar%time_var) then
        	arraycounter_time = 6
            if (.not.(programparams1%fault_loctimevar%spatial_var)) arraycounter_time = 3
			! time
	    	write(faultsim_values,"(A,f12.5,A)") &
               trim(faultsim_values),t_0+(randnumb_array(arraycounter_time)*t_range),SeparatorString
		endif 

		! strike
	    write(faultsim_values,"(A,f6.1,A)") trim(faultsim_values),randnumb_array(1)*360.0,SeparatorString        

		! dip angle 
		write (faultsim_values,"(A,f6.1,A)") trim(faultsim_values),acos(randnumb_array(2))*180.0/pi

        ! write the entire string to the temporary file
 		write(17, "(A)") trim(faultsim_values) 

  	end do
    
	! deallocate the array storing the random numbers,
    ! since it is no longer needed 
	deallocate(randnumb_array)
    
	rewind(17)


end subroutine generatesimulatedfaults

!-------------------------
! defines the output files 
! containing the results and the metadata

subroutine outputfiledef()

	integer :: ios
	character (len= 37) :: OutputFileName  
    
	do
		write (*,"(A)", ADVANCE="no") 'Enter output filename (without extension): '
		read (*,*) OutputFileName

		! creates sequential-access output file
		open (unit=18,file=trim(OutputFileName)//'_data.txt',status='NEW' &
  		, access='sequential', form='formatted', iostat=ios)

		if (ios /= 0) then
  			write (*,"(A,/,A)") 'Error with output file creation.','Change name'
  			cycle
		end if

		! creates sequential-access metadata file
		open (unit=19,file=trim(OutputFileName)//'_metadata.txt',status='NEW' &
  		, access='sequential', form='formatted', iostat=ios)

		if (ios /= 0) then
  			write (*,"(A,/,A)") 'Error with output file creation.','Change name'
  			cycle
		end if

        exit

	end do

    write (*,*)
    

end subroutine outputfiledef

!-------------------------


!------------------------

subroutine write_datafile_header(fault_loctimevar1)

type (record_variables) :: fault_loctimevar1

character(len=100) :: header_firstpart
character(len=150) :: header_result
    
select case (fault_loctimevar1%readcase)  
  case(0)
	write(header_firstpart,"(2A)") 'Id_rec',','
  case(1)
  if (fault_loctimevar1%spat_coord_system == 1) then
	write(header_firstpart,"(8A)") 'Id_rec',',','X',',','Y',',','Z',','
  elseif (fault_loctimevar1%spat_coord_system == 2) then
	write(header_firstpart,"(8A)") 'Id_rec',',','Lat',',','Lon',',','Z',','
  endif
  case(2)
  if (fault_loctimevar1%spat_coord_system == 1) then
	write(header_firstpart,"(10A)") 'Id_rec',',','X',',','Y',',','Z',',','Time',','
  elseif (fault_loctimevar1%spat_coord_system == 2) then
	write(header_firstpart,"(10A)") 'Id_rec',',','Lat',',','Lon',',','Z',',','Time',','
  endif
  case(3)
	write(header_firstpart,"(4A)") 'Id_rec',',','Time',','
end select 

write (header_result,"(19A)") 'Strike_rhr',',','Dip_angle',',','Th_rake',',','Th_trend_slick',',','Th_plunge_slick'   &
    ,',','Sigma_Magn',',','SigmaN_Magn',',','SigmaT_Magn',',','Mod_SlipTendency',',','Deformation_Index'

write (18,"(A)") trim(header_firstpart)//trim(header_result)

end subroutine write_datafile_header

!------------------------


!------------------------
! writes the stress analysis results in the data output file

subroutine write_datafile_result(fault_loctimevar1,faultrec1,stresssolution1)

type (record_variables) :: fault_loctimevar1
type(fault_datum) :: faultrec1
type(stress_solution) :: stresssolution1

character(len=100) :: result_firstpart
character(len=100) :: result_secondpart
character(len=200) :: result_thirdpart

select case (fault_loctimevar1%readcase)  
  case(0)
	write(result_firstpart,"(i9)") faultrec1%id 
  case(1)
	write(result_firstpart,"(i9,3(A,f15.5))") faultrec1%id,',',faultrec1%x,',',faultrec1%y,',',faultrec1%z
  case(2)
	write(result_firstpart,"(i9,4(A,f15.5))") faultrec1%id,',',faultrec1%x,',',faultrec1%y,',',faultrec1%z,',',faultrec1%t
  case(3)
	write(result_firstpart,"(i9,1(A,f15.5))") faultrec1%id,',',faultrec1%t
end select 

write(result_secondpart,"(2(A,f6.1))") ',',faultrec1%fltplane%strike_rhr,',',faultrec1%fltplane%dipangle	

if (stresssolution1%valid_solution) then     
	write (result_thirdpart,"(3(A,f6.1),5(A,f8.4))") ',',stresssolution1%theor_rake  &
    ,',',stresssolution1%theor_slickenline%trend,',',stresssolution1%theor_slickenline%plunge  &   
	,',',stresssolution1%tractionvect_magn,',',stresssolution1%normalstress_magn,',',stresssolution1%shearstress_magn  &
    ,',',stresssolution1%modif_sliptendency,',',stresssolution1%deformation_index 
else
	write (result_thirdpart,"(8(A))") ',',',',',',',',',',',',',',','
endif   

write (18,"(A)") trim(result_firstpart)//trim(result_secondpart)//trim(result_thirdpart)


end subroutine write_datafile_result

!------------------------


!------------------------

subroutine write_metadata(programparams1, stressparams1)

integer (kind=i2b) :: i,j
type(program_parameters) :: programparams1
type(stressanalysis_param) :: stressparams1
    
write(19,"(A)") 'Metadata for Homogeneous Stress Forward Simulation analysis' 
write(19,"(A)") 'Output of ForwardStress program'
write (19,"(A,1x,i4,A,i0,A,i0)", advance="no") 'Analysis of:',programparams1%time_initial(1),'/'    &
		,programparams1%time_initial(2),'/',programparams1%time_initial(3)
write (19,"(A,3(i0,A))") ' - starts at ',programparams1%time_initial(5)    &
				,':',programparams1%time_initial(6),':',programparams1%time_initial(7)
write (19,*)                
if (programparams1%input_case == 1) then
	write (19,"(3A)") 'Input: external file "'//trim(programparams1%inputfile_name)//'"'
elseif (programparams1%input_case == 2) then
	write (19,"(A)") 'Input: internal simulation'
endif
write (19,*)
write (19,"(A,i0)") 'Total number of faults: ', programparams1%faults_totnumb        
write (19,*)
write (19,*) 'Homogeneous stress field with:'
write (19,"(A,2(f6.1,2x))") '  S1 trend and plunge: '   &
	,stressparams1%stressprinccomp1%S1_pole%trend,stressparams1%stressprinccomp1%S1_pole%plunge
write (19,"(A,2(f6.1,2x))") '  S2 trend and plunge: '   &
	,stressparams1%stressprinccomp1%S2_pole%trend,stressparams1%stressprinccomp1%S2_pole%plunge
write (19,"(A,2(f6.1,2x))") '  S3 trend and plunge: '   &
	,stressparams1%stressprinccomp1%S3_pole%trend,stressparams1%stressprinccomp1%S3_pole%plunge
write (19,*)
write (19,"(A,3(f9.3,2x))") 'Sigma 1,2,3 magnitudes: '   &
	,stressparams1%stressprinccomp1%sigma1,stressparams1%stressprinccomp1%sigma2   &
    ,stressparams1%stressprinccomp1%sigma3
write (19,"(A,f7.5)") 'Phi value: ',stressparams1%stressprinccomp1%phi
write (19,*)
write (19,*) 'Rotation matrix components - R11,12,13,21,22,23,31,32,33'
do i=1,3
  do j=1,3
    write(19,"(f9.5)",advance="no") stressparams1%rot_matrix(i,j)
  end do
  write(19,*)
end do
write (19,*)
write (19,*) 'Stress tensor components - T11,12,13,21,22,23,31,32,33'
do i=1,3
  do j=1,3
    write(19,"(f9.5)",advance="no") stressparams1%tens(i,j)
  end do
  write(19,*)
end do

end subroutine write_metadata

!------------------------


end module inputoutput_defs

!------------------------------------------------ 



!------------------------------------------------ 

program ForwardStress


use inputoutput_defs

implicit none

integer (kind=i4b) :: i

type(program_parameters) :: programparams1
type(stressanalysis_param) :: stressparams1
type(fault_datum) :: faultrec1
type(stress_solution) :: stresssolution1
integer (kind=i4b) :: numfaults

pi = dacos(0.0d0)*2.0d0 ! pi radians
r2d = 180.d0/pi			! conversion from radians to degrees
d2r = pi/180.d0			! conversion from degrees to radians


!-------------------------------

! get initial time_date of analysis 
call DATE_AND_TIME(VALUES=programparams1%time_initial)


! print the initial message to screen
call initial_message()

                            
! user input of homogeneous stress field parameters
write(*,*) 'Definition of stress field characteristics'
stressparams1%stressprinccomp1 = define_princcompstresstensor()
write(*,*)

    
! define rotation matrix based on stress principal axes
! see Kuipers, 2002, eqs. 7.8 and 7.9
stressparams1%rot_matrix = rotmatr(stressparams1%stressprinccomp1)
! calculate stress tensor expressed in frame components
stressparams1%tens = stresstensorcalc(stressparams1%stressprinccomp1  &
							,stressparams1%rot_matrix)
  
                          
! define considered fault charateristics other than id, strike and dip angle
! i.e. spatial location and/or time
programparams1%fault_loctimevar = faultparam_define()


! manage data from input file or from simulation
programparams1%input_case = define_input_case()
select case (programparams1%input_case)  
  case(1)
	! define and open input file	
	call inputfiledef(programparams1)    
  case(2)
  	! define and generate simulation
    call simulationparamdefine(programparams1)
    write(*,*) ' generating simulation ....'    
    call generatesimulatedfaults(programparams1)    
	write(*,*) ' simulation successfully generated' 
end select 

write(*,*)

! define and open output file
call outputfiledef()

!  input definition phase completed
! result definition phase begins
    
! write header to output file 
call write_datafile_header(programparams1%fault_loctimevar)
    
! reading of data from input/temp file, calculation and output of results
write(*,*) ' calculating forward stress simulation ....'

! skip header
if (programparams1%inputfile_headerrownumb > 0) then
  do i=1,programparams1%inputfile_headerrownumb
	read (17,*)
  end do
endif

numfaults = 0
 
do

	select case (programparams1%fault_loctimevar%readcase)  
  	case(0)
   		read (17,*,iostat=ios) faultrec1%id,faultrec1%fltplane%strike_rhr,faultrec1%fltplane%dipangle
  	case(1)
		read (17,*,iostat=ios) faultrec1%id,faultrec1%x,faultrec1%y,faultrec1%z &
        ,faultrec1%fltplane%strike_rhr,faultrec1%fltplane%dipangle     	
  	case(2)
		read (17,*,iostat=ios) faultrec1%id,faultrec1%x,faultrec1%y,faultrec1%z,faultrec1%t &
        ,faultrec1%fltplane%strike_rhr,faultrec1%fltplane%dipangle
  	case(3)
		read (17,*,iostat=ios) faultrec1%id,faultrec1%t &
        ,faultrec1%fltplane%strike_rhr,faultrec1%fltplane%dipangle    	    
	end select

    ! if finished reading, program stops
	if (ios/=0) then                 
    	exit
    end if
    
	numfaults = numfaults + 1

	! calculation of fault dip direction
	call dipdir_calc(faultrec1%fltplane) 
       
	stresssolution1 = stresssolution_calc(stressparams1%tens,faultrec1)

	call write_datafile_result(programparams1%fault_loctimevar,faultrec1,stresssolution1)

end do

programparams1%faults_totnumb = numfaults

write(*,*) 'processing completed'

call write_metadata(programparams1, stressparams1)
  
! if source file is obtained through simulation,
! closes and deletes the temporary file
if (programparams1%input_case==2) then
	close(unit=17,status='delete')
endif

end program ForwardStress