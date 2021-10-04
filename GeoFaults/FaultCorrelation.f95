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

implicit none

type :: vector  
	sequence             
	real (kind=r8b)	:: x,y,z    
end type vector

type :: axis
	sequence 
	real (kind=r8b)	:: trend, plunge
end type axis

type :: orthonormal_triad
	sequence
	type(vector) :: X,Y,Z
end type orthonormal_triad

type :: triad_axes
	type(axis) :: axis_a, axis_b, axis_c    
end type triad_axes

type :: spherical_loc
	sequence
    real(kind=r8b) :: lat,lon, depth ! lat and lon in decimal degrees, depth in km 
end type spherical_loc
    
type :: plane
	sequence
	real (kind=r8b)	:: strike_rhr,dipdirection,dipangle
end type plane
    
type :: quaternion
	real(kind=r8b) :: q(0:3) ! quaternion components, the last is the rotation component
end type quaternion

type :: rotation
	type(axis) :: rot_axis
    real(kind=r8b) :: rot_angle
end type rotation

type :: fault_datum
	integer (kind=i4b) :: id
	type (plane) :: fltplane
	type (axis) :: slickenline    
	real (kind=r8b)	:: rake_aki
	type(vector) :: faultnorm_vect,slickenl_vect
	type(triad_axes) :: focmech ! 3 mutually orthogonal axes (as trend and plunge) representing T, P and B
    type(orthonormal_triad) :: fochmecv !X=T, Y=P, Z=B kinematic axes
    type(spherical_loc) :: sphercoord
    real(kind=r8b) :: spatloc(3) !x,y,z coordinates   
    real(kind=r8b) :: t
end type fault_datum


type :: fault_proc
	sequence
	integer (kind=i4b) :: id
    real(kind=r8b) :: time
    real(kind=r8b) :: spatloc(3) !x,y,z coordinates
	type(vector) :: faultnorm_vect,slickenl_vect
    real(kind=r8b) :: mom_tens(3,3)
    type(orthonormal_triad) :: fochmecv !X=T, Y=P, Z=B kinematic axes
end type fault_proc

    
! spatial information
type :: spatialdataformat
	logical :: present
	integer (i1b) :: cod   ! 0: cartesian (x,y,z); 1: spherical (lat,lon,depth)
end type spatialdataformat


! time information    
type :: timedataformat
	logical :: present
end type timedataformat


! orientation information    
type :: orientationdataformat
	integer (i1b) :: cod  ! 0 - strike dip rake; 1 - strike dip slick. trend & plunge; 2 - P-axis trend & plunge T-axis trend & plunge
end type orientationdataformat


! data format 
type :: dataformat
	type(spatialdataformat) :: spatial
    type(timedataformat) :: time
    type(orientationdataformat) :: orientation 
    character (len=6) :: cod   
end type dataformat


! pairwise or central statistics 
type :: statistics_type
	integer (i1b) :: cod  ! 0 - pairwise; 1 - central statistics
	integer (i1b) :: centrflt_orientcod  ! same codes as orientation information 
    character (len=6) :: centrflt_generalcod    
	type (fault_datum)  ::  centrfaultdat
end type statistics_type

 
! settings of analysis session
type :: analysistype
	type(statistics_type) :: statistics
	logical :: focmechrot		! focal mechanism rotation statistics
	logical :: kinematicaxesrot	! T-P-B axes rotations
	logical :: faultelementsrot	! fault sub-element rotation statistics
	logical :: coherence		! coherence statistics
	logical :: spatialsep		! spatial separation statistics
	logical :: timelag			! time lag statistics    
end type analysistype


! parameters of file
type :: file_params
	character (len=255) :: name
    integer (kind=i1b) :: statuscod ! 0: replace
	integer :: unitnum
    integer (kind=i2b) :: headerrownumber
	integer (kind=i4b) :: recsnumber, reclength
end type file_params


! analysis results
type :: analysis_results
	type(rotation) :: focmech_rots(4), faultpl_rot, slickenl_rot
	type(rotation) :: Taxis_rot, Paxis_rot, Baxis_rot
    real (kind=r8b) :: coherence, time_lag, spat_sep
end type analysis_results


! constants declaration
real (kind=r8b) :: pi ! pi radians
real (kind=r8b) :: r2d ! conversion from radians to degrees
real (kind=r8b) :: d2r ! conversion from degrees to radians

! declaration of Earth ellipsoid parameters
! from Longley et al., 2001, table 4.2, WGS84 ellipsoid case
real (kind=r8b), parameter :: Earth_major_radius = 6378.137 !km
! from http://home.online.no/~sigurdhu/WGS84_Eng.html, and also eq. 6.29, p. 225 in Leick, 1995 
real (kind=r8b), parameter :: eccentricity_sq = 0.00669437999014 



end module data_types


!------------------------------------------------    



!------------------------------------------------    

module vector_processing


use data_types

logical :: vect_is_normalized
real (kind=r8b):: scaledscalarproduct, vector1_magn, vector2_magn
real (kind=r8b),parameter :: unitary_vect_tolerance = 1.0e-6 ! maximum accepted value for unitary vector magnitude difference with unit
real (kind=r8b),parameter :: vect_normaliz_tolerance = 1.0e-6 ! minimum accepted value for vector magnitude to apply normalization


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
! vector1 projection on vector2

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




!----------------------------------------------


module geometric_processing

use vector_processing



contains


!------------------------
! calculates the cartesian components from the polar ones
! for a unit sphere

type(vector) function pole2cartesian(axis1) result(vector1)

!  QUATg_pole2cartesian
!  created 2002-12-30

type(axis):: axis1


vector1%x = cos(d2r*axis1%plunge)*cos(d2r*axis1%trend)
vector1%y  = cos(d2r*axis1%plunge)*sin(d2r*axis1%trend)
vector1%z  =  sin(d2r*axis1%plunge)


end function pole2cartesian

!------------------------


!------------------------
! calculates the cartesian components from the polar ones
! for an Earth-radius sphere

function pole2earthcartesian(spherloc1) result(spatloc1)

!  QUATg_pole2cartesian
!  created 2002-12-30

type(spherical_loc), intent(in):: spherloc1
real (kind=r8b) :: spatloc1(3)

real (kind=r8b) :: curvature_rad 


curvature_rad = Earth_major_radius/dsqrt(1-eccentricity_sq*(dsin(d2r*spherloc1%lat))**2)

spatloc1(1) = (curvature_rad-spherloc1%depth)*dcos(d2r*spherloc1%lat)*dcos(d2r*spherloc1%lon)
spatloc1(2) = (curvature_rad-spherloc1%depth)*dcos(d2r*spherloc1%lat)*dsin(d2r*spherloc1%lon)
spatloc1(3) = (curvature_rad*(1-eccentricity_sq)-spherloc1%depth)*dsin(d2r*spherloc1%lat)


end function pole2earthcartesian

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
    elseif (axis1%trend >= 360.0) then
		axis1%trend = axis1%trend - 360.0    
    endif
endif


end function cartesian2pole

!------------------------


!------------------------
! calculates the anti-pole to a given pole
type (axis) function antipole(axis1)

!  QUATg_pole2antipole
!  created 2005-03-26

type(axis) :: axis1


antipole%plunge = - axis1%plunge
antipole%trend = axis1%trend + 180.0

if (antipole%trend >= 360) then
  antipole%trend = antipole%trend - 360.0
endif

end function antipole

!------------------------


!------------------------
! calculates the down axis from the axis

type(axis) function axis2downaxis(axis1) result(axis2)

!  QUATg_axis2downaxis
!  created 2003-02-01

type(axis) :: axis1

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


!--------------------------
! calculates the vector normal to a given plane 
 
type(vector) function plane_normal(plane1) result(normaltoplane1)

type(plane) :: plane1

! formulas from Aki and Richards, 1980
normaltoplane1%x = -sin(d2r*plane1%dipangle)*sin(d2r*plane1%strike_rhr)
normaltoplane1%y = sin(d2r*plane1%dipangle)*cos(d2r*plane1%strike_rhr)
normaltoplane1%z = -cos(d2r*plane1%dipangle)

end function plane_normal

!--------------------------


!--------------------------
! checks if an axis (given as trend & plunge) lies on plane

logical function isaxisonplane(plane1,axis1) 

type(plane), intent(in) :: plane1
type(axis), intent(in) :: axis1
type(vector) :: plane_normal_vect1, axis_vect1

plane_normal_vect1 = plane_normal(plane1)
axis_vect1 = pole2cartesian(axis1)

anglebetweenaxes = r2d*axes_angle_rad(plane_normal_vect1,axis_vect1)
if (anglebetweenaxes >89.0) then
	isaxisonplane = .true.
else
	isaxisonplane = .false.
endif

end function isaxisonplane

!--------------------------


!--------------------------
! calculates the 3D euclidean distance between two points 
  
real(kind=r8b) function euclidean_distance(point1, point2)

real (kind=r8b), dimension(3) :: point1, point2

euclidean_distance = dsqrt((point2(1)-point1(1))**2 + (point2(2)-point1(2))**2 + (point2(3)-point1(3))**2)


end function euclidean_distance

!--------------------------



end module geometric_processing


!----------------------------------------------




!------------------------------------------------    


module quaternion_processing

use geometric_processing

implicit none

logical :: quat_is_normalized
real (kind=r8b) :: quat_sqnorm, quat_norm
real (kind=r8b),parameter :: quat_normaliz_tolerance = 1.0e-6


type (triad_axes) ::  focmech_1, focmech_2, central_focalmechanism
real (kind=r8b) :: focmec1_matrix(3,3), focmec2_matrix(3,3)



contains


!--------------------------


type(quaternion) function quat_product(quat1,quat2) result(quat_prod)

! QUATq_product
! quaternion product
! Avenue created: 2005-02-11

type(quaternion) :: quat1,quat2


quat_prod%q(0) =  (quat1%q(0)*quat2%q(0))   &
                 -(quat1%q(1)*quat2%q(1))	&
                 -(quat1%q(2)*quat2%q(2))	&
                 -(quat1%q(3)*quat2%q(3))
                      
quat_prod%q(1) =  (quat1%q(0)*quat2%q(1))	&
                 +(quat1%q(1)*quat2%q(0))	&
                 +(quat1%q(2)*quat2%q(3))	&
                 -(quat1%q(3)*quat2%q(2))

quat_prod%q(2) =  (quat1%q(0)*quat2%q(2))	&
                 -(quat1%q(1)*quat2%q(3))	&
                 +(quat1%q(2)*quat2%q(0))	&
                 +(quat1%q(3)*quat2%q(1))

quat_prod%q(3) =  (quat1%q(0)*quat2%q(3))	&
                 +(quat1%q(1)*quat2%q(2))	&
                 -(quat1%q(2)*quat2%q(1))	&
                 +(quat1%q(3)*quat2%q(0))              
         

end function quat_product

!--------------------------


!--------------------------


type(quaternion) function quat_conjugate(quat1) result(quat_conj)

! QUATq_conjugate
! created 2005-02-10

type(quaternion) :: quat1

quat_conj%q(0) =   quat1%q(0)
quat_conj%q(1) = - quat1%q(1)
quat_conj%q(2) = - quat1%q(2)
quat_conj%q(3) = - quat1%q(3)

end function quat_conjugate

!--------------------------


!--------------------------


real (kind=r8b) function quat_squarednorm(quat1)

!  QUATq_squarednorm
!  created 2005-02-12

type(quaternion) :: quat1

quat_squarednorm = (quat1%q(0)**2)+(quat1%q(1)**2)+(quat1%q(2)**2)+(quat1%q(3)**2)

end function quat_squarednorm

!--------------------------


!--------------------------


type(quaternion) function quat_scalardivision(quat1, scaldiv) result(quat_scaldiv)

! QUATq_divisionbyscalar
! quaternion division by a scalar
! created: 2005-02-12

integer (kind=i1b) :: i
type(quaternion), intent(in) :: quat1
real (kind=r8b), intent(in) :: scaldiv

do i=0,3
	quat_scaldiv%q(i) = quat1%q(i)/scaldiv
end do

end function quat_scalardivision

!--------------------------


!--------------------------

type(quaternion) function quat_inverse(quat1) result(quat_inv)

!  quaternion inverse
!  created: 2005-02-19

type(quaternion) :: quat1, quat_conj

quat_conj = quat_conjugate(quat1)
quat_sqnorm = quat_squarednorm(quat1)                  
quat_inv = quat_scalardivision(quat_conj, quat_sqnorm) 


end function quat_inverse

!--------------------------


!--------------------------

logical function quat_normaliztest(quat1) result(quat_is_normalized)

! QUATq_normalizedquaterniontest
! created 2005-02-12

real (kind=r8b) :: abs_diff
type(quaternion) :: quat1

quat_sqnorm = quat_squarednorm(quat1)                                       
abs_diff = dabs(1-quat_sqnorm)
                           
if (abs_diff > quat_normaliz_tolerance) then
  quat_is_normalized = .false.
else
  quat_is_normalized = .true.
endif

end function quat_normaliztest

!--------------------------


!--------------------------

type(quaternion) function quat_normalization(quat1) result(quat_normalized)

!  QUATq_normalization
!  transformation from quaternion to normalized quaternion
!  created: 2005-02-14

type(quaternion) :: quat1


quat_sqnorm = quat_squarednorm(quat1)
quat_norm = dsqrt(quat_sqnorm)
quat_normalized = quat_scalardivision(quat1, quat_norm)

quat_is_normalized = quat_normaliztest(quat_normalized)
if (.not.quat_is_normalized) then
  write(*,*) 'Error in quaternion normalization. Hit any key to stop'
  read(*,*)
  stop
endif

end function quat_normalization

!--------------------------



!--------------------------
! converts orthonormal triad components into a 3x3 matrix

function cartesmatrix(orthotriad1) result (matrix33)

! QUATq_TPBcomp2cartesmatrix
! created: 2005-02-17

real (kind=r8b):: matrix33(3,3) 

type(orthonormal_triad) :: orthotriad1


matrix33(1,1) = orthotriad1%X%x    
matrix33(1,2) = orthotriad1%Y%x 
matrix33(1,3) = orthotriad1%Z%x
matrix33(2,1) = orthotriad1%X%y
matrix33(2,2) = orthotriad1%Y%y
matrix33(2,3) = orthotriad1%Z%y
matrix33(3,1) = orthotriad1%X%z
matrix33(3,2) = orthotriad1%Y%z
matrix33(3,3) = orthotriad1%Z%z


end function cartesmatrix

!--------------------------


!--------------------------
! calculates a quaternion from a 3x3 matrix

type(quaternion) function quaternfromcartmatr(focmec1_matrix) result(quat1)

! QUATq_TPBcartesmatrix2quaternion
! modified 2005-02-17


real (kind=r8b) :: focmec1_matrix(3,3)
real (kind=r8b) :: Q0, Q1, Q2, Q3
real (kind=r8b) :: Q0Q1,Q0Q2,Q0Q3,Q1Q2,Q1Q3,Q2Q3

! myR11 = t1 = focmec1_matrix(1,1)
! myR21 = t2 = focmec1_matrix(2,1)
! myR31 = t3 = focmec1_matrix(3,1)
! myR12 = p1 = focmec1_matrix(1,2)
! myR22 = p2 = focmec1_matrix(2,2)
! myR32 = p3 = focmec1_matrix(3,2)
! myR13 = b1 = focmec1_matrix(1,3)
! myR23 = b2 = focmec1_matrix(2,3)
! myR33 = b3 = focmec1_matrix(3,3)


Q0 = 0.5*(dsqrt(1+focmec1_matrix(1,1)+focmec1_matrix(2,2)+focmec1_matrix(3,3)))
Q1 = 0.5*(dsqrt(1+focmec1_matrix(1,1)-focmec1_matrix(2,2)-focmec1_matrix(3,3)))
Q2 = 0.5*(dsqrt(1-focmec1_matrix(1,1)+focmec1_matrix(2,2)-focmec1_matrix(3,3)))
Q3 = 0.5*(dsqrt(1-focmec1_matrix(1,1)-focmec1_matrix(2,2)+focmec1_matrix(3,3)))

Q0Q1 = 0.25*(focmec1_matrix(3,2) - focmec1_matrix(2,3))
Q0Q2 = 0.25*(focmec1_matrix(1,3) - focmec1_matrix(3,1))
Q0Q3 = 0.25*(focmec1_matrix(2,1) - focmec1_matrix(1,2))
Q1Q2 = 0.25*(focmec1_matrix(1,2) + focmec1_matrix(2,1))
Q1Q3 = 0.25*(focmec1_matrix(1,3) + focmec1_matrix(3,1))
Q2Q3 = 0.25*(focmec1_matrix(2,3) + focmec1_matrix(3,2))

if((3*Q0)>(Q1+Q2+Q3)) then
  Q1 = Q0Q1/Q0
  Q2 = Q0Q2/Q0
  Q3 = Q0Q3/Q0
elseif ((3*Q1)>(Q0+Q2+Q3)) then
  Q0 = Q0Q1/Q1
  Q2 = Q1Q2/Q1
  Q3 = Q1Q3/Q1
elseif ((3*Q2)>(Q0+Q1+Q3)) then
  Q0 = Q0Q2/Q2
  Q1 = Q1Q2/Q2
  Q3 = Q2Q3/Q2
else
  Q0 = Q0Q3/Q3
  Q1 = Q1Q3/Q3
  Q2 = Q2Q3/Q3
end if

quat1%q(0)= Q0
quat1%q(1)= Q1
quat1%q(2)= Q2
quat1%q(3)= Q3

end function quaternfromcartmatr

!--------------------------


!--------------------------
! calculates the rotation from a quaternion

type(rotation) function rotationaxis(quat1) result(rot1)

!  QUATq_quaternion2rotationpole
!  transformation from quaternion to rotation pole
!  created: 2005-02-12
!  modified: 2005-02-15


type(vector) :: vector1
type(quaternion) :: quat1

quat1 = quat_normalization(quat1)

vector1%x = quat1%q(1)
vector1%y = quat1%q(2)
vector1%z = quat1%q(3)

vector1_magn = vector_magn(vector1)

if (vector1_magn < 0.000001) then   
	rot1%rot_axis%trend = 0.0
	rot1%rot_axis%plunge = 0.0
	rot1%rot_angle = 0.0 
else 
	rot1%rot_angle = 2*r2d*dacos(quat1%q(0))
	if (rot1%rot_angle > 180.0) then
		rot1%rot_angle = -(360.0 - rot1%rot_angle)
	endif

	vector1 = vector_normalization(vector1)
	rot1%rot_axis = cartesian2pole(vector1)

	if(rot1%rot_axis%plunge < 0.0) then
		rot1%rot_axis = axis2downaxis(rot1%rot_axis)
		rot1%rot_angle = -rot1%rot_angle
	endif
endif

  
end function rotationaxis

!------------------------


end module quaternion_processing


!----------------------------------------------



!----------------------------------------------


module fault_processing

use quaternion_processing

implicit none

contains


!--------------------------
! calculates slikenline vector from strike, dip, rake

type(vector) function faultrake2slick_vector(strike_rhr,dipangle,rake_aki) result(slick_vect) 

! QUATf_faultpole2faultvector

real (kind=r8b) :: strike_rhr,dipangle,rake_aki

! formulas from Aki and Richards, 1980
slick_vect%x = cos(d2r*rake_aki)*cos(d2r*strike_rhr)+sin(d2r*rake_aki)*cos(d2r*dipangle)*sin(d2r*strike_rhr)
slick_vect%y = cos(d2r*rake_aki)*sin(d2r*strike_rhr)-sin(d2r*rake_aki)*cos(d2r*dipangle)*cos(d2r*strike_rhr)
slick_vect%z = -sin(d2r*rake_aki)*sin(d2r*dipangle)


end function faultrake2slick_vector

!--------------------------


!--------------------------
! converts from slickenline (vector) to rake
 
subroutine slickenline2rake(fault_rec1)

type(fault_datum), intent(inout) :: fault_rec1

type(axis) :: strike_axis
type(vector) :: strike_vect
real (kind=r8b) :: lambda_angle

strike_axis%trend = fault_rec1%fltplane%strike_rhr
strike_axis%plunge = 0

strike_vect = pole2cartesian(strike_axis)

lambda_angle = r2d*vector_angle_rad(strike_vect,fault_rec1%slickenl_vect)

if (fault_rec1%slickenl_vect%z < 0.0) then
  fault_rec1%rake_aki = lambda_angle
else
  fault_rec1%rake_aki = -lambda_angle
endif

end subroutine slickenline2rake

!--------------------------


!--------------------------
! calculates T-P-B components from fault normal and slickenline vectors

type(orthonormal_triad) function TPBvectors_calc1(faultnormal_vect1, slickenline_vect1) result(TPB_vectors1)

	! QUATf_faultvector2TPBcomp


	type(vector), intent(in) :: faultnormal_vect1, slickenline_vect1


	TPB_vectors1%X = vector_sum(faultnormal_vect1, slickenline_vect1)                    
	TPB_vectors1%X = vector_normalization(TPB_vectors1%X)

	TPB_vectors1%Y = vector_diff(faultnormal_vect1, slickenline_vect1)                    
	TPB_vectors1%Y = vector_normalization(TPB_vectors1%Y)

	TPB_vectors1%Z = vector_vectprod(TPB_vectors1%X,TPB_vectors1%Y)

    

end function TPBvectors_calc1

!--------------------------


!------------------------
!  incremental strain tensor calculation from strike, dip and rake

function incremstraintensorcalc0(strike_rhr1,dipangle1,rake_aki1) result(incrstraintens)

!  QUATc_momenttensorcalc
!  Created: 2004-04-03

real (kind=r8b), intent(in) :: strike_rhr1,dipangle1,rake_aki1
real (kind=r8b) :: incrstraintens(3,3)

!  formulas taken from: Aki & Richards, 1980, vol. 1, Box 4.4

incrstraintens(1,1) = - sin(d2r*dipangle1)*cos(d2r*rake_aki1)*sin(2*d2r*strike_rhr1) &
					  - sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)*sin(d2r*strike_rhr1)**2
incrstraintens(1,2) = sin(d2r*dipangle1)*cos(d2r*rake_aki1)*cos(2*d2r*strike_rhr1) &
					+ 0.5*sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)*sin(2*d2r*strike_rhr1) 
incrstraintens(1,3) = - cos(d2r*dipangle1)*cos(d2r*rake_aki1)*cos(d2r*strike_rhr1)  &
					 - cos(2*d2r*dipangle1)*sin(d2r*rake_aki1)*sin(d2r*strike_rhr1)
incrstraintens(2,1) = incrstraintens(1,2) 
incrstraintens(2,2) = sin(d2r*dipangle1)*cos(d2r*rake_aki1)*sin(2*d2r*strike_rhr1) &
					- sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)*cos(d2r*strike_rhr1)**2
incrstraintens(2,3) = - cos(d2r*dipangle1)*cos(d2r*rake_aki1)*sin(d2r*strike_rhr1) &
						+ cos(2*d2r*dipangle1)*sin(d2r*rake_aki1)*cos(d2r*strike_rhr1)
incrstraintens(3,1) = incrstraintens(1,3)
incrstraintens(3,2) = incrstraintens(2,3)
incrstraintens(3,3) = sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)

end function incremstraintensorcalc0

!------------------------


!------------------------
!  incremental strain tensor calculation from focal mechanisms

function incremstraintensorcalc1(Tvect1,Pvect1) result(incrstraintens)

!  Created: 2008-01-27


type(vector), intent(in) :: Tvect1, Pvect1
real (kind=r8b) :: incrstraintens(3,3)

!  formulas taken from: Kagan & Knopoff, 1985b, p. 649 

incrstraintens(1,1) = (Tvect1%x + Pvect1%x)*(Tvect1%x - Pvect1%x)
incrstraintens(1,2) = (Tvect1%x + Pvect1%x)*(Tvect1%y - Pvect1%y)
incrstraintens(1,3) = (Tvect1%x + Pvect1%x)*(Tvect1%z - Pvect1%z)
incrstraintens(2,1) = (Tvect1%y + Pvect1%y)*(Tvect1%x - Pvect1%x) 
incrstraintens(2,2) = (Tvect1%y + Pvect1%y)*(Tvect1%y - Pvect1%y)
incrstraintens(2,3) = (Tvect1%y + Pvect1%y)*(Tvect1%z - Pvect1%z)
incrstraintens(3,1) = (Tvect1%z + Pvect1%z)*(Tvect1%x - Pvect1%x)
incrstraintens(3,2) = (Tvect1%z + Pvect1%z)*(Tvect1%y - Pvect1%y)
incrstraintens(3,3) = (Tvect1%z + Pvect1%z)*(Tvect1%z - Pvect1%z)

end function incremstraintensorcalc1

!------------------------

  
!------------------------
! calculates the focal mechanisms vectors from a strike-dip-rake input

type(orthonormal_triad) function fochalmechanismvectorcalc0(strike_rhr1, dipangle1, rake_aki1) result(TPBvect1)

type(vector) :: faultnormal_vect1, slickenline_vect1
type(plane) :: faultplane1
real (kind=r8b) :: strike_rhr1, dipangle1, rake_aki1

faultplane1%strike_rhr = strike_rhr1
faultplane1%dipangle = dipangle1


faultnormal_vect1 = plane_normal(faultplane1)
slickenline_vect1 = faultrake2slick_vector(strike_rhr1, dipangle1, rake_aki1) 

TPBvect1 =  TPBvectors_calc1(faultnormal_vect1, slickenline_vect1)

end function fochalmechanismvectorcalc0

!------------------------


end module fault_processing



!------------------------------------------------ 


module similarity_processing

use fault_processing

implicit none



contains


!------------------------
! calculates coherence between fault/focal mechanism pairs

real (kind=r8b) function coherence_solution(mom_tens1,mom_tens2) result(coherence)

!  formulas taken from: Kagan & Knopoff, 1985a,b

real(kind=r8b), dimension(3,3) :: mom_tens1,mom_tens2

coherence = sum(mom_tens1*mom_tens2)


end function coherence_solution

!------------------------

    
!------------------------
! calculate the rotation solutions between a couple of triadic vectors
    
function triadvectors_rotsolution(orthotriadvect_1,orthotriadvect_2) result(rotation_solution)
    

type (orthonormal_triad), intent(in) ::  orthotriadvect_1, orthotriadvect_2
type(rotation) :: rotation_solution(4)

integer (kind=i1b) :: i
real (kind=r8b) :: angle_between_X_axes,angle_between_Y_axes
type(quaternion) :: fm1_quaternion, fm2_quaternion
type(quaternion) :: fm1_inversequatern, fm2_inversequatern
type(quaternion) :: rotationquatern(4)
type(quaternion) :: suppl_quat(3), suppl_prodquat(3), suppl_prod2quat(3)


! processing of equal focal mechanisms
angle_between_X_axes = r2d*axes_angle_rad(orthotriadvect_1%X,orthotriadvect_2%X)
angle_between_Y_axes = r2d*axes_angle_rad(orthotriadvect_1%Y,orthotriadvect_2%Y)

if((angle_between_X_axes<0.5).and.(angle_between_Y_axes<0.5)) then
  do i = 1,4
	rotation_solution(i)%rot_angle = 0.0
    rotation_solution(i)%rot_axis%trend = 0.0
    rotation_solution(i)%rot_axis%plunge = 0.0
  end do
  return
endif

  
! transformation of XYZ axes cartesian components (fm1,2) into quaternions q1,2

focmec1_matrix = cartesmatrix(orthotriadvect_1)
focmec2_matrix = cartesmatrix(orthotriadvect_2)

fm1_quaternion = quaternfromcartmatr(focmec1_matrix)
fm2_quaternion = quaternfromcartmatr(focmec2_matrix)


! calculation of quaternion inverse q1,2[-1]

fm1_inversequatern = quat_inverse(fm1_quaternion) 
fm2_inversequatern = quat_inverse(fm2_quaternion) 


! calculation of rotation quaternion : q' = q2*q1[-1]

rotationquatern(1) = quat_product(fm2_quaternion, fm1_inversequatern) 


! calculation of secondary rotation pure quaternions: a(i,j,k) = q2*(i,j,k)*q2[-1]
do i=1,3
	suppl_quat(i)%q = 0.0
end do
do i=1,3
	suppl_quat(i)%q(i) = 1.0
end do

do i=1,3
	suppl_prodquat(i) = quat_product(suppl_quat(i), fm2_inversequatern)
end do

do i=1,3
	suppl_prod2quat(i) = quat_product(fm2_quaternion, suppl_prodquat(i))
end do

! calculation of the other 3 rotation quaternions: q'(i,j,k) = a(i,j,k)*q'
do i=2,4
	rotationquatern(i) = quat_product(suppl_prod2quat(i-1), rotationquatern(1))
end do

! calculation of 4 quaternion polar parameters
do i=1,4
	rotation_solution(i) = rotationaxis(rotationquatern(i))
end do

! determination of minimum rotation solution
rotation_solution = sortrotationsolutions(rotation_solution)


end function triadvectors_rotsolution

!------------------------


!------------------------
! sorting of rotation solutions based on the magnitudes of the rotation angles

function sortrotationsolutions(rotationsolution) result(sortedrotationsolution)

type(rotation), intent(in) :: rotationsolution(4)
type(rotation) :: sortedrotationsolution(4)

integer (kind=i1b):: i, sortindex(4), sortndx(1)


real (kind=r8b) :: absrotang(4)


do i =1,4
	absrotang(i) = abs(rotationsolution(i)%rot_angle)
end do

do i =1,4
	sortndx = minloc(absrotang)
	sortindex(i) = sortndx(1)
    absrotang(sortindex(i)) = 999.9
end do

do i =1,4
	sortedrotationsolution(i)%rot_angle = rotationsolution(sortindex(i))%rot_angle
	sortedrotationsolution(i)%rot_axis%trend = rotationsolution(sortindex(i))%rot_axis%trend
	sortedrotationsolution(i)%rot_axis%plunge = rotationsolution(sortindex(i))%rot_axis%plunge
end do


end function sortrotationsolutions

!------------------------




end module similarity_processing


!------------------------------------------------    




!------------------------------------------------    


module inputoutput_defs


use similarity_processing

implicit none

type(fault_proc) :: faultproc1
 

contains


!------------------------
! prints initial message to screen
subroutine initial_message()

    write (*,*)
	write (*,"(10x,A)")'----------------------------------------------------------'
	write (*,"(25x,A)") 'Similarity Analysis'
	write (*,"(25x,A)") 'vers. 1.0, 2008-07'
	write (*,"(10x,A)")'----------------------------------------------------------'
    write (*,*)
    write (*,"(5x,A)") ' Program for analysing similarity between faults/focal mechanisms'
    write (*,"(5x,A)") ' Based on Kagan and Knopoff (1985a,b), Kagan (1990,1991,1992a,b)'      
    write (*,*)
	write (*,"(10x,A)")'----------------------------------------------------------'
    write (*,*)   
    
end subroutine initial_message

!------------------------

!------------------------
! define input spatial information

subroutine defineinput_spatialinfo(inputspatialformat)

implicit none

type (spatialdataformat), intent(out) :: inputspatialformat
character :: dummy_input


write(*,*) 
    
do
  	
	write(*,"(A)", ADVANCE="no") 'Spatial information present? (y/n): '
	read(*,*) dummy_input

    if (dummy_input == 'y') then
       inputspatialformat%present = .true.
    elseif (dummy_input == 'n') then
       inputspatialformat%present = .false.
    else
      cycle
    endif

    exit

end do

! select acceptable input format for spatial information
if (inputspatialformat%present) then

  write(*,*)	
  write(*,"(A)") 'Spatial information:'
  write(*,"(1x,A)") '0 - x, y, z (all values in the same units, e.g. m, km)'
  write(*,"(1x,A)") '1 - lat, lon, depth (lat & lon in decimal degrees, depth in km)' 
  do
	write(*,"(A)",advance= "no") 'Enter code for spatial data format (0/1): '
	read (*,*) inputspatialformat%cod
	if(inputspatialformat%cod == 0 .or. inputspatialformat%cod == 1) exit    
  end do

endif 
  
  
end subroutine defineinput_spatialinfo

!------------------------


!------------------------
! define input time information

subroutine defineinput_timeinfo(inputtimeformat)

implicit none

type (timedataformat), intent(out) :: inputtimeformat
character :: dummy_input

write(*,*) 

do
  
	write(*,"(A)", ADVANCE="no") 'Time information present? (y/n): '
	read(*,*) dummy_input

    if (dummy_input == 'y') then
       inputtimeformat%present = .true.
    elseif (dummy_input == 'n') then
       inputtimeformat%present = .false.
    else
      cycle
    endif

    exit

end do
 
  
end subroutine defineinput_timeinfo

!------------------------


!------------------------
! orientation information

subroutine defineinput_orientationinfo(inputorientationformat)

implicit none

type (orientationdataformat), intent(out) :: inputorientationformat

write(*,*) 
write(*,"(A)") 'Orientation data'
write(*,"(4x,A)") 'codes:'
write(*,"(5x,A)") '0 - strike dip rake (e.g. 312 32 -142)'
write(*,"(5x,A)") '1 - strike dip slick. trend & plunge (e.g. 312 32 132 0)'
write(*,"(5x,A)") '2 - P-axis trend & plunge T-axis trend & plunge (e.g. 200 37 290 0)'
    
do 
  
	write(*,"(A)",advance="no") 'Enter code (0/1/2): '  
	read (*,*) inputorientationformat%cod
	if(inputorientationformat%cod == 0 .or. inputorientationformat%cod == 1 .or. inputorientationformat%cod == 2) exit    

end do


end subroutine defineinput_orientationinfo

!------------------------


!------------------------
! program codes the input file format

subroutine inputfile_coding(inputdataformat)

type (dataformat), intent(inout) :: inputdataformat

character (len=2) :: space_code, time_code, orientation_code


! read code for spatial info
if (inputdataformat%spatial%present) then
	write(space_code,"(A,i1)") '1',inputdataformat%spatial%cod 
else
	space_code = "00"  
endif

! read code for time info
if (inputdataformat%time%present) then
	time_code = "10" 
else
	time_code = "00"  
endif

! read code for orientation info
write(orientation_code,"(A,i1)") '1',inputdataformat%orientation%cod 

inputdataformat%cod = space_code//time_code//orientation_code


end subroutine inputfile_coding

!------------------------



!------------------------
!define if pairwise or central statistics calculation

subroutine defineoutput_statisticsinfo(inputdataformat, analysissettings)

implicit none

type(dataformat), intent(in) :: inputdataformat
type (analysistype), intent(inout) :: analysissettings

	write(*,"(A)")
	write(*,"(A)") 'Code for statistic analysis: '
	write(*,"(5x,A)") '0: pairwise'
	write(*,"(5x,A)") '1: central statistics'
    
do
	write(*,"(A)",advance="no") 'Enter code (0/1): '
    read(*,*) analysissettings%statistics%cod  
	if (analysissettings%statistics%cod==0 .or. analysissettings%statistics%cod==1) exit
end do

if (analysissettings%statistics%cod==1) then

 call definecentralfaultdatum(inputdataformat, analysissettings)
  
endif

end subroutine defineoutput_statisticsinfo

!------------------------


!------------------------
! define the central fault datum  

subroutine definecentralfaultdatum(inputdataformat, analysissettings)

type(dataformat), intent(in) :: inputdataformat
type (analysistype), intent(inout) :: analysissettings

real (kind=r8b) :: angle_between_t0p0
type (vector) :: vector_a, vector_b


! id code for central datum (automatic definition)
analysissettings%statistics%centrfaultdat%id = 0


! user defines the orientation of the central datum
write (*,*)
write (*,"(1x,A)") 'Central statistic: parameters of central fault datum'
  

! user inputs the spatial coordinates for the central datum
if (analysissettings%spatialsep) then

	if (inputdataformat%spatial%cod == 0) then
		write (*,"(3x,A)", ADVANCE="no") 'X-Y-Z coordinates: '
		read (*,*) analysissettings%statistics%centrfaultdat%spatloc
    elseif (inputdataformat%spatial%cod == 1) then
		write (*,"(3x,A)", ADVANCE="no") 'Lat-lon-depth values: '
		read (*,*) analysissettings%statistics%centrfaultdat%sphercoord
    endif

endif


! user inputs the time coordinates for the central datum
if (analysissettings%timelag) then

		write (*,"(3x,A)", ADVANCE="no") 'Time: '
		read (*,*) analysissettings%statistics%centrfaultdat%t

endif



! choice of the format for central datum

    !if input fault orientation different from P-T axes, gives the user the possibility to choose P-T axes for central datum input
if (analysissettings%faultelementsrot) then
  do
  	write (*,"(5x,A)") '0 : strike-dip-rake'
    write (*,"(5x,A)") '1 : strike - dip - slickenline trend & plunge'
	read (*,*) analysissettings%statistics%centrflt_orientcod
    if (analysissettings%statistics%centrflt_orientcod == 0 .or. analysissettings%statistics%centrflt_orientcod == 1) exit
  end do
elseif (inputdataformat%orientation%cod /= 2) then 
	do
		write (*,"(5x,A)") '0 : strike-dip-rake'
		write (*,"(5x,A)") '1 : strike - dip - slickenline trend & plunge'
		write (*,"(5x,A)") '2 : P-T axes'
		read (*,*) analysissettings%statistics%centrflt_orientcod
        if (analysissettings%statistics%centrflt_orientcod == 0   .or.  &
            analysissettings%statistics%centrflt_orientcod == 1   .or.  &
            analysissettings%statistics%centrflt_orientcod == 2) then
            exit
        endif
    end do

else

  analysissettings%statistics%centrflt_orientcod = 2

endif


! the program formats the data code for the central datum as the data one
! (it will be updated for different orientation input format in the next if ..  
analysissettings%statistics%centrflt_generalcod=inputdataformat%cod


! user enters specific values for the orientation of the central datum,
! and the corresponding format code for the orientation is updated 
if (analysissettings%statistics%centrflt_orientcod == 0) then  !  0  -  strike dip rake

	analysissettings%statistics%centrflt_generalcod(6:6)='0'
    
	write (*,"(3x,A)", ADVANCE="no") 'Strike: '
	read (*,*) analysissettings%statistics%centrfaultdat%fltplane%strike_rhr        

	write (*,"(3x,A)", ADVANCE="no") 'Dip angle: '
	read (*,*) analysissettings%statistics%centrfaultdat%fltplane%dipangle  
    
	write (*,"(3x,A)", ADVANCE="no") 'Rake angle: '
	read (*,*) analysissettings%statistics%centrfaultdat%rake_aki
    
elseif  (analysissettings%statistics%centrflt_orientcod == 1) then  !  1  -  strike dip slick. trend & plunge

    analysissettings%statistics%centrflt_generalcod(6:6)='1'
    
	faultandslickenline_loop: do
		write (*,"(3x,A)", ADVANCE="no") 'Strike: '
		read (*,*) analysissettings%statistics%centrfaultdat%fltplane%strike_rhr        

		write (*,"(3x,A)", ADVANCE="no") 'Dip angle: '
		read (*,*) analysissettings%statistics%centrfaultdat%fltplane%dipangle  
    
		write (*,"(3x,A)", ADVANCE="no") 'Slickenline trend and plunge: '
		read (*,*) analysissettings%statistics%centrfaultdat%slickenline
	
		if(isaxisonplane(analysissettings%statistics%centrfaultdat%fltplane &
    			,analysissettings%statistics%centrfaultdat%slickenline)) then
       		exit faultandslickenline_loop
    	else
      		write(*,*) "slickenline doesn't lie on plane"
    	endif
    end do faultandslickenline_loop

elseif (analysissettings%statistics%centrflt_orientcod == 2) then  !  2  -  P-axis trend & plunge T-axis trend & plunge

	analysissettings%statistics%centrflt_generalcod(6:6)='2'
    
	TandPaxes_loop: do
		write (*,"(3x,A)", ADVANCE="no") 'T central axis trend and plunge: '
		read (*,*) analysissettings%statistics%centrfaultdat%focmech%axis_a
        
		write (*,"(3x,A)", ADVANCE="no") 'P central axis trend and plunge: '
		read (*,*) analysissettings%statistics%centrfaultdat%focmech%axis_b

		vector_a = pole2cartesian(analysissettings%statistics%centrfaultdat%focmech%axis_a)
		vector_b = pole2cartesian(analysissettings%statistics%centrfaultdat%focmech%axis_b) 

        angle_between_t0p0 = r2d*axes_angle_rad(vector_a,vector_b)
        if (angle_between_t0p0 < 89) then
			write(*,*) 'The two axes are not perpendicular. Re-enter values'
		else
			exit TandPaxes_loop
        end if
        
	end do TandPaxes_loop

endif


    
end subroutine definecentralfaultdatum

!------------------------



!------------------------
! define if specific statistic

logical function defineoutput_statistic() result(analyse)

character :: dummy_input

write(*,"(5x,A)",advance="no") 'calculate (y/n)?: '
read(*,*) dummy_input  

do
	if (dummy_input == 'y' .or. dummy_input == 'Y') then
		analyse = .true.
        exit
	elseif (dummy_input == 'n' .or. dummy_input == 'N') then
		analyse = .false. 
        exit
    endif 
end do

end function defineoutput_statistic

!------------------------  


!-------------------------
! program reads data from input file (according to the different input formats)
!  and writes them to a temp file

type(fault_datum)  function readinputdata(unitnum, dataformatcod) result(fault_rec1)

integer, intent(in) :: unitnum
character (len=6), intent(in)  :: dataformatcod

   
select case (dataformatcod)

! only orientation
  case("000010")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%rake_aki
  
  case("000011")  
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%slickenline%trend,fault_rec1%slickenline%plunge 
  
  case("000012")
  read(unitnum,*) fault_rec1%id  &     
  ,fault_rec1%focmech%axis_a%trend,fault_rec1%focmech%axis_a%plunge  &
  ,fault_rec1%focmech%axis_b%trend,fault_rec1%focmech%axis_b%plunge 
    
! time missing
  case("100010")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3)  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%rake_aki

  case("100011")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3)  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%slickenline%trend,fault_rec1%slickenline%plunge 

  case("100012")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3)  &      
  ,fault_rec1%focmech%axis_a%trend,fault_rec1%focmech%axis_a%plunge  &
  ,fault_rec1%focmech%axis_b%trend,fault_rec1%focmech%axis_b%plunge 

  case("110010")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3)  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%rake_aki
  
  case("110011")  
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%sphercoord%lat,fault_rec1%sphercoord%lon,fault_rec1%sphercoord%depth  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%slickenline%trend,fault_rec1%slickenline%plunge 

  case("110012")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%sphercoord%lat,fault_rec1%sphercoord%lon,fault_rec1%sphercoord%depth  &
  ,fault_rec1%focmech%axis_a%trend,fault_rec1%focmech%axis_a%plunge  &
  ,fault_rec1%focmech%axis_b%trend,fault_rec1%focmech%axis_b%plunge 
    
! space missing
  case("001010")
  read(unitnum,*) fault_rec1%id,fault_rec1%t  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%rake_aki
  
  case("001011")
  read(unitnum,*) fault_rec1%id,fault_rec1%t  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%slickenline%trend,fault_rec1%slickenline%plunge 
  
  case("001012")
  read(unitnum,*) fault_rec1%id,fault_rec1%t  & 
  ,fault_rec1%focmech%axis_a%trend,fault_rec1%focmech%axis_a%plunge  &
  ,fault_rec1%focmech%axis_b%trend,fault_rec1%focmech%axis_b%plunge 

! space, time present
  case("101010")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3),fault_rec1%t  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%rake_aki

  case("101011")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3),fault_rec1%t  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%slickenline%trend,fault_rec1%slickenline%plunge 

  case("101012")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%spatloc(1),fault_rec1%spatloc(2),fault_rec1%spatloc(3),fault_rec1%t  &
  ,fault_rec1%focmech%axis_a%trend,fault_rec1%focmech%axis_a%plunge  &
  ,fault_rec1%focmech%axis_b%trend,fault_rec1%focmech%axis_b%plunge 

  case("111010")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%sphercoord%lat,fault_rec1%sphercoord%lon,fault_rec1%sphercoord%depth,fault_rec1%t  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%rake_aki
  
  case("111011")  
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%sphercoord%lat,fault_rec1%sphercoord%lon,fault_rec1%sphercoord%depth,fault_rec1%t  & 
  ,fault_rec1%fltplane%strike_rhr,fault_rec1%fltplane%dipangle  &
  ,fault_rec1%slickenline%trend,fault_rec1%slickenline%plunge 

  case("111012")
  read(unitnum,*) fault_rec1%id  &
  ,fault_rec1%sphercoord%lat,fault_rec1%sphercoord%lon,fault_rec1%sphercoord%depth,fault_rec1%t  &
  ,fault_rec1%focmech%axis_a%trend,fault_rec1%focmech%axis_a%plunge  &
  ,fault_rec1%focmech%axis_b%trend,fault_rec1%focmech%axis_b%plunge 

  
end select


end function readinputdata

!------------------------


!------------------------
! process the fault record, for the successive writing into a temporary file
type(fault_proc) function processinputdata(fault_rec1,dataformatcod,analysissettings) result(fault_proc1)


character (len=6), intent(in)  :: dataformatcod
type(analysistype), intent(in) :: analysissettings
type(fault_datum) :: fault_rec1


! common processings

! 1 - writes id
fault_proc1%id = fault_rec1%id


! 2 - calculates slickenline vector and rake
if (dataformatcod(6:6) == "1") then  
  	fault_rec1%slickenl_vect = pole2cartesian(fault_rec1%slickenline)
	call slickenline2rake(fault_rec1)
endif


! spatial separation statistics
if (analysissettings%spatialsep) then

	! if spherical coordinates, calculate cartesian ones
	! convert spherical coordinates to cartesian ones
	if (dataformatcod(2:2) == "1") then
		fault_proc1%spatloc = pole2earthcartesian(fault_rec1%sphercoord)
	else
		fault_proc1%spatloc = fault_rec1%spatloc   
	endif

endif


! time lag statistics
if (analysissettings%timelag) then
  
	fault_proc1%time = fault_rec1%t   

endif


! coherence statistics and/or focal mechanism rotation statistics
if (analysissettings%coherence .or. analysissettings%focmechrot .or. analysissettings%kinematicaxesrot) then

	if (dataformatcod(6:6) == "0" .or. dataformatcod(6:6) == "1") then  

	  fault_proc1%fochmecv = fochalmechanismvectorcalc0(fault_rec1%fltplane%strike_rhr  &
                               ,fault_rec1%fltplane%dipangle,fault_rec1%rake_aki)

	  fault_proc1%mom_tens = incremstraintensorcalc0(fault_rec1%fltplane%strike_rhr  &
                               ,fault_rec1%fltplane%dipangle,fault_rec1%rake_aki)
                               
	elseif (dataformatcod(6:6) == "2") then

      fault_proc1%fochmecv%X = pole2cartesian(fault_rec1%focmech%axis_a)
	  fault_proc1%fochmecv%Y = pole2cartesian(fault_rec1%focmech%axis_b)
		
	  fault_proc1%fochmecv%Z = vector_vectprod(fault_proc1%fochmecv%X,fault_proc1%fochmecv%Y)
	  fault_proc1%fochmecv%Y = vector_vectprod(fault_proc1%fochmecv%Z,fault_proc1%fochmecv%X)      
                               
	  fault_proc1%mom_tens = incremstraintensorcalc1(fault_proc1%fochmecv%X,fault_proc1%fochmecv%Y)

    endif
         
endif


! fault sub-element rotation statistics
if (analysissettings%faultelementsrot) then

	fault_proc1%faultnorm_vect = plane_normal(fault_rec1%fltplane)
    
	if (dataformatcod(6:6) == "0") then  
		fault_proc1%slickenl_vect = faultrake2slick_vector(fault_rec1%fltplane%strike_rhr &
    						,fault_rec1%fltplane%dipangle,fault_rec1%rake_aki)

	elseif (dataformatcod(6:6) == "1") then  
		fault_proc1%slickenl_vect = pole2cartesian(fault_rec1%slickenline)

	endif
      
endif


end function processinputdata

!------------------------


!------------------------

subroutine writeheaderoutputdatafile(unitnum, analysissettings)


type(analysistype), intent(in) :: analysissettings
integer, intent(in) :: unitnum 

character (len=500) :: headerstring


write(headerstring, "(A)") 'RecA_id,RecB_id'


! spatial separation statistics
if (analysissettings%spatialsep) then
	write(headerstring, "(2A)") trim(headerstring),',RecDist'
endif


! time lag statistics
if (analysissettings%timelag) then
	write(headerstring, "(2A)") trim(headerstring),',TimeLag' 
endif


! coherence statistics 
if (analysissettings%coherence) then
	write(headerstring, "(2A)") trim(headerstring),',Coherence'          
endif


! focal mechanism rotation statistics
if (analysissettings%focmechrot) then  
  write(headerstring, "(2A)") trim(headerstring),',FMSol1_Tr,FMSol1_Pl,FMSol1_An'
  write(headerstring, "(2A)") trim(headerstring),',FMSol2_Tr,FMSol2_Pl,FMSol2_An'
  write(headerstring, "(2A)") trim(headerstring),',FMSol3_Tr,FMSol3_Pl,FMSol3_An' 
  write(headerstring, "(2A)") trim(headerstring),',FMSol4_Tr,FMSol4_Pl,FMSol4_An'   
endif


! T-P-B axes rotations
if (analysissettings%kinematicaxesrot) then

	write(headerstring, "(2A)") trim(headerstring),',TSol_Tr,TSol_Pl,TSol_An'  
	write(headerstring, "(2A)") trim(headerstring),',PSol_Tr,PSol_Pl,PSol_An'  
	write(headerstring, "(2A)") trim(headerstring),',BSol_Tr,BSol_Pl,BSol_An' 

endif


! fault sub-element rotation statistics
if ( analysissettings%faultelementsrot)  then
    
	write(headerstring, "(2A)") trim(headerstring),',FPSol_Tr,FPSol_Pl,FPSol_An'  
  	write(headerstring, "(2A)") trim(headerstring),',SlSol_Tr,SlSol_Pl,SlSol_An'  

endif


write(unitnum, "(A)") trim(headerstring)


end subroutine writeheaderoutputdatafile

!------------------------


!------------------------
! write the results in the output file

subroutine write_results(outputfileparams,analysissettings,faultproc_1,faultproc_2,faultpair_res1)	


	type(file_params), intent(in) :: outputfileparams
	type(analysistype), intent(in) :: analysissettings
	type(fault_proc), intent(in) :: faultproc_1,faultproc_2
	type(analysis_results), intent(in) :: faultpair_res1

	character (len=500) :: resultstring
	integer (kind=i4b) :: i

    

	write(resultstring, "(i0,A,i0)") faultproc_1%id,',',faultproc_2%id


	! spatial separation statistics
	if (analysissettings%spatialsep) then
		write(resultstring, "(2A,f10.2)") trim(resultstring),',',faultpair_res1%spat_sep
	endif


	! time lag statistics
	if (analysissettings%timelag) then
		write(resultstring, "(2A,f10.3)") trim(resultstring),',',faultpair_res1%time_lag 
	endif


	! coherence statistics 
	if (analysissettings%coherence) then
		write(resultstring, "(2A,f10.3)") trim(resultstring),',',faultpair_res1%coherence          
	endif


	! focal mechanism rotation statistics
	if (analysissettings%focmechrot) then
			
		if (abs(faultpair_res1%focmech_rots(1)%rot_angle) < 1.0) then
 			do i = 1,4           
				write(resultstring, "(2A)") trim(resultstring),',,,0.0'
         	end do               
        else
 			do i = 1,4           
			 write(resultstring, "(A,3(A,f10.2))") trim(resultstring),','   &
              ,faultpair_res1%focmech_rots(i)%rot_axis%trend,',',faultpair_res1%focmech_rots(i)%rot_axis%plunge,','  &
              ,faultpair_res1%focmech_rots(i)%rot_angle
         	end do               
        endif    
    
	endif

	
	! T-P-B axes rotations
	if (analysissettings%kinematicaxesrot) then
		
		! T axis
	  if (abs(faultpair_res1%Taxis_rot%rot_angle) < 1.0) then
		write(resultstring, "(2A,f10.2)") trim(resultstring),',,,',faultpair_res1%Taxis_rot%rot_angle
      else      
  	    write(resultstring, "(A,3(A,f10.2))") trim(resultstring),',',faultpair_res1%Taxis_rot%rot_axis%trend   &
        ,',',faultpair_res1%Taxis_rot%rot_axis%plunge,',',faultpair_res1%Taxis_rot%rot_angle 
      endif 

		! P axis
	  if (abs(faultpair_res1%Paxis_rot%rot_angle) < 1.0) then
		write(resultstring, "(2A,f10.2)") trim(resultstring),',,,',faultpair_res1%Paxis_rot%rot_angle
      else 
  	  write(resultstring, "(A,3(A,f10.2))") trim(resultstring),',',faultpair_res1%Paxis_rot%rot_axis%trend   &
        ,',',faultpair_res1%Paxis_rot%rot_axis%plunge,',',faultpair_res1%Paxis_rot%rot_angle 
      endif 

		! B axis
	  if (abs(faultpair_res1%Baxis_rot%rot_angle) < 1.0) then
		write(resultstring, "(2A,f10.2)") trim(resultstring),',,,',faultpair_res1%Baxis_rot%rot_angle
      else              
  	  write(resultstring, "(A,3(A,f10.2))") trim(resultstring),',',faultpair_res1%Baxis_rot%rot_axis%trend   &
        ,',',faultpair_res1%Baxis_rot%rot_axis%plunge,',',faultpair_res1%Baxis_rot%rot_angle  
      endif
      
    endif


	! fault sub-element rotation statistics
     if (analysissettings%faultelementsrot) then

	  if (abs(faultpair_res1%faultpl_rot%rot_angle) < 1.0) then
		write(resultstring, "(2A,f10.2)") trim(resultstring),',,,',faultpair_res1%faultpl_rot%rot_angle
      else  
		write(resultstring, "(A,3(A,f10.2))") trim(resultstring),',',faultpair_res1%faultpl_rot%rot_axis%trend   &
         	,',',faultpair_res1%faultpl_rot%rot_axis%plunge,',',faultpair_res1%faultpl_rot%rot_angle  
      endif   

	  if (abs(faultpair_res1%slickenl_rot%rot_angle) < 1.0) then
		write(resultstring, "(2A,f10.2)") trim(resultstring),',,,',faultpair_res1%slickenl_rot%rot_angle
      else           
		write(resultstring, "(A,3(A,f10.2))") trim(resultstring),',',faultpair_res1%slickenl_rot%rot_axis%trend   &
         	,',',faultpair_res1%slickenl_rot%rot_axis%plunge,',',faultpair_res1%slickenl_rot%rot_angle
      endif
            
	endif


	write(outputfileparams%unitnum, "(A)") resultstring



end subroutine write_results

!------------------------


!------------------------
! write metadata 

subroutine write_metadata(time_initial, time_final, metadatafileparams  &
		,inputfileparams,inputdataformat,outputfileparams,analysissettings)


integer (kind=i4b), intent(in) :: time_initial(8),time_final(8)
type(file_params), intent(in) :: metadatafileparams, inputfileparams, outputfileparams
type(dataformat), intent(in):: inputdataformat
type(analysistype), intent(in) :: analysissettings


write(metadatafileparams%unitnum,"(A)") 'Kinematic correlation metadata'
write(metadatafileparams%unitnum,"(A)") 'Output of FaultCorrelation_01.f95 program'

write (metadatafileparams%unitnum,"(A,1x,i0,A,i0,A,i0,3x,3(i0,A))") 'Analysis starts at:',time_initial(1),'-'    &
 ,time_initial(2),'-',time_initial(3),time_initial(5),':',time_initial(6),':',time_initial(7)
write (metadatafileparams%unitnum,"(A,1x,i0,A,i0,A,i0,3x,3(i0,A))") 'Analysis ends at:',time_final(1),'-'    &
 ,time_final(2),'-',time_final(3),time_final(5),':',time_final(6),':',time_final(7)

write (metadatafileparams%unitnum,*)

write(metadatafileparams%unitnum,"(A,A)") 'Input file:    ',inputfileparams%name
write(metadatafileparams%unitnum,"(A,A)") 'Output file:   ',outputfileparams%name

write(metadatafileparams%unitnum,*)

write(metadatafileparams%unitnum,"(A,i0)") 'Number of records:  ',inputfileparams%recsnumber

write(metadatafileparams%unitnum,*)

write(metadatafileparams%unitnum,"(2A)") 'Types of statistical analysis: '

if (analysissettings%focmechrot) then
	write(metadatafileparams%unitnum,"(4x,A)") 'focal mechanism rotation:  true'
else
	write(metadatafileparams%unitnum,"(4x,A)") 'focal mechanism rotation:  false'
endif

if (analysissettings%kinematicaxesrot) then
	write(metadatafileparams%unitnum,"(4x,A)") 'single P-T-B axes rotation:  true'
else
	write(metadatafileparams%unitnum,"(4x,A)") 'single P-T-B axes rotation:  false'
endif

if (analysissettings%faultelementsrot) then
	write(metadatafileparams%unitnum,"(4x,A)") 'fault plane and slickenline axes rotation:  true'
else
	write(metadatafileparams%unitnum,"(4x,A)") 'fault plane and slickenline axes rotation:  false'
endif

if (analysissettings%coherence) then
	write(metadatafileparams%unitnum,"(4x,A)") 'coherence:  true'
else
	write(metadatafileparams%unitnum,"(4x,A)") 'coherence:  false'
endif

if (analysissettings%spatialsep) then
	write(metadatafileparams%unitnum,"(4x,A)") 'spatial separation:  true'
else
	write(metadatafileparams%unitnum,"(4x,A)") 'spatial separation:  false'
endif

if (analysissettings%timelag) then
	write(metadatafileparams%unitnum,"(4x,A)") 'time lag:  true'
else
	write(metadatafileparams%unitnum,"(4x,A)") 'time lag:  false'
endif



write(metadatafileparams%unitnum,"(4x,A)") 'statistics type:'


if (analysissettings%statistics%cod == 0) then

	write(metadatafileparams%unitnum,"(6x,A)") 'pairwise'

elseif (analysissettings%statistics%cod == 1) then
  
	write(metadatafileparams%unitnum,"(6x,A)") 'central'


	write(metadatafileparams%unitnum,"(6x,A)") 'parameters:'

    
	! the spatial coordinates of the central datum
	if (analysissettings%spatialsep) then
    
    	if (inputdataformat%spatial%cod == 0) then
    		write(metadatafileparams%unitnum,"(6x,A,3(f10.5,2x))") 'location x, y, z: ' &
            							, analysissettings%statistics%centrfaultdat%spatloc
    	elseif (inputdataformat%spatial%cod == 1) then
			write(metadatafileparams%unitnum,"(6x,A,3(f10.5,2x))") 'location lat-lon-depth: '  &
            							, analysissettings%statistics%centrfaultdat%sphercoord
    	endif

    endif


	! the time coordinates of the central datum
	if (analysissettings%timelag) then

		write(metadatafileparams%unitnum,"(6x,A,f10.5))") 'time: ', analysissettings%statistics%centrfaultdat%t

	endif


	! the orientation of the central datum	
	write(metadatafileparams%unitnum,"(6x,A)") 'orientation: ' 
	centralstat_orientation: select case (analysissettings%statistics%centrflt_orientcod)
    
    	case(0) !0 - strike dip rake
          
			write(metadatafileparams%unitnum,"(6x,A,f6.1,3x,A,f5.1,3x,A,f6.1)")    &
       		'fault strike: ',analysissettings%statistics%centrfaultdat%fltplane%strike_rhr  &
       		,' ; fault dip: ',analysissettings%statistics%centrfaultdat%fltplane%dipangle   &
       		,' ; slickenline rake: ',analysissettings%statistics%centrfaultdat%rake_aki
    
   
		case(1) ! 1 - strike dip slick. trend & plunge
  
			write(metadatafileparams%unitnum,"(6x,A,f6.1,3x,A,f5.1,2(3x,A,f6.1))")    &
      		'fault strike: ',analysissettings%statistics%centrfaultdat%fltplane%strike_rhr  &
     		,' ; fault dip: ',analysissettings%statistics%centrfaultdat%fltplane%dipangle   &
     		,' ; slickenline trend: ',analysissettings%statistics%centrfaultdat%slickenline%trend  &
     		,' ; slickenline plunge: ',analysissettings%statistics%centrfaultdat%slickenline%plunge
          
		case(2) ! 2 - P-axis trend & plunge T-axis trend & plunge
          
			write(metadatafileparams%unitnum,"(6x,2(2A,f6.1,3x,A,f5.1))")    &
       		'T axis trend: ',analysissettings%statistics%centrfaultdat%focmech%axis_a%trend    &
       		,' ; T axis plunge: ',analysissettings%statistics%centrfaultdat%focmech%axis_a%plunge  &
       		,' ; P axis trend: ',analysissettings%statistics%centrfaultdat%focmech%axis_b%trend   &
       		,' ; P axis plunge: ',analysissettings%statistics%centrfaultdat%focmech%axis_b%plunge
    
	end select centralstat_orientation

end if

            
end subroutine write_metadata

!------------------------



end module inputoutput_defs

!------------------------------------------------ 




!------------------------------------------------ 

module file_management

use inputoutput_defs

contains



!-------------------------
! user chooses existing file name and program opens it

function sequentialfile_choice(unitnum) result(filename)

	integer, intent(in) :: unitnum
	integer :: ios
	character (len= 255) :: filename  
 	logical :: exists    


	do
		write (*,"(A)", ADVANCE="no") 'Enter filename: '
		read (*,*) filename

		filename = trim(filename)
        
		inquire(file=filename,exist=exists)
		if (.not.exists) then
			write(*,*) 'File not found.'
			cycle        
		end if
    
		! open sequential-access file
		open (unit=unitnum,file=filename,status='old' &
  		, access='sequential', form='formatted', iostat=ios)

		if (ios /= 0) then
  			write (*,"(A)") 'Error with file opening. Program stops'            
  			stop
		end if

        exit

	end do
 

end function sequentialfile_choice

!-------------------------


!-------------------------
! defines new file 
! containing the results and the metadata

subroutine newsequentialfile(newfileparams)

	type (file_params), intent(inout) :: newfileparams
	integer :: ios
	character (len= 50) :: filename  
 	logical :: exists 


        
	do
		write (*,"(A)", ADVANCE="no") 'Enter filename: '
		read (*,*) filename

        newfileparams%name = trim(filename)
        
		inquire(file=newfileparams%name,exist=exists)
		if (exists) then
			write(*,*) 'File already existing. Change name'
			cycle        
		end if

		! creates sequential-access output file
		open (unit=newfileparams%unitnum,file=newfileparams%name,status='new' &
  		, access='sequential', form='formatted', iostat=ios)

		if (ios /= 0) then
  			write (*,"(A,/,A)") 'Error with file creation.','Change name'
  			cycle
		end if

        exit

	end do

    write (*,*)
    

end subroutine newsequentialfile

!-------------------------


!-------------------------
! program creates direct-access file
subroutine directfilecreation(directfileparams)

type(file_params) :: directfileparams
integer  :: ios
character (len=50) :: status_name

if (directfileparams%statuscod == 0) then
  status_name = 'replace'
endif

open (unit=directfileparams%unitnum,file=directfileparams%name,status=status_name &
 , access='direct', form='unformatted',recl=directfileparams%reclength, iostat=ios)

if (ios /= 0) then
  write (*,"(A)") 'Error with file creation. Program will stop'
  read(*,*)
  stop
end if


end subroutine directfilecreation

!-------------------------


!------------------------
! calculates the number of records in sequential file    
    
integer (kind=i4b) function rowsnumbercount(fileunitnumber)

    integer, intent(in) :: fileunitnumber
    integer (kind=i4b) :: j 
    integer :: ios
    character :: dummy_character
  
      
	j=0
    
	do
		j=j+1
		read(fileunitnumber,"(A)",iostat=ios) dummy_character   
        if (ios /= 0) exit
	end do 
        
	rowsnumbercount = j-1

    
    rewind(unit=fileunitnumber)

       
end function rowsnumbercount

!------------------------


end module file_management


!------------------------------------------------ 



!------------------------------------------------ 

module statistical_processing

use file_management

contains

!------------------------
! calculates the similarity statistics for the pairwise case

subroutine calculatepairwisestatistics(tempfileparams  &
  								,analysissettings, outputfileparams) 


type(analysistype), intent(in) :: analysissettings
type(file_params), intent(in) :: tempfileparams, outputfileparams

integer  :: ios
integer (kind=i4b) :: k,l

type(fault_proc) :: faultproc_1, faultproc_2

type(analysis_results) :: faultpair_res1

do k = 1, tempfileparams%recsnumber-1
! k: cycle from 1 to (n-1)
	
	! read from TEMP file: rec k code and ptb axes coordinates
	read(tempfileparams%unitnum,rec=k,iostat=ios) faultproc_1

	do l = k+1, tempfileparams%recsnumber
    ! l: cycle from (k+1) to n    
       
		! read from TEMP file: rec l code and ptb axes coordinates
		read(tempfileparams%unitnum,rec=l,iostat=ios) faultproc_2

		write(*,*) 'analysing records ', faultproc_1%id, ' ,',faultproc_2%id
        
		faultpair_res1 = calculatefaultpairstatistics(analysissettings,faultproc_1,faultproc_2)
              
		! write the results in the output file
		call write_results(outputfileparams,analysissettings,faultproc_1,faultproc_2,faultpair_res1)

	! end cycle l
	end do
    
! end cycle k
end do 

end subroutine calculatepairwisestatistics


!------------------------


!------------------------
! calculates the similarity statistics for the central datum case

subroutine calculatecentralstatistics(tempfileparams,analysissettings, outputfileparams)


type(analysistype), intent(in) :: analysissettings
type(file_params), intent(in) :: tempfileparams, outputfileparams

integer :: ios
integer (kind=i4b) :: k


type(fault_proc) :: faultproc_0, faultproc_1
type(analysis_results) :: faultpair_res1


faultproc_0 = processinputdata(analysissettings%statistics%centrfaultdat  &
		,analysissettings%statistics%centrflt_generalcod,analysissettings)


do k = 1, tempfileparams%recsnumber
! k: cycle from 1 to n
	
	! read from TEMP file: rec k code and ptb axes coordinates
	read(tempfileparams%unitnum,rec=k,iostat=ios) faultproc_1

	write(*,*) 'analysing records ', faultproc_0%id, ' ,',faultproc_1%id
        
	faultpair_res1 = calculatefaultpairstatistics(analysissettings,faultproc_0,faultproc_1)
              
	! write the results in the output file
	call write_results(outputfileparams,analysissettings,faultproc_0,faultproc_1,faultpair_res1)		

end do
    

end subroutine calculatecentralstatistics

!------------------------



!------------------------
! calculates the statistics for a given fault pair
type (analysis_results) function calculatefaultpairstatistics(analysissettings,faultproc_1,faultproc_2) result(faultpair_res1)


type(analysistype), intent(in) :: analysissettings
type(fault_proc), intent(in) :: faultproc_1, faultproc_2

   
 if (analysissettings%focmechrot) then ! focal mechanism rotation statistics
  ! given the two different ptb sets, calculates rotation quaternion 
  faultpair_res1%focmech_rots = triadvectors_rotsolution(faultproc_1%fochmecv,faultproc_2%fochmecv)
 endif

            
 if (analysissettings%kinematicaxesrot) then ! T-P-B axes rotations
  faultpair_res1%Taxis_rot = axis_rotsolution(faultproc_1%fochmecv%X,faultproc_2%fochmecv%X) ! T axis
  faultpair_res1%Paxis_rot = axis_rotsolution(faultproc_1%fochmecv%Y,faultproc_2%fochmecv%Y) ! P axis  
  faultpair_res1%Baxis_rot = axis_rotsolution(faultproc_1%fochmecv%Z,faultproc_2%fochmecv%Z) ! B axis
 endif

  
 if (analysissettings%faultelementsrot) then ! fault sub-element rotation statistics
  faultpair_res1%faultpl_rot = axis_rotsolution(faultproc_1%faultnorm_vect,faultproc_2%faultnorm_vect) !fault plane
  faultpair_res1%slickenl_rot = vector_rotsolution(faultproc_1%slickenl_vect,faultproc_2%slickenl_vect) ! slickenline
 endif

             
 if (analysissettings%coherence) then ! coherence statistics
  faultpair_res1%coherence = coherence_solution(faultproc_1%mom_tens,faultproc_2%mom_tens)
 endif

  
 if (analysissettings%spatialsep) then ! spatial separation statistics
  faultpair_res1%spat_sep = euclidean_distance(faultproc_1%spatloc,faultproc_2%spatloc)
 endif

  
 if (analysissettings%timelag) then ! time lag statistics
  faultpair_res1%time_lag = faultproc_2%time - faultproc_1%time
 endif 


end function calculatefaultpairstatistics

!------------------------


!------------------------
! calculates rotation solution for a vector pair
type (rotation) function vector_rotsolution(vector1,vector2) result(rotsolution)

! QUATr_axesrotation
! created: 2005-03-26
!  original name: "r_axesrotation"
!  2007-09-16: modified to handle two parallel vectors

implicit none

type(vector), intent(in) :: vector1,vector2

type(vector) :: vector3

rotsolution%rot_angle = r2d*vector_angle_rad(vector1,vector2)

! case of sub-parallel vectors
if(rotsolution%rot_angle < 1.0) then 
  	rotsolution%rot_axis%trend = 0.0
  	rotsolution%rot_axis%plunge = 0.0
else ! defined rotation pole
	vector3 = vector_vectprod(vector1,vector2)
	vector3 = vector_normalization(vector3)
    
    rotsolution%rot_axis = cartesian2pole(vector3)
    
    if(rotsolution%rot_axis%plunge < 0.0) then
    	rotsolution%rot_axis = axis2downaxis(rotsolution%rot_axis)
		rotsolution%rot_angle = - rotsolution%rot_angle
	endif

    
endif

end function vector_rotsolution

!------------------------



!------------------------
! calculates rotation solution for a vector pair (considered as axes)
type (rotation) function axis_rotsolution(vector1,vector2) result(rotsolution)

! QUATr_axesrotation
! created: 2005-03-26
!  original name: "r_axesrotation"
!  2007-09-16: modified to handle two parallel vectors

implicit none

type(vector), intent(in) :: vector1,vector2

type(vector) :: vector3

rotsolution%rot_angle = r2d*vector_angle_rad(vector1,vector2)

! case of sub-parallel vectors
if(rotsolution%rot_angle < 1.0) then 
  	rotsolution%rot_axis%trend = 0.0
  	rotsolution%rot_axis%plunge = 0.0
else ! defined rotation pole
	vector3 = vector_vectprod(vector1,vector2)
	vector3 = vector_normalization(vector3)
    
    rotsolution%rot_axis = cartesian2pole(vector3)
    
    if (rotsolution%rot_angle > 90.0) then
		rotsolution%rot_axis = antipole(rotsolution%rot_axis)
		rotsolution%rot_angle = 180.0 - rotsolution%rot_angle
    endif

    if(rotsolution%rot_axis%plunge < 0.0) then
    	rotsolution%rot_axis = axis2downaxis(rotsolution%rot_axis)
		rotsolution%rot_angle = - rotsolution%rot_angle
	endif

    
endif

end function axis_rotsolution

!------------------------



end module statistical_processing


!------------------------------------------------ 





!------------------------------------------------ 


program FaultCorrelation


use statistical_processing

implicit none

type(dataformat) :: inputdataformat
type(analysistype) :: analysissettings
type(file_params) :: inputfileparams, tempfileparams, outputfileparams, metadatafileparams
   
type(fault_datum) :: fault_rec1
type(fault_proc) :: fault_proc1

integer (kind=i4b) :: time_initial(8),time_final(8)


integer (kind=i4b) :: dummy_i4b, i
integer :: ios

pi = dacos(0.0d0)*2.0d0 ! pi radians
r2d = 180.d0/pi			! conversion from radians to degrees
d2r = pi/180.d0			! conversion from degrees to radians

! input file
inputfileparams%unitnum = 17

! temporary file
tempfileparams%name = 'temp.dat'
tempfileparams%statuscod = 0 ! 'replace'
tempfileparams%unitnum = 18
inquire(iolength=dummy_i4b) faultproc1
tempfileparams%reclength  = dummy_i4b

! output file
outputfileparams%unitnum = 19


! metadata file
metadatafileparams%unitnum = 20


!------------------------------------------------------


! gets the initial time_date of analysis 
call DATE_AND_TIME(VALUES=time_initial)

! print initial message
call initial_message()

!-- user definitions of input information/format

write(*,*)
write(*,"(10x,A)") '-------------        Input data format        ------------'
write(*,*)

! user defines spatial information
call defineinput_spatialinfo(inputdataformat%spatial)

! user defines time information
call defineinput_timeinfo(inputdataformat%time)

! user defines orientation information
call defineinput_orientationinfo(inputdataformat%orientation)

! program codes the input file format
! based on previous user choices 
call inputfile_coding(inputdataformat)

!-- user definitions of output information/format

write(*,*)
write(*,*)
write(*,"(10x,A)") '-------------      Data analysis settings     -------------'
write(*,*)

! initialize the settings for analysis
analysissettings%focmechrot = .false.
analysissettings%kinematicaxesrot = .false.
analysissettings%faultelementsrot = .false.
analysissettings%coherence = .false.
analysissettings%spatialsep = .false.
analysissettings%timelag = .false.


! user defines spatial separation statistics
if (inputdataformat%spatial%present) then
  	write(*,*) 
	write(*,"(A)") 'Spatial separation between events' 
	analysissettings%spatialsep = defineoutput_statistic()
endif


! user defines time lag statistics
if (inputdataformat%time%present) then
  	write(*,*) 
  	write(*,"(A)") 'Time lag between events'  
	analysissettings%timelag = defineoutput_statistic()
endif


! user defines coherence statistics
write(*,*) 
write(*,"(A)") 'Coherence statistics'
analysissettings%coherence = defineoutput_statistic()


! user defines focal mechanism rotation statistics
write(*,*) 
write(*,"(A)") 'Rotation statistics for focal mechanism'
analysissettings%focmechrot = defineoutput_statistic()

  
! user defines single P-T-B axes rotation statistics
write(*,*) 
write(*,"(A)") 'Rotation statistics for each kinematic axis'
analysissettings%kinematicaxesrot = defineoutput_statistic()
 
  
! user defines fault plane and slickenline axes rotation statistics
if (inputdataformat%orientation%cod /= 2) then
	write(*,*) 
	write(*,"(A)") 'Rotation statistics for fault plane and slickenline'
	write(*,"(4x,A)") '(if you want to perform central statistics, note that you will not'
    write(*,"(4x,A)") 'be able to enter central datum orientation as T-P axes)'  
	analysissettings%faultelementsrot = defineoutput_statistic()
endif

 
! user defines pairwise or central statistics
call defineoutput_statisticsinfo(inputdataformat, analysissettings)


write(*,*)
write(*,*)
write(*,"(10x,A)") '-------------          Input/output           -------------'
write(*,*)


! user defines input file name and program opens it
write(*,"(5x,A)") 'Input file'
inputfileparams%name = sequentialfile_choice(inputfileparams%unitnum)

write(*,"(2x,A)", ADVANCE="no") 'Enter number of header rows in input file (default: 1): '
read(*,*) inputfileparams%headerrownumber

! user defines output file and program opens it
write(*,*) 
write(*,"(5x,A)") 'Output data file'
call newsequentialfile(outputfileparams)
write(*,*) 
write(*,"(5x,A)") 'Metadata file'
call newsequentialfile(metadatafileparams)


! - end of user interaction -


! program counts the number of records in input file
inputfileparams%recsnumber = rowsnumbercount(inputfileparams%unitnum) - inputfileparams%headerrownumber
tempfileparams%recsnumber = inputfileparams%recsnumber


! program creates the temporary file that stores the data to be analyzed
call directfilecreation(tempfileparams)


! program reads data from input file (according to the different input formats)
!  and writes them to the temp file

! skip header
if (inputfileparams%headerrownumber > 0) then
  do i=1,inputfileparams%headerrownumber
	read(inputfileparams%unitnum,*)
  end do
endif

! read data, process them and write into temporary file
do i = 1, inputfileparams%recsnumber
	fault_rec1 = readinputdata(inputfileparams%unitnum, inputdataformat%cod)
    fault_proc1 = processinputdata(fault_rec1,inputdataformat%cod,analysissettings)
    write(tempfileparams%unitnum,rec=i,iostat=ios) fault_proc1
end do


! program write header in data output file
call writeheaderoutputdatafile(outputfileparams%unitnum, analysissettings)


! program calculates and output results
statistic_case: select case(analysissettings%statistics%cod)
  case(0) 
   call calculatepairwisestatistics(tempfileparams  &
  								,analysissettings, outputfileparams)
  case(1) 
   call calculatecentralstatistics(tempfileparams  &
  								,analysissettings, outputfileparams)
end select statistic_case


! program deletes the temporary file
close(unit=tempfileparams%unitnum,status='delete')        

call DATE_AND_TIME(VALUES=time_final)

call write_metadata(time_initial, time_final, metadatafileparams,inputfileparams,inputdataformat, outputfileparams,analysissettings)

write(*,*)
write(*,"(2x,A)", ADVANCE="no") 'Analysis completed. Press any key to close the window'
read(*,*) 

 
end program FaultCorrelation