module faults

    implicit none

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

    ! pairwise or central statistics
    type :: statistics_type
        integer (i1b) :: cod  ! 0 - pairwise; 1 - central statistics
        integer (i1b) :: centrflt_orientcod  ! same codes as orientation information
        character (len=6) :: centrflt_generalcod
        type (fault_datum)  ::  centrfaultdat
    end type statistics_type

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


end module faults