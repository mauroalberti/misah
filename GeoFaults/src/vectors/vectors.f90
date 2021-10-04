module vectors

    implicit none

    type :: vector
        sequence
        real (kind=r8b)	:: x,y,z
    end type vector

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



end module vectors