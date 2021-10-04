module orientations

    use kinds_mod

    implicit none

    type :: axis
        sequence
        real (kind=r8b)	:: trend, plunge
    end type axis

    type :: plane
        sequence
        real (kind=r8b)	:: strike_rhr,dipdirection,dipangle
    end type plane

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



end module orientations