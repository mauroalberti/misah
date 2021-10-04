module focal_mechanisms

    implicit none


    contains


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


end module focal_mechanisms