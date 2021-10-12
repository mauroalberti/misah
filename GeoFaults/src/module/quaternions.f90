module quaternions

    use iso_c_binding

    implicit none

    type, bind(c) :: quaternion
        real(c_double) :: q(0:3) ! quaternion components, the last is the rotation component
    end type quaternion

    real (c_double), parameter :: quat_normaliz_tolerance = 1.0e-6

contains

    !--------------------------

    subroutine quat_product(quat1, quat2, product)

        ! QUATq_product
        ! quaternion product
        ! Avenue created: 2005-02-11

        type(quaternion) :: quat1, quat2, product


        product%q(0) = (quat1%q(0)*quat2%q(0))   &
                  -(quat1%q(1)*quat2%q(1))	&
                  -(quat1%q(2)*quat2%q(2))	&
                  -(quat1%q(3)*quat2%q(3))

        product%q(1) = (quat1%q(0)*quat2%q(1))	&
                  +(quat1%q(1)*quat2%q(0))	&
                  +(quat1%q(2)*quat2%q(3))	&
                  -(quat1%q(3)*quat2%q(2))

        product%q(2) = (quat1%q(0)*quat2%q(2))	&
                  -(quat1%q(1)*quat2%q(3))	&
                  +(quat1%q(2)*quat2%q(0))	&
                  +(quat1%q(3)*quat2%q(1))

        product%q(3) = (quat1%q(0)*quat2%q(3))	&
                  +(quat1%q(1)*quat2%q(2))	&
                  -(quat1%q(2)*quat2%q(1))	&
                  +(quat1%q(3)*quat2%q(0))

    end subroutine quat_product

    !--------------------------

    subroutine quat_conjugate(quat, conjugate)

        ! QUATq_conjugate
        ! created 2005-02-10

        type(quaternion) :: quat, conjugate

        conjugate%q(0) =   quat%q(0)
        conjugate%q(1) = - quat%q(1)
        conjugate%q(2) = - quat%q(2)
        conjugate%q(3) = - quat%q(3)

    end subroutine quat_conjugate

    !--------------------------

    subroutine quat_squarednorm(quat, squarednorm)

        !  QUATq_squarednorm
        !  created 2005-02-12

        type(quaternion) :: quat
        real (c_double) :: squarednorm

        squarednorm = (quat%q(0) ** 2) + (quat%q(1) ** 2) + (quat%q(2) ** 2) + (quat%q(3) ** 2)

    end subroutine quat_squarednorm

    !--------------------------

    subroutine quat_scalardivision(quat, scaldiv, quat_scaldiv)

        ! QUATq_divisionbyscalar
        ! quaternion division by a scalar
        ! created: 2005-02-12

        integer (c_int) :: i
        type(quaternion), intent(in) :: quat, quat_scaldiv
        real (c_double), intent(in) :: scaldiv

        do i=0,3
            quat_scaldiv%q(i) = quat%q(i) / scaldiv
        end do

    end subroutine quat_scalardivision

    !--------------------------

    subroutine quat_inverse(quat, inverse)

        !  quaternion inverse
        !  created: 2005-02-19

        type(quaternion) :: quat, inverse
        real (c_double) :: squarednorm

        squarednorm = quat_squarednorm(quat_conjugate(quat))
        inverse = quat_scalardivision(quat_conj, squarednorm)

    end subroutine quat_inverse

    !--------------------------

    subroutine is_quat_normalized(quat, is_normalized)

        ! QUATq_normalizedquaterniontest
        ! created 2005-02-12

        real (c_double) :: squarednorm, abs_diff
        type(quaternion) :: quat
        logical(c_bool) :: is_normalized

        squarednorm = quat_squarednorm(quat)
        abs_diff = dabs(1 - squarednorm)

        if (abs_diff > quat_normaliz_tolerance) then
            is_normalized = .false.
        else
            is_normalized = .true.
        endif

    end subroutine is_quat_normalized

    !--------------------------

    subroutine quat_normalization(quat, normalized)

        !  QUATq_normalization
        !  transformation from quaternion to normalized quaternion
        !  created: 2005-02-14

        type(quaternion) :: quat, normalized
        real (c_double) :: norm

        norm = dsqrt(quat_squarednorm(quat))
        normalized = quat_scalardivision(quat, norm)

    end subroutine quat_normalization

    !--------------------------
    ! calculates a quaternion from a 3x3 matrix

    subroutine quaternfromcartmatr(focmec_matrix, quat)

        ! QUATq_TPBcartesmatrix2quaternion
        ! modified 2005-02-17

        real (c_double) :: focmec_matrix(3, 3)
        real (c_double) :: Q0, Q1, Q2, Q3
        real (c_double) :: Q0Q1, Q0Q2, Q0Q3, Q1Q2, Q1Q3, Q2Q3
        type(quaternion) :: quat

        ! myR11 = t1 = focmec_matrix(1,1)
        ! myR21 = t2 = focmec_matrix(2,1)
        ! myR31 = t3 = focmec_matrix(3,1)
        ! myR12 = p1 = focmec_matrix(1,2)
        ! myR22 = p2 = focmec_matrix(2,2)
        ! myR32 = p3 = focmec_matrix(3,2)
        ! myR13 = b1 = focmec_matrix(1,3)
        ! myR23 = b2 = focmec_matrix(2,3)
        ! myR33 = b3 = focmec_matrix(3,3)

        Q0 = 0.5*(dsqrt(1.0 + focmec_matrix(1, 1) + focmec_matrix(2, 2) + focmec_matrix(3, 3)))
        Q1 = 0.5*(dsqrt(1.0 + focmec_matrix(1, 1) - focmec_matrix(2, 2) - focmec_matrix(3, 3)))
        Q2 = 0.5*(dsqrt(1.0 - focmec_matrix(1, 1) + focmec_matrix(2, 2) - focmec_matrix(3, 3)))
        Q3 = 0.5*(dsqrt(1.0 - focmec_matrix(1, 1) - focmec_matrix(2, 2) + focmec_matrix(3, 3)))

        Q0Q1 = 0.25*(focmec_matrix(3, 2) - focmec_matrix(2, 3))
        Q0Q2 = 0.25*(focmec_matrix(1, 3) - focmec_matrix(3, 1))
        Q0Q3 = 0.25*(focmec_matrix(2, 1) - focmec_matrix(1, 2))
        Q1Q2 = 0.25*(focmec_matrix(1, 2) + focmec_matrix(2, 1))
        Q1Q3 = 0.25*(focmec_matrix(1, 3) + focmec_matrix(3, 1))
        Q2Q3 = 0.25*(focmec_matrix(2, 3) + focmec_matrix(3, 2))

        if((3 * Q0) > (Q1 + Q2 + Q3)) then
            Q1 = Q0Q1/Q0
            Q2 = Q0Q2/Q0
            Q3 = Q0Q3/Q0
        elseif ((3*Q1) > (Q0 + Q2 + Q3)) then
            Q0 = Q0Q1/Q1
            Q2 = Q1Q2/Q1
            Q3 = Q1Q3/Q1
        elseif ((3*Q2) > (Q0 + Q1 + Q3)) then
            Q0 = Q0Q2/Q2
            Q1 = Q1Q2/Q2
            Q3 = Q2Q3/Q2
        else
            Q0 = Q0Q3/Q3
            Q1 = Q1Q3/Q3
            Q2 = Q2Q3/Q3
        end if

        quat%q(0)= Q0
        quat%q(1)= Q1
        quat%q(2)= Q2
        quat%q(3)= Q3

    end subroutine quaternfromcartmatr

    !--------------------------

end module quaternions