module quaternions

    use iso_c_binding

    implicit none

    type, bind(c) :: quaternion
        real(c_double) :: q(0:3) ! quaternion components, the last is the rotation component
    end type quaternion

    logical(c_bool) :: is_q_normalized
    real (c_double) :: q_sqnorm, q_norm
    real (c_double), parameter :: q_norm_tolerance = 1.0e-6

    real (c_double) :: focmec_matrix(3,3)

contains

    !--------------------------

    type(quaternion) function q_product(quat1,quat2)

        ! QUATq_product
        ! quaternion product
        ! Avenue created: 2005-02-11

        type(quaternion) :: quat1,quat2

        q_product%q(0) =  (quat1%q(0)*quat2%q(0))   &
                -(quat1%q(1)*quat2%q(1))	&
                -(quat1%q(2)*quat2%q(2))	&
                -(quat1%q(3)*quat2%q(3))

        q_product%q(1) =  (quat1%q(0)*quat2%q(1))	&
                +(quat1%q(1)*quat2%q(0))	&
                +(quat1%q(2)*quat2%q(3))	&
                -(quat1%q(3)*quat2%q(2))

        q_product%q(2) =  (quat1%q(0)*quat2%q(2))	&
                -(quat1%q(1)*quat2%q(3))	&
                +(quat1%q(2)*quat2%q(0))	&
                +(quat1%q(3)*quat2%q(1))

        q_product%q(3) =  (quat1%q(0)*quat2%q(3))	&
                +(quat1%q(1)*quat2%q(2))	&
                -(quat1%q(2)*quat2%q(1))	&
                +(quat1%q(3)*quat2%q(0))

    end function q_product

    !--------------------------

    !--------------------------

    type(quaternion) function q_conjugate(quat)

        ! QUATq_conjugate
        ! created 2005-02-10

        type(quaternion) :: quat

        q_conjugate%q(0) =   quat%q(0)
        q_conjugate%q(1) = - quat%q(1)
        q_conjugate%q(2) = - quat%q(2)
        q_conjugate%q(3) = - quat%q(3)

    end function q_conjugate

    !--------------------------

    !--------------------------

    real (c_double) function q_sqrdnorm(quat)

        !  QUATq_squarednorm
        !  created 2005-02-12

        type(quaternion) :: quat

        q_sqrdnorm = (quat%q(0) ** 2) + (quat%q(1) ** 2) + (quat%q(2) ** 2) + (quat%q(3) ** 2)

    end function q_sqrdnorm

    !--------------------------

    !--------------------------

    type(quaternion) function q_scalardiv(quat, scaldiv)

        ! QUATq_divisionbyscalar
        ! quaternion division by a scalar
        ! created: 2005-02-12

        integer (c_int) :: i
        type(quaternion), intent(in) :: quat
        real (c_double), intent(in) :: scaldiv

        do i=0,3
            q_scalardiv%q(i) = quat%q(i)/scaldiv
        end do

    end function q_scalardiv

    !--------------------------

    !--------------------------

    type(quaternion) function q_inverse(quat)

        !  quaternion inverse
        !  created: 2005-02-19

        type(quaternion) :: quat, quat_conj

        quat_conj = q_conjugate(quat)
        q_sqnorm = q_sqrdnorm(quat)
        q_inverse = q_scalardiv(quat_conj, q_sqnorm)

    end function q_inverse

    !--------------------------

    !--------------------------

    logical(c_bool) function q_is_normalized(quat)

        ! QUATq_normalizedquaterniontest
        ! created 2005-02-12

        real (c_double) :: abs_diff
        type(quaternion) :: quat

        q_sqnorm = q_sqrdnorm(quat)
        abs_diff = dabs(1 - q_sqnorm)

        if (abs_diff > q_norm_tolerance) then
            q_is_normalized = .false.
        else
            q_is_normalized = .true.
        endif

    end function q_is_normalized

    !--------------------------

    !--------------------------

    type(quaternion) function q_normalized(quat)

        !  QUATq_normalization
        !  transformation from quaternion to normalized quaternion
        !  created: 2005-02-14

        type(quaternion) :: quat

        q_sqnorm = q_sqrdnorm(quat)
        q_norm = dsqrt(q_sqnorm)
        q_normalized = q_scalardiv(quat, q_norm)

    end function q_normalized

    !--------------------------

    !--------------------------
    ! calculates a quaternion from a 3x3 matrix

    type(quaternion) function q_from_matrix(focmec_matrix)

        ! QUATq_TPBcartesmatrix2quaternion
        ! modified 2005-02-17

        real (c_double) :: focmec_matrix(3,3)
        real (c_double) :: Q0, Q1, Q2, Q3
        real (c_double) :: Q0Q1,Q0Q2,Q0Q3,Q1Q2,Q1Q3,Q2Q3

        ! myR11 = t1 = focmec_matrix(1,1)
        ! myR21 = t2 = focmec_matrix(2,1)
        ! myR31 = t3 = focmec_matrix(3,1)
        ! myR12 = p1 = focmec_matrix(1,2)
        ! myR22 = p2 = focmec_matrix(2,2)
        ! myR32 = p3 = focmec_matrix(3,2)
        ! myR13 = b1 = focmec_matrix(1,3)
        ! myR23 = b2 = focmec_matrix(2,3)
        ! myR33 = b3 = focmec_matrix(3,3)

        Q0 = 0.5*(dsqrt(1+focmec_matrix(1,1)+focmec_matrix(2,2)+focmec_matrix(3,3)))
        Q1 = 0.5*(dsqrt(1+focmec_matrix(1,1)-focmec_matrix(2,2)-focmec_matrix(3,3)))
        Q2 = 0.5*(dsqrt(1-focmec_matrix(1,1)+focmec_matrix(2,2)-focmec_matrix(3,3)))
        Q3 = 0.5*(dsqrt(1-focmec_matrix(1,1)-focmec_matrix(2,2)+focmec_matrix(3,3)))

        Q0Q1 = 0.25*(focmec_matrix(3,2) - focmec_matrix(2,3))
        Q0Q2 = 0.25*(focmec_matrix(1,3) - focmec_matrix(3,1))
        Q0Q3 = 0.25*(focmec_matrix(2,1) - focmec_matrix(1,2))
        Q1Q2 = 0.25*(focmec_matrix(1,2) + focmec_matrix(2,1))
        Q1Q3 = 0.25*(focmec_matrix(1,3) + focmec_matrix(3,1))
        Q2Q3 = 0.25*(focmec_matrix(2,3) + focmec_matrix(3,2))

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

        q_from_matrix%q(0)= Q0
        q_from_matrix%q(1)= Q1
        q_from_matrix%q(2)= Q2
        q_from_matrix%q(3)= Q3

    end function q_from_matrix

    !--------------------------

end module quaternions