module quaternions

    implicit none

    type :: quaternion
        real(kind=r8b) :: q(0:3) ! quaternion components, the last is the rotation component
    end type quaternion

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


end module quaternions