module rotations

    implicit none

    type :: rotation
        type(axis) :: rot_axis
        real(kind=r8b) :: rot_angle
    end type rotation

    contains


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


end module rotations