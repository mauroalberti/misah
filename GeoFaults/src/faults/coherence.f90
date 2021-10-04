module coherence

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



end module coherence