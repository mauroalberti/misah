module points

    implicit none



    contains


        !--------------------------
        ! calculates the 3D euclidean distance between two points

        real(kind=r8b) function euclidean_distance(point1, point2)

            real (kind=r8b), dimension(3) :: point1, point2

            euclidean_distance = dsqrt((point2(1)-point1(1))**2 + (point2(2)-point1(2))**2 + (point2(3)-point1(3))**2)


        end function euclidean_distance

        !--------------------------


end module points