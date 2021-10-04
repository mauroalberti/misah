module strain

    implicit none

    contains


        !------------------------
        !  incremental strain tensor calculation from strike, dip and rake

        function incremstraintensorcalc0(strike_rhr1,dipangle1,rake_aki1) result(incrstraintens)

            !  QUATc_momenttensorcalc
            !  Created: 2004-04-03

            real (kind=r8b), intent(in) :: strike_rhr1,dipangle1,rake_aki1
            real (kind=r8b) :: incrstraintens(3,3)

            !  formulas taken from: Aki & Richards, 1980, vol. 1, Box 4.4

            incrstraintens(1,1) = - sin(d2r*dipangle1)*cos(d2r*rake_aki1)*sin(2*d2r*strike_rhr1) &
                    - sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)*sin(d2r*strike_rhr1)**2
            incrstraintens(1,2) = sin(d2r*dipangle1)*cos(d2r*rake_aki1)*cos(2*d2r*strike_rhr1) &
                    + 0.5*sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)*sin(2*d2r*strike_rhr1)
            incrstraintens(1,3) = - cos(d2r*dipangle1)*cos(d2r*rake_aki1)*cos(d2r*strike_rhr1)  &
                    - cos(2*d2r*dipangle1)*sin(d2r*rake_aki1)*sin(d2r*strike_rhr1)
            incrstraintens(2,1) = incrstraintens(1,2)
            incrstraintens(2,2) = sin(d2r*dipangle1)*cos(d2r*rake_aki1)*sin(2*d2r*strike_rhr1) &
                    - sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)*cos(d2r*strike_rhr1)**2
            incrstraintens(2,3) = - cos(d2r*dipangle1)*cos(d2r*rake_aki1)*sin(d2r*strike_rhr1) &
                    + cos(2*d2r*dipangle1)*sin(d2r*rake_aki1)*cos(d2r*strike_rhr1)
            incrstraintens(3,1) = incrstraintens(1,3)
            incrstraintens(3,2) = incrstraintens(2,3)
            incrstraintens(3,3) = sin(2*d2r*dipangle1)*sin(d2r*rake_aki1)

        end function incremstraintensorcalc0

        !------------------------


        !------------------------
        !  incremental strain tensor calculation from focal mechanisms

        function incremstraintensorcalc1(Tvect1,Pvect1) result(incrstraintens)

            !  Created: 2008-01-27


            type(vector), intent(in) :: Tvect1, Pvect1
            real (kind=r8b) :: incrstraintens(3,3)

            !  formulas taken from: Kagan & Knopoff, 1985b, p. 649

            incrstraintens(1,1) = (Tvect1%x + Pvect1%x)*(Tvect1%x - Pvect1%x)
            incrstraintens(1,2) = (Tvect1%x + Pvect1%x)*(Tvect1%y - Pvect1%y)
            incrstraintens(1,3) = (Tvect1%x + Pvect1%x)*(Tvect1%z - Pvect1%z)
            incrstraintens(2,1) = (Tvect1%y + Pvect1%y)*(Tvect1%x - Pvect1%x)
            incrstraintens(2,2) = (Tvect1%y + Pvect1%y)*(Tvect1%y - Pvect1%y)
            incrstraintens(2,3) = (Tvect1%y + Pvect1%y)*(Tvect1%z - Pvect1%z)
            incrstraintens(3,1) = (Tvect1%z + Pvect1%z)*(Tvect1%x - Pvect1%x)
            incrstraintens(3,2) = (Tvect1%z + Pvect1%z)*(Tvect1%y - Pvect1%y)
            incrstraintens(3,3) = (Tvect1%z + Pvect1%z)*(Tvect1%z - Pvect1%z)

        end function incremstraintensorcalc1

        !------------------------

end module strain