module formats

    implicit none

    ! spatial information
    type :: spatialdataformat
        logical :: present
        integer (i1b) :: cod   ! 0: cartesian (x,y,z); 1: spherical (lat,lon,depth)
    end type spatialdataformat

    ! time information
    type :: timedataformat
        logical :: present
    end type timedataformat

    ! orientation information
    type :: orientationdataformat
        integer (i1b) :: cod  ! 0 - strike dip rake; 1 - strike dip slick. trend & plunge; 2 - P-axis trend & plunge T-axis trend & plunge
    end type orientationdataformat

    ! data format
    type :: dataformat
        type(spatialdataformat) :: spatial
        type(timedataformat) :: time
        type(orientationdataformat) :: orientation
        character (len=6) :: cod
    end type dataformat


contains


end module formats