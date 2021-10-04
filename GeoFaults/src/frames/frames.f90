module frames

    implicit none

    type :: orthonormal_triad
        sequence
        type(vector) :: X,Y,Z
    end type orthonormal_triad

    type :: triad_axes
        type(axis) :: axis_a, axis_b, axis_c
    end type triad_axes


    contains


end module frames