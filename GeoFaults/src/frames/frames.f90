module frames

    implicit none

    type :: orthonormal_triad
        sequence
        type(vector) :: X,Y,Z
    end type orthonormal_triad

    type :: triad_axes
        type(axis) :: axis_a, axis_b, axis_c
    end type triad_axes

    type (triad_axes) ::  focmech_1, focmech_2, central_focalmechanism
    real (kind=r8b) :: focmec1_matrix(3,3), focmec2_matrix(3,3)

    real (kind=r8b) :: focmec1_matrix(3,3), focmec2_matrix(3,3)

    contains

end module frames