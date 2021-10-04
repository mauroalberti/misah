module analyses

    implicit none

    ! settings of analysis session
    type :: analysistype
        type(statistics_type) :: statistics
        logical :: focmechrot		! focal mechanism rotation statistics
        logical :: kinematicaxesrot	! T-P-B axes rotations
        logical :: faultelementsrot	! fault sub-element rotation statistics
        logical :: coherence		! coherence statistics
        logical :: spatialsep		! spatial separation statistics
        logical :: timelag			! time lag statistics
    end type analysistype



    ! parameters of file
    type :: file_params
        character (len=255) :: name
        integer (kind=i1b) :: statuscod ! 0: replace
        integer :: unitnum
        integer (kind=i2b) :: headerrownumber
        integer (kind=i4b) :: recsnumber, reclength
    end type file_params


    ! analysis results
    type :: analysis_results
        type(rotation) :: focmech_rots(4), faultpl_rot, slickenl_rot
        type(rotation) :: Taxis_rot, Paxis_rot, Baxis_rot
        real (kind=r8b) :: coherence, time_lag, spat_sep
    end type analysis_results


contains


end module analyses