module tracer_constants
    ! Value constants shared across libtracer: the missing-value sentinel and
    ! error flags. Kept separate from tracer_precision, which holds only kinds.

    use tracer_precision

    implicit none

    ! Missing value aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp
    real(wp), parameter :: MV     = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV_INT = int(MISSING_VALUE_DEFAULT)

    ! Error values
    real(wp), parameter :: ERR_DIST = 1E8_dp
    integer,  parameter :: ERR_IND  = -1

end module tracer_constants
