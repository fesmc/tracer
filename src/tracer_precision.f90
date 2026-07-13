module tracer_precision
    ! Numerical kind parameters used throughout libtracer.
    ! Value constants (missing value, error flags) live in tracer_constants.

    implicit none

    integer,  parameter :: sp  = kind(1.0)
    integer,  parameter :: dp  = kind(1.d0)

    ! Working precision (formerly `prec`)
    integer,  parameter :: wp        = sp
    integer,  parameter :: prec_time = sp
    integer,  parameter :: prec_wrt  = sp

end module tracer_precision
