module tracer_stats
    ! Gridded (Eulerian) statistics of the Lagrangian tracer cloud.
    !
    ! Two products, split by when they carry meaning:
    !
    !   transient monitor  (tracer_write_stats)   -- appended along a time axis.
    !       Per (xc, yc, depth_norm): the number of tracers surviving in each
    !       depth band and the mean/sd of their deposition tags. A quick ncview
    !       check that a coupled run is behaving.
    !
    !   present-day product (tracer_write_product) -- one snapshot, at t = 0 ka BP.
    !       Adds the isochrone block (depth of each time_iso layer) and full grid
    !       metadata, so it is a drop-in comparison for a radiostratigraphy
    !       dataset such as MacGregor et al. (2015).
    !
    ! The statistics are formed in a single pass over the tracers
    ! (calc_tracer_stats). Sums are accumulated in double precision: dep_time
    ! runs to ~1e5 years and a single-precision one-pass variance would lose the
    ! whole signal to cancellation.
    !
    ! This module reads tracer state as plain arrays, never the tracer3D derived
    ! types, so it does not depend on tracer3D (which embeds tracer_stats_class).

    use tracer_precision
    use tracer_constants
    use ncio
    use coords, only : grid_class, grid_write

    implicit none

    private

    ! Deposition-tag fields carried in the depth-layer block. All are averaged
    ! over the tracers in a (cell, depth band); a field a run does not supply is
    ! all-MV and comes out missing.
    integer, parameter :: N_TAG   = 7
    integer, parameter :: F_DEPTIME = 1  ! deposition time (absolute, ka BP)
    integer, parameter :: F_DEPZ    = 2  ! deposition elevation (m)
    integer, parameter :: F_LON     = 3
    integer, parameter :: F_LAT     = 4
    integer, parameter :: F_T2M     = 5  ! annual mean temperature
    integer, parameter :: F_PR      = 6  ! annual precipitation
    integer, parameter :: F_D18O    = 7  ! annual d18O

    type stat_field_class
        character(len=32)  :: name
        character(len=32)  :: units
        character(len=128) :: long_name
        real(prec_wrt), allocatable :: mean(:,:,:)
        real(prec_wrt), allocatable :: sd(:,:,:)
    end type

    type tracer_stats_class
        ! Axes
        real(prec_wrt), allocatable :: x(:), y(:)         ! grid axes [m]
        real(prec_wrt), allocatable :: depth_norm(:)      ! normalized depth [1]
        real(prec_wrt), allocatable :: time_iso(:)        ! isochrone deposition times [ka]
        real(prec_wrt) :: dt_iso                          ! isochrone half-width [ka]

        ! Grid metadata for the output files (from the host model, e.g. Yelmo).
        logical          :: has_grid = .FALSE.
        type(grid_class) :: grid

        ! Depth-layer block (nx, ny, ndepth)
        integer,        allocatable :: count(:,:,:)       ! tracers in band (0 if none)
        type(stat_field_class)      :: fld(N_TAG)

        ! Isochrone block (nx, ny, niso)
        integer,        allocatable :: count_iso(:,:,:)
        real(prec_wrt), allocatable :: depth_iso(:,:,:)
        real(prec_wrt), allocatable :: depth_iso_sd(:,:,:)
        real(prec_wrt), allocatable :: dep_z_iso(:,:,:)
    end type

    public :: tracer_stats_class
    public :: tracer_stats_init, tracer_stats_end
    public :: calc_tracer_stats
    public :: tracer_write_stats_init, tracer_write_stats, tracer_write_product

contains

    subroutine tracer_stats_init(st, x, y, depth_norm, time_iso, dt_iso, grid)
        ! Allocate the gridded-stats object on the (x, y) grid, with the given
        ! depth and isochrone deposition-time axes. grid, if present, supplies the projection
        ! metadata (lon/lat/area/crs) written into the output files; without it
        ! the files carry bare xc/yc axes.

        implicit none

        type(tracer_stats_class), intent(INOUT) :: st
        real(wp),       intent(IN) :: x(:), y(:)
        real(prec_wrt),   intent(IN) :: depth_norm(:), time_iso(:)
        real(prec_wrt),   intent(IN) :: dt_iso
        type(grid_class), intent(IN), optional :: grid

        integer :: nx, ny, ndepth, niso, f

        call tracer_stats_end(st)

        nx = size(x); ny = size(y)
        ndepth = size(depth_norm); niso = size(time_iso)

        allocate(st%x(nx), st%y(ny))
        st%x = real(x, prec_wrt)
        st%y = real(y, prec_wrt)

        allocate(st%depth_norm(ndepth)); st%depth_norm = depth_norm
        allocate(st%time_iso(niso));     st%time_iso   = time_iso
        st%dt_iso = dt_iso

        st%has_grid = present(grid)
        if (present(grid)) st%grid = grid

        allocate(st%count(nx,ny,ndepth))
        do f = 1, N_TAG
            allocate(st%fld(f)%mean(nx,ny,ndepth), st%fld(f)%sd(nx,ny,ndepth))
        end do

        allocate(st%count_iso(nx,ny,niso))
        allocate(st%depth_iso(nx,ny,niso), st%depth_iso_sd(nx,ny,niso))
        allocate(st%dep_z_iso(nx,ny,niso))

        call set_field_meta(st)

        return

    end subroutine tracer_stats_init

    subroutine set_field_meta(st)
        ! Name, units and long name of each depth-layer tag field, for output.

        implicit none

        type(tracer_stats_class), intent(INOUT) :: st

        call meta(st%fld(F_DEPTIME), "dep_time",    "ka",            "Deposition time (age BP)")
        call meta(st%fld(F_DEPZ),    "dep_z",       "m",             "Deposition elevation")
        call meta(st%fld(F_LON),     "dep_lon",     "degrees_east",  "Deposition longitude")
        call meta(st%fld(F_LAT),     "dep_lat",     "degrees_north", "Deposition latitude")
        call meta(st%fld(F_T2M),     "dep_t2m_ann", "K",             "Deposition annual mean temperature")
        call meta(st%fld(F_PR),      "dep_pr_ann",  "m/a",           "Deposition annual precipitation")
        call meta(st%fld(F_D18O),    "dep_d18O_ann","permil",        "Deposition annual d18O")

        return

    contains

        subroutine meta(fld, name, units, long_name)
            type(stat_field_class), intent(INOUT) :: fld
            character(len=*), intent(IN) :: name, units, long_name
            fld%name = name; fld%units = units; fld%long_name = long_name
        end subroutine meta

    end subroutine set_field_meta

    subroutine tracer_stats_end(st)

        implicit none

        type(tracer_stats_class), intent(INOUT) :: st
        integer :: f

        if (allocated(st%x))            deallocate(st%x)
        if (allocated(st%y))            deallocate(st%y)
        if (allocated(st%depth_norm))   deallocate(st%depth_norm)
        if (allocated(st%time_iso))     deallocate(st%time_iso)
        if (allocated(st%count))        deallocate(st%count)
        do f = 1, N_TAG
            if (allocated(st%fld(f)%mean)) deallocate(st%fld(f)%mean)
            if (allocated(st%fld(f)%sd))   deallocate(st%fld(f)%sd)
        end do
        if (allocated(st%count_iso))    deallocate(st%count_iso)
        if (allocated(st%depth_iso))    deallocate(st%depth_iso)
        if (allocated(st%depth_iso_sd)) deallocate(st%depth_iso_sd)
        if (allocated(st%dep_z_iso))    deallocate(st%dep_z_iso)

        st%has_grid = .FALSE.

        return

    end subroutine tracer_stats_end

    subroutine calc_tracer_stats(st, active, px, py, dpth, H, &
                                 dep_time, dep_z, dep_lon, dep_lat, &
                                 dep_t2m_ann, dep_pr_ann, dep_d18O_ann)
        ! One pass over the tracers, binning them by grid cell and either
        ! normalized depth (the depth-layer block) or deposition time (the
        ! isochrone block), and accumulating the mean/sd of each field. Only
        ! deposited tracers (active == 2) count. A tag equal to MV is skipped, so
        ! a field a run does not carry comes out all-missing rather than polluting
        ! the mean. All arrays are indexed by tracer slot.

        implicit none

        type(tracer_stats_class), intent(INOUT) :: st
        integer,    intent(IN) :: active(:)
        real(wp), intent(IN) :: px(:), py(:), dpth(:), H(:)
        real(wp), intent(IN) :: dep_time(:), dep_z(:), dep_lon(:), dep_lat(:)
        real(wp), intent(IN) :: dep_t2m_ann(:), dep_pr_ann(:), dep_d18O_ann(:)

        integer :: nx, ny, ndepth, niso
        integer :: p, i, j, q, qi, f
        real(prec_wrt) :: sig
        real(dp) :: dep_time_ka, v
        real(wp) :: vals(N_TAG)

        ! Double-precision accumulators (freed on return)
        real(dp), allocatable :: acc_sum(:,:,:,:), acc_sq(:,:,:,:)
        integer,  allocatable :: fcount(:,:,:,:)
        real(dp), allocatable :: iso_dsum(:,:,:), iso_dsq(:,:,:), iso_zsum(:,:,:)

        nx = size(st%x); ny = size(st%y)
        ndepth = size(st%depth_norm); niso = size(st%time_iso)

        allocate(acc_sum(nx,ny,ndepth,N_TAG), acc_sq(nx,ny,ndepth,N_TAG))
        allocate(fcount(nx,ny,ndepth,N_TAG))
        allocate(iso_dsum(nx,ny,niso), iso_dsq(nx,ny,niso), iso_zsum(nx,ny,niso))

        st%count     = 0
        st%count_iso = 0
        acc_sum = 0.0_dp; acc_sq = 0.0_dp; fcount = 0
        iso_dsum = 0.0_dp; iso_dsq = 0.0_dp; iso_zsum = 0.0_dp

        do p = 1, size(active)

            if (active(p) .ne. 2) cycle
            if (H(p) .le. 0.0)    cycle

            i = nearest_index(st%x, real(px(p), prec_wrt), bounded=.TRUE.)
            j = nearest_index(st%y, real(py(p), prec_wrt), bounded=.TRUE.)
            if (i .eq. 0 .or. j .eq. 0) cycle

            ! Assemble the tag values in their output units, MV preserved.
            vals(F_DEPTIME) = to_ka(dep_time(p))
            vals(F_DEPZ)    = dep_z(p)
            vals(F_LON)     = dep_lon(p)
            vals(F_LAT)     = dep_lat(p)
            vals(F_T2M)     = dep_t2m_ann(p)
            vals(F_PR)      = dep_pr_ann(p)
            vals(F_D18O)    = dep_d18O_ann(p)

            ! -- depth-layer binning (nearest level; shallow tracers fall in band 1)
            sig = real(dpth(p)/H(p), prec_wrt)
            q = nearest_index(st%depth_norm, sig, bounded=.FALSE.)
            if (q .ge. 1 .and. q .le. ndepth) then
                st%count(i,j,q) = st%count(i,j,q) + 1
                do f = 1, N_TAG
                    if (vals(f) .ne. MV) then
                        v = real(vals(f), dp)
                        acc_sum(i,j,q,f) = acc_sum(i,j,q,f) + v
                        acc_sq(i,j,q,f)  = acc_sq(i,j,q,f)  + v*v
                        fcount(i,j,q,f)  = fcount(i,j,q,f)  + 1
                    end if
                end do
            end if

            ! -- isochrone binning (deposition time within dt_iso of a target; may match none)
            if (dep_time(p) .ne. MV) then
                dep_time_ka = real(dep_time(p), dp)*1e-3_dp
                do qi = 1, niso
                    if (abs(dep_time_ka - real(st%time_iso(qi), dp)) .le. real(st%dt_iso, dp)) then
                        st%count_iso(i,j,qi) = st%count_iso(i,j,qi) + 1
                        iso_dsum(i,j,qi) = iso_dsum(i,j,qi) + real(dpth(p), dp)
                        iso_dsq(i,j,qi)  = iso_dsq(i,j,qi)  + real(dpth(p), dp)**2
                        if (dep_z(p) .ne. MV) iso_zsum(i,j,qi) = iso_zsum(i,j,qi) + real(dep_z(p), dp)
                    end if
                end do
            end if

        end do

        ! -- finalize depth-layer block
        do f = 1, N_TAG
            call finalize_meansd(fcount(:,:,:,f), acc_sum(:,:,:,f), acc_sq(:,:,:,f), &
                                 st%fld(f)%mean, st%fld(f)%sd)
        end do

        ! -- finalize isochrone block
        call finalize_meansd(st%count_iso, iso_dsum, iso_dsq, st%depth_iso, st%depth_iso_sd)
        where (st%count_iso .gt. 0)
            st%dep_z_iso = real(iso_zsum / max(st%count_iso, 1), prec_wrt)
        elsewhere
            st%dep_z_iso = MV
        end where

        deallocate(acc_sum, acc_sq, fcount, iso_dsum, iso_dsq, iso_zsum)

        return

    contains

        elemental function to_ka(t) result(t_ka)
            ! Convert a deposition time [years] to ka, preserving MV.
            real(wp), intent(IN) :: t
            real(wp) :: t_ka
            if (t .eq. MV) then
                t_ka = MV
            else
                t_ka = t*1e-3
            end if
        end function to_ka

    end subroutine calc_tracer_stats

    subroutine finalize_meansd(cnt, ssum, ssq, mean, sd)
        ! Turn accumulated count / sum / sum-of-squares into mean and (sample)
        ! standard deviation. Empty cells (cnt == 0) are missing; a single sample
        ! has sd = 0. The variance is clamped at zero against roundoff.

        implicit none

        integer,        intent(IN)  :: cnt(:,:,:)
        real(dp),       intent(IN)  :: ssum(:,:,:), ssq(:,:,:)
        real(prec_wrt), intent(OUT) :: mean(:,:,:), sd(:,:,:)

        integer  :: i, j, k
        real(dp) :: m, var

        do k = 1, size(cnt,3)
        do j = 1, size(cnt,2)
        do i = 1, size(cnt,1)

            if (cnt(i,j,k) .le. 0) then
                mean(i,j,k) = MV
                sd(i,j,k)   = MV
            else
                m = ssum(i,j,k) / cnt(i,j,k)
                mean(i,j,k) = real(m, prec_wrt)
                if (cnt(i,j,k) .eq. 1) then
                    sd(i,j,k) = 0.0
                else
                    var = (ssq(i,j,k) - ssum(i,j,k)**2 / cnt(i,j,k)) / (cnt(i,j,k) - 1)
                    sd(i,j,k) = real(sqrt(max(var, 0.0_dp)), prec_wrt)
                end if
            end if

        end do
        end do
        end do

        return

    end subroutine finalize_meansd

    pure function nearest_index(axis, val, bounded) result(idx)
        ! Index of the axis point nearest to val (axis ascending). With
        ! bounded=.TRUE., returns 0 if val lies more than half a cell beyond
        ! either end; with bounded=.FALSE., clamps to the end points (so a
        ! tracer shallower than depth_norm(1) falls in band 1).

        implicit none

        real(prec_wrt), intent(IN) :: axis(:), val
        logical,        intent(IN) :: bounded
        integer :: idx, n, lo, hi, mid

        n = size(axis)
        idx = 0
        if (n .lt. 1) return
        if (n .eq. 1) then; idx = 1; return; end if

        if (bounded) then
            if (val .lt. axis(1) - 0.5*(axis(2)-axis(1)))   return
            if (val .gt. axis(n) + 0.5*(axis(n)-axis(n-1))) return
        end if

        lo = 1; hi = n
        do while (hi - lo .gt. 1)
            mid = (lo + hi)/2
            if (axis(mid) .le. val) then
                lo = mid
            else
                hi = mid
            end if
        end do

        if (abs(val - axis(lo)) .le. abs(axis(hi) - val)) then
            idx = lo
        else
            idx = hi
        end if

        return

    end function nearest_index

    subroutine tracer_write_stats_init(st, fldr, filename)
        ! Create the transient monitor file: grid axes (with projection metadata
        ! if a grid was supplied), the depth_norm axis, and an unlimited time
        ! axis for tracer_write_stats to append to.

        implicit none

        type(tracer_stats_class), intent(IN) :: st
        character(len=*), intent(IN) :: fldr, filename

        character(len=512) :: path

        path = trim(fldr)//"/"//trim(filename)

        call create_grid(st, path)
        call nc_write_dim(path, "depth_norm", x=st%depth_norm, units="1")
        call nc_write_dim(path, "time", x=real(MV,prec_wrt), unlimited=.TRUE., units="ka")

        return

    end subroutine tracer_write_stats_init

    subroutine tracer_write_stats(st, time, fldr, filename)
        ! Append one time record of the depth-layer block to the monitor file.

        implicit none

        type(tracer_stats_class), intent(IN) :: st
        real(prec_time),  intent(IN) :: time
        character(len=*), intent(IN) :: fldr, filename

        character(len=512) :: path
        integer, allocatable :: dims(:)
        integer :: nt, nx, ny, nd, f
        real(prec_wrt) :: time_ka, time_in

        path = trim(fldr)//"/"//trim(filename)
        nx = size(st%x); ny = size(st%y); nd = size(st%depth_norm)
        time_ka = real(time, prec_wrt)*1e-3

        ! Next time index (mirrors tracer_write's append idiom)
        call nc_dims(path, "time", dims=dims)
        nt = dims(1)
        call nc_read(path, "time", time_in, start=[nt], count=[1])
        if (time_in .ne. MV .and. abs(time_ka - time_in) .gt. 1e-5) nt = nt + 1

        call nc_write(path, "time", time_ka, dim1="time", start=[nt], count=[1], missing_value=real(MV,prec_wrt))

        call nc_write(path, "count", st%count, dim1="xc", dim2="yc", dim3="depth_norm", dim4="time", &
                      start=[1,1,1,nt], count=[nx,ny,nd,1], missing_value=int(MV), long_name="Tracer count")

        do f = 1, N_TAG
            call write_layer(path, trim(st%fld(f)%name),        st%fld(f)%mean, &
                             trim(st%fld(f)%units), trim(st%fld(f)%long_name), nx, ny, nd, nt)
            call write_layer(path, trim(st%fld(f)%name)//"_sd", st%fld(f)%sd, &
                             trim(st%fld(f)%units), trim(st%fld(f)%long_name)//" (sd)", nx, ny, nd, nt)
        end do


        return

    contains

        subroutine write_layer(path, name, var, units, long_name, nx, ny, nd, nt)
            character(len=*), intent(IN) :: path, name, units, long_name
            real(prec_wrt),   intent(IN) :: var(:,:,:)
            integer,          intent(IN) :: nx, ny, nd, nt
            call nc_write(path, name, var, dim1="xc", dim2="yc", dim3="depth_norm", dim4="time", &
                          start=[1,1,1,nt], count=[nx,ny,nd,1], missing_value=real(MV,prec_wrt), &
                          units=units, long_name=long_name)
        end subroutine write_layer

    end subroutine tracer_write_stats

    subroutine tracer_write_product(st, fldr, filename, H_ice)
        ! Write the present-day product: the depth-layer block, the isochrone
        ! block, grid metadata, and (optionally) ice thickness. One snapshot, no
        ! time axis. Intended to be called once, at t = 0 ka BP.

        implicit none

        type(tracer_stats_class), intent(IN) :: st
        character(len=*), intent(IN) :: fldr, filename
        real(wp), intent(IN), optional :: H_ice(:,:)

        character(len=512) :: path
        integer :: f

        path = trim(fldr)//"/"//trim(filename)

        call create_grid(st, path)
        call nc_write_dim(path, "depth_norm", x=st%depth_norm, units="1")
        call nc_write_dim(path, "time_iso",   x=st%time_iso,   units="ka")

        ! -- depth-layer block
        call nc_write(path, "count", st%count, dim1="xc", dim2="yc", dim3="depth_norm", &
                      missing_value=int(MV), long_name="Tracer count")
        do f = 1, N_TAG
            call nc_write(path, trim(st%fld(f)%name), st%fld(f)%mean, &
                          dim1="xc", dim2="yc", dim3="depth_norm", missing_value=real(MV,prec_wrt), &
                          units=trim(st%fld(f)%units), long_name=trim(st%fld(f)%long_name))
            call nc_write(path, trim(st%fld(f)%name)//"_sd", st%fld(f)%sd, &
                          dim1="xc", dim2="yc", dim3="depth_norm", missing_value=real(MV,prec_wrt), &
                          units=trim(st%fld(f)%units), long_name=trim(st%fld(f)%long_name)//" (sd)")
        end do
        ! -- isochrone block
        call nc_write(path, "depth_iso", st%depth_iso, dim1="xc", dim2="yc", dim3="time_iso", &
                      missing_value=real(MV,prec_wrt), units="m", long_name="Isochrone depth")
        call nc_write(path, "depth_iso_err", st%depth_iso_sd, dim1="xc", dim2="yc", dim3="time_iso", &
                      missing_value=real(MV,prec_wrt), units="m", long_name="Isochrone depth (sd)")
        call nc_write(path, "dep_z_iso", st%dep_z_iso, dim1="xc", dim2="yc", dim3="time_iso", &
                      missing_value=real(MV,prec_wrt), units="m", long_name="Isochrone deposition elevation")
        call nc_write(path, "count_iso", st%count_iso, dim1="xc", dim2="yc", dim3="time_iso", &
                      missing_value=int(MV), long_name="Tracer count (isochrones)")

        if (present(H_ice)) &
            call nc_write(path, "H_ice", real(H_ice,prec_wrt), dim1="xc", dim2="yc", &
                          missing_value=real(MV,prec_wrt), units="m", long_name="Ice thickness")

        return

    end subroutine tracer_write_product

    subroutine create_grid(st, path)
        ! Create the output file and write its horizontal axes. With a grid
        ! object this is the full projection metadata (x2D, lon2D, area, crs);
        ! without one, bare xc/yc in km.

        implicit none

        type(tracer_stats_class), intent(IN) :: st
        character(len=*),         intent(IN) :: path

        if (st%has_grid) then
            call grid_write(st%grid, path, "xc", "yc", create=.TRUE.)
        else
            call nc_create(path)
            call nc_write_dim(path, "xc", x=st%x*1e-3, units="km")
            call nc_write_dim(path, "yc", x=st%y*1e-3, units="km")
        end if

        return

    end subroutine create_grid

end module tracer_stats
