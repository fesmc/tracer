
program tracertest
    ! 3D Greenland benchmark: advect Lagrangian tracers through a Yelmo ice-sheet
    ! state on the 16 km Greenland grid, forced offline from a restart file.
    !
    ! The forcing is read at an explicit index along the restart's time
    ! dimension. A restart holds a single snapshot, but an offline transient run
    ! is the normal use case for this driver, so the time index is a parameter
    ! and the read is written to slice one record rather than to assume there is
    ! only ever one.

    use nml
    use ncio
    use coords, only : grid_class, grid_init
    use tracer

    implicit none

    type(tracer_class) :: trc1, trc2
    type(grid_class)   :: grd

    character(len=512) :: filename_nml
    character(len=512) :: file_topo, file_restart
    character(len=512) :: fldr, filename, filename_stats, filename_pd

    integer :: nx, ny, nz, nz_ac
    integer :: time_index
    integer, allocatable :: dims(:)

    real(wp), allocatable :: xc(:), yc(:), zeta(:), zeta_ac(:)
    real(wp), allocatable :: xc_ac(:), yc_ac(:)
    real(wp), allocatable :: lon2D(:,:), lat2D(:,:)
    real(wp), allocatable :: z_srf(:,:), H_ice(:,:)

    ! Yelmo staggering: ux lives on the right border of H_ice(i,j) (acx nodes),
    ! uy on the top border (acy nodes), both on the zeta levels. uz lives at the
    ! cell centre horizontally but on the zeta_ac levels (layer interfaces).
    ! These are passed to tracer_update on their native axes (xc_ac, yc_ac,
    ! zeta_ac) — no host-side destaggering to aa-nodes.
    real(wp), allocatable :: ux_acx(:,:,:), uy_acy(:,:,:), uz_ac(:,:,:)

    integer    :: k, nstep
    real(wp) :: time, time_start, time_end
    logical    :: dep_now, write_now, stats_now

    filename_nml = "Greenland.nml"

    call nml_read(filename_nml,"grl_par","file_topo",   file_topo)
    call nml_read(filename_nml,"grl_par","file_restart",file_restart)
    call nml_read(filename_nml,"grl_par","time_index",  time_index)
    call nml_read(filename_nml,"grl_par","time_start",  time_start)
    call nml_read(filename_nml,"grl_par","time_end",    time_end)

    ! === Load the grid and the forcing =========================

    call nc_dims(file_restart,"ux",dims=dims)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    call nc_dims(file_restart,"uz",dims=dims)
    nz_ac = dims(3)

    if (nz_ac .ne. nz+1) then
        write(0,*) "test_greenland:: error: expected size(zeta_ac) == size(zeta)+1"
        write(0,*) "nz = ", nz, ", nz_ac = ", nz_ac
        error stop
    end if

    allocate(xc(nx),yc(ny),zeta(nz),zeta_ac(nz_ac))
    allocate(lon2D(nx,ny),lat2D(nx,ny))
    allocate(z_srf(nx,ny),H_ice(nx,ny))
    allocate(ux_acx(nx,ny,nz),uy_acy(nx,ny,nz),uz_ac(nx,ny,nz_ac))
    allocate(xc_ac(nx),yc_ac(ny))

    call nc_read(file_restart,"xc",xc)
    call nc_read(file_restart,"yc",yc)
    call nc_read(file_restart,"zeta",zeta)
    call nc_read(file_restart,"zeta_ac",zeta_ac)

    ! Axes are stored in km; the tracer works in m.
    xc = xc*1e3
    yc = yc*1e3

    ! Staggered horizontal axes: acx nodes sit half a cell to the right of the
    ! aa nodes, acy nodes half a cell up. tracer_update samples ux/uy there.
    xc_ac = xc + 0.5*(xc(2)-xc(1))
    yc_ac = yc + 0.5*(yc(2)-yc(1))

    ! lon/lat come from the topography file, which defines the grid. Check it is
    ! the same grid as the restart rather than trusting the file names.
    call check_same_grid(file_topo,xc,yc)
    call nc_read(file_topo,"lon2D",lon2D)
    call nc_read(file_topo,"lat2D",lat2D)

    ! The projected grid (with lon/lat/area/crs) for the gridded-stats output.
    ! In a Yelmo coupling this object is owned by Yelmo and passed straight in.
    call grid_init(grd,filename=file_topo)

    ! Geometry and velocity at one record of the restart's time axis.
    call nc_read(file_restart,"z_srf",z_srf, start=[1,1,time_index],  count=[nx,ny,1])
    call nc_read(file_restart,"H_ice",H_ice, start=[1,1,time_index],  count=[nx,ny,1])
    call nc_read(file_restart,"ux",ux_acx,   start=[1,1,1,time_index],count=[nx,ny,nz,1])
    call nc_read(file_restart,"uy",uy_acy,   start=[1,1,1,time_index],count=[nx,ny,nz,1])
    call nc_read(file_restart,"uz",uz_ac,    start=[1,1,1,time_index],count=[nx,ny,nz_ac,1])

    ! Yelmo's uz is already positive upward, so it is used as read. (The former
    ! GRISLI forcing was positive downward and had to be negated here.) The
    ! velocities keep their native staggered locations; tracer_update handles
    ! the interpolation to each particle point.

    ! === Run ===================================================

    fldr            = "output/GRL-16KM"
    filename        = "tracer.nc"
    filename_stats  = "tracer-stats.nc"
    filename_pd     = "tracer-pd.nc"

    call tracer_init(trc1,filename_nml,time=real(time_start,prec_time), &
                     x=xc,y=yc,is_sigma=.TRUE.,grid=grd)
    call tracer_write_init(trc1,fldr,filename)
    call tracer_write_stats_init(trc1%stats,fldr,filename_stats)

    nstep = int((time_end-time_start)/trc1%par%dt)

    do k = 0, nstep

        time = time_start + trc1%par%dt*k

        dep_now   = (mod(time,trc1%par%dt_dep)         .eq. 0.0)
        write_now = (mod(time,trc1%par%dt_write)       .eq. 0.0)
        stats_now = (mod(time,trc1%par%dt_write_stats) .eq. 0.0)

        ! Climate tagging fields (t2m, pr, d18O) are not carried by the restart,
        ! so they are omitted and recorded as missing. lon/lat are supplied.
        call tracer_update(trc1,time=real(time,prec_time), &
                           x=xc,y=yc,z=zeta,z_srf=z_srf,H=H_ice, &
                           ux=ux_acx,uy=uy_acy,uz=uz_ac, &
                           x_ux=xc_ac,y_uy=yc_ac,z_uz=zeta_ac, &
                           lon=lon2D,lat=lat2D, &
                           dep_now=dep_now,stats_now=stats_now)

        if (write_now) then
            write(*,"(a,f12.1,a,i8)") " time = ", time, "   n_active = ", trc1%par%n_active
            call tracer_write(trc1,real(time,prec_time),fldr,filename)
        end if

        ! Append a gridded-statistics snapshot (the transient monitor).
        if (stats_now) call tracer_write_stats(trc1%stats,real(time,prec_time),fldr,filename_stats)

    end do

    ! Present-day product: isochrones, depth-layer stats and grid metadata.
    ! calc_tracer_stats ran on the final step (stats_now true at time_end = 0).
    call tracer_write_product(trc1%stats,fldr,filename_pd,H_ice=H_ice)

    if (trc1%par%n_active .eq. 0) then
        write(0,*) "test_greenland:: error: no tracers were deposited."
        error stop
    end if

    write(*,*) "test_greenland:: done. n_active = ", trc1%par%n_active

    ! === Restart round-trip self-check =========================
    ! Re-read the final record straight back into a fresh tracer object and
    ! confirm the state matches the in-memory object it was written from. A
    ! clean restart must reproduce every active tracer to within storage tol.
    call tracer_init(trc2,filename_nml,time=real(time_end,prec_time), &
                     x=xc,y=yc,is_sigma=.TRUE.,grid=grd)
    call tracer_read(trc2,trim(fldr)//"/"//trim(filename),time=real(time_end,prec_time))
    call check_restart(trc1,trc2)

contains

    subroutine check_restart(trc,trc_r)
        ! Assert that a restart-read object (trc_r) matches the original (trc)
        ! on every active tracer, to within single-precision round-trip tol.

        type(tracer_class), intent(IN) :: trc, trc_r

        real(wp) :: tol
        integer    :: n_mismatch

        tol = 1.0_wp        ! 1 m / 1 (m/a) — well above single-precision noise

        if (count(trc_r%now%active .ne. trc%now%active) .gt. 0) then
            write(0,*) "check_restart:: error: active mask differs after restart."
            error stop
        end if

        n_mismatch = 0
        call chk("x",     trc%now%x,     trc_r%now%x,     trc%now%active, tol, n_mismatch)
        call chk("y",     trc%now%y,     trc_r%now%y,     trc%now%active, tol, n_mismatch)
        call chk("z",     trc%now%z,     trc_r%now%z,     trc%now%active, tol, n_mismatch)
        call chk("dpth",  trc%now%dpth,  trc_r%now%dpth,  trc%now%active, tol, n_mismatch)
        call chk("z_srf", trc%now%z_srf, trc_r%now%z_srf, trc%now%active, tol, n_mismatch)
        call chk("ux",    trc%now%ux,    trc_r%now%ux,    trc%now%active, tol, n_mismatch)
        call chk("uy",    trc%now%uy,    trc_r%now%uy,    trc%now%active, tol, n_mismatch)
        call chk("uz",    trc%now%uz,    trc_r%now%uz,    trc%now%active, tol, n_mismatch)
        call chk("H",     trc%now%H,     trc_r%now%H,     trc%now%active, tol, n_mismatch)
        call chk("dep_time", trc%dep%time, trc_r%dep%time, trc%now%active, tol, n_mismatch)
        call chk("dep_z",    trc%dep%z,    trc_r%dep%z,    trc%now%active, tol, n_mismatch)

        if (count(trc_r%now%id     .ne. trc%now%id     .and. trc%now%active.gt.0) .gt. 0 .or. &
            count(trc_r%now%parent .ne. trc%now%parent .and. trc%now%active.gt.0) .gt. 0 .or. &
            count(trc_r%now%n_cloned .ne. trc%now%n_cloned .and. trc%now%active.gt.0) .gt. 0) then
            write(0,*) "check_restart:: error: lineage (id/parent/n_cloned) differs after restart."
            error stop
        end if

        if (n_mismatch .gt. 0) then
            write(0,*) "check_restart:: error: ", n_mismatch, " field(s) exceeded tol after restart."
            error stop
        end if

        write(*,*) "test_greenland:: restart round-trip OK (", &
                   count(trc%now%active.gt.0), "active tracers matched)."

        return

    end subroutine check_restart

    subroutine chk(name,a,b,active,tol,n_mismatch)
        character(len=*), intent(IN)    :: name
        real(wp),       intent(IN)    :: a(:), b(:)
        integer,          intent(IN)    :: active(:)
        real(wp),       intent(IN)    :: tol
        integer,          intent(INOUT) :: n_mismatch

        real(wp) :: dmax

        dmax = 0.0
        if (any(active.gt.0)) dmax = maxval(abs(a-b),mask=active.gt.0)
        if (dmax .gt. tol) then
            write(0,"(a,a,a,g14.6)") " check_restart:: ", trim(name), " max|diff| = ", dmax
            n_mismatch = n_mismatch + 1
        end if

        return

    end subroutine chk

    subroutine check_same_grid(filename,xc_ref,yc_ref)
        ! The topography and restart files must describe the same grid, since
        ! lon/lat are taken from one and the ice state from the other.

        implicit none

        character(len=*), intent(IN) :: filename
        real(wp),       intent(IN) :: xc_ref(:), yc_ref(:)

        real(wp), allocatable :: xc_chk(:), yc_chk(:)
        real(wp), parameter   :: tol = 1.0    ! m

        allocate(xc_chk(size(xc_ref)),yc_chk(size(yc_ref)))

        call nc_read(filename,"xc",xc_chk)
        call nc_read(filename,"yc",yc_chk)

        xc_chk = xc_chk*1e3
        yc_chk = yc_chk*1e3

        if (maxval(abs(xc_chk-xc_ref)) .gt. tol .or. &
            maxval(abs(yc_chk-yc_ref)) .gt. tol) then
            write(0,*) "check_same_grid:: error: grid mismatch with "//trim(filename)
            error stop
        end if

        return

    end subroutine check_same_grid

end program tracertest
