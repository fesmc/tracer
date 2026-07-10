
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
    use tracer

    implicit none

    type(tracer_class) :: trc1

    character(len=512) :: filename_nml
    character(len=512) :: file_topo, file_restart
    character(len=512) :: fldr, filename, filename_stats, filename_slice

    integer :: nx, ny, nz, nz_ac
    integer :: time_index
    integer, allocatable :: dims(:)

    real(prec), allocatable :: xc(:), yc(:), zeta(:), zeta_ac(:)
    real(prec), allocatable :: lon2D(:,:), lat2D(:,:)
    real(prec), allocatable :: z_srf(:,:), H_ice(:,:)

    ! Yelmo staggering: ux lives on the right border of H_ice(i,j) (acx nodes),
    ! uy on the top border (acy nodes), both on the zeta levels. uz lives at the
    ! cell centre horizontally but on the zeta_ac levels (layer interfaces).
    ! tracer_update wants all three at cell centres on a single z axis.
    real(prec), allocatable :: ux_acx(:,:,:), uy_acy(:,:,:), uz_ac(:,:,:)
    real(prec), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)

    integer    :: k, nstep
    real(prec) :: time, time_start, time_end
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
    allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz))

    call nc_read(file_restart,"xc",xc)
    call nc_read(file_restart,"yc",yc)
    call nc_read(file_restart,"zeta",zeta)
    call nc_read(file_restart,"zeta_ac",zeta_ac)

    ! Axes are stored in km; the tracer works in m.
    xc = xc*1e3
    yc = yc*1e3

    ! lon/lat come from the topography file, which defines the grid. Check it is
    ! the same grid as the restart rather than trusting the file names.
    call check_same_grid(file_topo,xc,yc)
    call nc_read(file_topo,"lon2D",lon2D)
    call nc_read(file_topo,"lat2D",lat2D)

    ! Geometry and velocity at one record of the restart's time axis.
    call nc_read(file_restart,"z_srf",z_srf, start=[1,1,time_index],  count=[nx,ny,1])
    call nc_read(file_restart,"H_ice",H_ice, start=[1,1,time_index],  count=[nx,ny,1])
    call nc_read(file_restart,"ux",ux_acx,   start=[1,1,1,time_index],count=[nx,ny,nz,1])
    call nc_read(file_restart,"uy",uy_acy,   start=[1,1,1,time_index],count=[nx,ny,nz,1])
    call nc_read(file_restart,"uz",uz_ac,    start=[1,1,1,time_index],count=[nx,ny,nz_ac,1])

    ! Yelmo's uz is already positive upward, so it is used as read. (The former
    ! GRISLI forcing was positive downward and had to be negated here.)
    call destagger_acx_to_aa(ux_acx,ux)
    call destagger_acy_to_aa(uy_acy,uy)
    call interp_zeta_ac_to_zeta(uz_ac,zeta_ac,zeta,uz)

    ! === Run ===================================================

    fldr           = "output/GRL-16KM"
    filename       = "GRL-16KM_trc1.nc"
    filename_stats = "GRL-16KM_trc1-stats.nc"
    filename_slice = "GRL-16KM_trc1-slice.nc"

    call tracer_init(trc1,filename_nml,time=real(time_start,prec_time), &
                     x=xc,y=yc,is_sigma=.TRUE.)
    call tracer_write_init(trc1,fldr,filename)

    nstep = int((time_end-time_start)/trc1%par%dt)

    do k = 0, nstep

        time = time_start + trc1%par%dt*k

        dep_now   = (mod(time,trc1%par%dt_dep)   .eq. 0.0)
        write_now = (mod(time,trc1%par%dt_write) .eq. 0.0)
        stats_now = (k .eq. nstep)

        ! Climate tagging fields (t2m, pr, d18O) are not carried by the restart,
        ! so they are omitted and recorded as missing. lon/lat are supplied.
        call tracer_update(trc1,time=real(time,prec_time), &
                           x=xc,y=yc,z=zeta,z_srf=z_srf,H=H_ice, &
                           ux=ux,uy=uy,uz=uz,lon=lon2D,lat=lat2D, &
                           dep_now=dep_now,stats_now=stats_now)

        if (write_now) then
            write(*,"(a,f12.1,a,i8)") " time = ", time, "   n_active = ", trc1%par%n_active
            call tracer_write(trc1,real(time,prec_time),fldr,filename)
        end if

    end do

    call tracer_write_stats(trc1,real(time_end,prec_time),fldr,filename_stats)

    ! A compact snapshot of just the active tracers at the final time.
    call tracer_write_slice(trc1,real(time_end,prec_time),fldr,filename_slice)

    if (trc1%par%n_active .eq. 0) then
        write(0,*) "test_greenland:: error: no tracers were deposited."
        error stop
    end if

    write(*,*) "test_greenland:: done. n_active = ", trc1%par%n_active

contains

    subroutine check_same_grid(filename,xc_ref,yc_ref)
        ! The topography and restart files must describe the same grid, since
        ! lon/lat are taken from one and the ice state from the other.

        implicit none

        character(len=*), intent(IN) :: filename
        real(prec),       intent(IN) :: xc_ref(:), yc_ref(:)

        real(prec), allocatable :: xc_chk(:), yc_chk(:)
        real(prec), parameter   :: tol = 1.0    ! m

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

    subroutine destagger_acx_to_aa(var_acx,var_aa)
        ! var_acx(i,j) sits on the right border of cell (i,j). The cell-centred
        ! value is the mean of the two bordering faces; at i=1 only the right
        ! face exists.

        implicit none

        real(prec), intent(IN)  :: var_acx(:,:,:)
        real(prec), intent(OUT) :: var_aa(:,:,:)

        integer :: i

        var_aa(1,:,:) = var_acx(1,:,:)

        do i = 2, size(var_acx,1)
            var_aa(i,:,:) = 0.5*(var_acx(i-1,:,:) + var_acx(i,:,:))
        end do

        return

    end subroutine destagger_acx_to_aa

    subroutine destagger_acy_to_aa(var_acy,var_aa)
        ! var_acy(i,j) sits on the top border of cell (i,j).

        implicit none

        real(prec), intent(IN)  :: var_acy(:,:,:)
        real(prec), intent(OUT) :: var_aa(:,:,:)

        integer :: j

        var_aa(:,1,:) = var_acy(:,1,:)

        do j = 2, size(var_acy,2)
            var_aa(:,j,:) = 0.5*(var_acy(:,j-1,:) + var_acy(:,j,:))
        end do

        return

    end subroutine destagger_acy_to_aa

    subroutine interp_zeta_ac_to_zeta(var_ac,z_ac,z,var)
        ! Linearly interpolate a field defined on the zeta_ac levels (layer
        ! interfaces) onto the zeta levels. Both axes are ascending and share
        ! their endpoints (0 = base, 1 = surface).

        implicit none

        real(prec), intent(IN)  :: var_ac(:,:,:)
        real(prec), intent(IN)  :: z_ac(:), z(:)
        real(prec), intent(OUT) :: var(:,:,:)

        integer    :: k, k_ac, n_ac
        real(prec) :: wt

        n_ac = size(z_ac)

        do k = 1, size(z)

            ! Bracket z(k) in z_ac: z_ac(k_ac) <= z(k) <= z_ac(k_ac+1)
            k_ac = 1
            do while (k_ac .lt. n_ac-1 .and. z_ac(k_ac+1) .lt. z(k))
                k_ac = k_ac + 1
            end do

            wt = (z(k) - z_ac(k_ac)) / (z_ac(k_ac+1) - z_ac(k_ac))

            var(:,:,k) = (1.0-wt)*var_ac(:,:,k_ac) + wt*var_ac(:,:,k_ac+1)

        end do

        return

    end subroutine interp_zeta_ac_to_zeta

end program tracertest
