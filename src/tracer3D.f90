
module tracer3D 

    use tracer_precision
    use tracer_constants
    use tracer_interp
    use tracer_stats
    use bspline_module, only : bspline_3d
    use nml
    use coords, only : grid_class

    implicit none

    type tracer_par_trans_class
        integer :: nt 

        real(wp), allocatable :: time(:)
        real(wp), allocatable :: H_min_dep(:)
        real(wp), allocatable :: dt_dep(:)
        integer,    allocatable :: n_max_dep(:) 
        real(wp), allocatable :: dt_write(:)
        
    end type 

    type tracer_par_class 
        integer :: n, n_active, n_max_dep, id_max
        logical :: is_sigma                     ! Is the defined z-axis in sigma coords
        logical :: is_profile                   ! Is this a 2D (x-z) domain with a ghost y-axis
        real(prec_time) :: dt, dt_dep, dt_write
        real(wp) :: H_min                     ! Minimum ice thickness to track (m)
        real(wp) :: depth_max                 ! Maximum depth of tracer (fraction)
        real(wp) :: U_max                     ! Maximum horizontal velocity of tracer to track (m/a)
        real(wp) :: U_max_dep                 ! Maximum horizontal velocity allowed for tracer deposition (m/a)
        real(wp) :: H_min_dep                 ! Minimum ice thickness for tracer deposition (m)
        real(wp) :: alpha                     ! Slope of probability function ("linear"/"quadratic" weights)
        character(len=56) :: weight             ! Deposition prob. distribution: vel, linear, quadratic, rand
        logical    :: noise                     ! Add noise to gridded deposition location
        integer    :: seed                      ! RNG seed; <= 0 => nondeterministic
        character(len=56) :: interp_method

        ! Particle cloning: spawn offset copies of deep, slow, well-aged tracers
        ! to raise the odds that old deposition times survive to present day.
        logical    :: clone                     ! Enable cloning
        integer    :: n_clones                  ! Clones spawned per eligible tracer
        real(wp) :: clone_depth_frac          ! Min depth (fraction of H) to be eligible
        real(wp) :: clone_U_max               ! Max horizontal velocity to be eligible (m/a)
        real(wp) :: clone_dep_time_min        ! Deposition-time window, lower bound (years)
        real(wp) :: clone_dep_time_max        ! Deposition-time window, upper bound (years)
        real(wp) :: clone_offset_xy           ! Half-width of uniform horizontal offset (m)
        real(wp) :: clone_offset_z            ! Max upward vertical offset (m)

        ! Transient parameters
        character(len=512) :: par_trans_file
        logical            :: use_par_trans
        type(tracer_par_trans_class) :: tpar

        ! Gridded statistics (Eulerian output). Read only when stats = .TRUE.;
        ! see tracer_stats. depth_norm is built as n_depth uniform levels.
        logical            :: stats
        real(prec_time)    :: dt_write_stats
        integer            :: n_depth
        real(prec_wrt), allocatable :: time_iso(:)
        real(prec_wrt)     :: dt_iso

    end type

    type tracer_state_class
        real(prec_time) :: time, time_old
        real(prec_time) :: time_dep, time_write 
        real(prec_time) :: dt  
        integer, allocatable :: active(:), id(:)
        integer, allocatable :: parent(:)            ! id of the parent tracer if this is a clone, 0 if an original
        integer, allocatable :: n_cloned(:)          ! number of clones already spawned from this tracer
        real(wp), allocatable :: x(:), y(:), z(:), sigma(:)
        real(wp), allocatable :: ux(:), uy(:), uz(:)
        real(wp), allocatable :: ax(:), ay(:), az(:)
        real(wp), allocatable :: dpth(:), z_srf(:)
        real(wp), allocatable :: thk(:)            ! Tracer thickness (for compression)
        real(wp), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
        real(wp), allocatable :: H(:)

    end type 

    type tracer_dep_class
        ! Standard deposition information (time and place)
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: H(:)
        real(wp), allocatable :: x(:), y(:), z(:)
        real(wp), allocatable :: lon(:), lat(:)

        ! Monthly climate/isotope tags at deposition, shape (n, nmon). These are
        ! the inputs the host supplies; the annual quantities below are derived
        ! from them at deposition time.
        real(wp), allocatable :: t2m(:,:), pr(:,:), d18O(:,:)

        ! Annual quantities derived from the monthly tags at deposition:
        !   t2m_ann    - mean of monthly t2m
        !   pr_ann     - mean of monthly pr
        !   t2m_prann  - precip-weighted mean temperature
        !   d18O_ann   - mean of monthly d18O
        !   d18O_prann - precip-weighted mean d18O
        real(wp), allocatable :: t2m_ann(:), pr_ann(:), t2m_prann(:)
        real(wp), allocatable :: d18O_ann(:), d18O_prann(:)

    end type

    type tracer_class 
        type(tracer_par_class)   :: par 
        type(tracer_state_class) :: now 
        type(tracer_dep_class)   :: dep
        type(tracer_stats_class) :: stats 

    end type 

    type(lin3_interp_par_type) :: par_lin
    type(bspline_3d)           :: bspline3d_ux, bspline3d_uy, bspline3d_uz

    ! Number of months in the deposition climate tags
    integer, parameter :: nmon = 12

    private

    ! For other tracer modules
    public :: nmon
    public :: tracer_par_class
    public :: tracer_state_class
    public :: tracer_dep_class
    public :: tracer_stats_class

    public :: tracer_reshape1D_vec
    public :: tracer_reshape2D_field
    public :: tracer_reshape3D_field
    public :: which

    ! General public 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 
    public :: tracer_end 

contains 

    subroutine tracer_init(trc,filename,time,x,y,is_sigma,grid,time_iso)

        implicit none

        type(tracer_class),   intent(OUT) :: trc
        character(len=*),     intent(IN)  :: filename
        real(wp), intent(IN) :: x(:), y(:)
        logical,    intent(IN) :: is_sigma
        real(prec_time), intent(IN) :: time
        ! Grid metadata for the gridded-stats output; only meaningful when the
        ! namelist sets stats = .TRUE., and even then optional (x/y is the
        ! fallback). Supplied by the host model (e.g. Yelmo owns the grid).
        type(grid_class), intent(IN), optional :: grid
        ! Isochrone deposition-time targets [ka] for the gridded stats. When the
        ! host owns the isochrone axis (e.g. Yelmo's ytrc), it passes them here
        ! and they take precedence over the namelist time_iso, so the two never
        ! need to be kept in sync. Passed in wp (matching x/y); stored as prec_wrt.
        real(wp), intent(IN), optional :: time_iso(:)

        ! Load the parameters
        call tracer_par_load(trc%par,filename,is_sigma)

        ! Update the transient parameters
        if (trc%par%use_par_trans) then
            call tracer_par_update(trc%par,trc%par%tpar,time)
        end if

        ! Allocate the state variables
        call tracer_allocate(trc%now,trc%dep,n=trc%par%n)

        ! Gridded statistics are opt-in (namelist stats flag). depth_norm is
        ! n_depth uniform levels; dt_iso, dt_write_stats come from par. The
        ! isochrone targets come from the host (time_iso argument) when present,
        ! otherwise from the namelist (par%time_iso).
        if (trc%par%stats) then
            if (present(time_iso)) then
                if (allocated(trc%par%time_iso)) deallocate(trc%par%time_iso)
                allocate(trc%par%time_iso(size(time_iso)))
                trc%par%time_iso = real(time_iso, prec_wrt)
            end if
            if (.not. allocated(trc%par%time_iso)) then
                write(*,*) "tracer_init:: Error: stats=.TRUE. requires isochrone targets, &
                           &either via the namelist time_iso or the time_iso argument."
                stop "Program stopped."
            end if
            call tracer_stats_init(trc%stats,x,y, &
                                   depth_norm=uniform_depth(trc%par%n_depth), &
                                   time_iso=trc%par%time_iso, dt_iso=trc%par%dt_iso, &
                                   grid=grid)
        end if

        ! Initialize state
        trc%now%active    = 0 

        trc%now%id        = mv
        trc%now%parent    = mv
        trc%now%n_cloned  = mv
        trc%now%x         = mv
        trc%now%y         = mv
        trc%now%z         = mv
        trc%now%sigma     = mv
        trc%now%z_srf     = mv
        trc%now%dpth      = mv 
        trc%now%ux        = mv 
        trc%now%uy        = mv 
        trc%now%uz        = mv 
        trc%now%ax        = mv 
        trc%now%ay        = mv 
        trc%now%az        = mv 
        trc%now%thk       = mv 
        trc%now%T         = mv 
        trc%now%H         = mv 

        trc%dep%time      = mv 
        trc%dep%H         = mv 
        trc%dep%x         = mv 
        trc%dep%y         = mv 
        trc%dep%z         = mv 

        trc%par%id_max    = 0 

        ! Initialize the time (to one older than now)
        trc%now%time   = time - 1000.0_dp
        trc%now%time_dep   = time - 1000.0 
        trc%now%time_write = time - 1000.0 

        ! Initialize random number generator
        call tracer_set_seed(trc%par%seed)

        return

    end subroutine tracer_init

    subroutine tracer_update(trc,time,x,y,z,z_srf,H,ux,uy,uz,              &
                             x_ux,y_uy,z_uz,                                &
                             lon,lat,t2m,pr,d18O,d18O_ann,                  &
                             dep_now,stats_now,order,sigma_srf)

        implicit none

        type(tracer_class), intent(INOUT) :: trc
        real(prec_time), intent(IN) :: time
        real(wp), intent(IN) :: x(:), y(:), z(:)
        real(wp), intent(IN) :: z_srf(:,:), H(:,:)
        real(wp), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)

        ! Native staggered velocity axes (all optional). When present, the
        ! matching component is interpolated directly on its Arakawa-C location
        ! instead of requiring the host to pre-destagger it to aa-nodes:
        !   x_ux - acx x-axis for ux (ux staggered in x; y and sigma stay aa)
        !   y_uy - acy y-axis for uy (uy staggered in y; x and sigma stay aa)
        !   z_uz - zeta_ac sigma/height axis for uz (uz staggered in z; x and y
        !          stay aa). These interface levels are typically one longer
        !          than z, so size(uz,3) may differ from size(ux,3).
        ! Any axis left absent falls back to the aa axis (x, y, z), so the pure
        ! aa-node call is unchanged.
        real(wp), intent(IN), optional :: x_ux(:), y_uy(:), z_uz(:)

        ! Deposition tagging fields. All optional: a caller that only wants to
        ! advect particles need not synthesize climate forcing. Any field left
        ! out is recorded as MV in trc%dep.
        !   t2m, pr - monthly grids (nx,ny,nmon); the annual means and the
        !             precip-weighted temperature are derived at deposition.
        !   d18O    - monthly grid (nx,ny,nmon). If only the annual d18O_ann
        !             (nx,ny) is supplied instead, all months are set to it.
        real(wp), intent(IN), optional :: lon(:,:), lat(:,:)
        real(wp), intent(IN), optional :: t2m(:,:,:), pr(:,:,:), d18O(:,:,:)
        real(wp), intent(IN), optional :: d18O_ann(:,:)

        logical, intent(IN) :: dep_now, stats_now
        character(len=*), intent(IN), optional :: order 
        real(wp), intent(IN), optional :: sigma_srf     ! Value at surface by default (1 or 0?)

        ! Local variables  
        character(len=3) :: idx_order
        integer    :: i, j, k, nx, ny, nz, nz_uz
        logical    :: rev_z, rev_z_uz, has_srf
        real(wp) :: sigma_srf_val
        real(wp), allocatable :: x1(:), y1(:), z1(:)
        real(wp), allocatable :: x_ux1(:), y_uy1(:), z_uz1(:)
        real(wp), allocatable :: z_srf1(:,:), H1(:,:)
        real(wp), allocatable :: ux1(:,:,:), uy1(:,:,:), uz1(:,:,:)
        real(wp), allocatable :: usig1(:,:,:)
        real(wp), allocatable :: ux_srf_aa(:,:), uy_srf_aa(:,:)
        real(wp), allocatable :: lon1(:,:), lat1(:,:)
        real(wp), allocatable :: t2m1(:,:,:), pr1(:,:,:), d18O1(:,:,:)   ! monthly (nx,ny,nmon)
        real(wp) :: pr_sum_pt                                            ! sum of monthly precip at a point
        real(wp) :: ux0, uy0, uz0
        real(wp) :: sigp                        ! Normalized sigma of a particle (z relative to its column)
        type(lin3_interp_par_type) :: par_pt      ! Thread-local bilinear weights (kept off the shared par_lin for OpenMP)
        real(wp) :: dt

        ! Update the transient parameters
        if (trc%par%use_par_trans) then 
            call tracer_par_update(trc%par,trc%par%tpar,time)
        end if 

        ! Update current time and time step
        trc%now%time_old = trc%now%time 
        trc%now%time     = time 
        trc%now%dt       = real(dble(trc%now%time) - dble(trc%now%time_old),prec_time)

        ! Update record of last deposition time if dep_now
        if (dep_now) trc%now%time_dep = trc%now%time 

        ! Determine order of indices (default ijk)
        idx_order = "ijk"
        if (present(order)) idx_order = trim(order)

        ! Surface-orientation flag (sigma_srf) as a plain value, so the axis
        ! preparation can be shared by the aa z-axis and the uz interface axis.
        has_srf = present(sigma_srf)
        sigma_srf_val = 0.0
        if (has_srf) sigma_srf_val = sigma_srf

        ! Note: GRISLI (nx,ny,nz): sigma goes from 1 to 0, so sigma(1)=1 [surface], sigma(nz)=0 [base]
        !       SICO (nz,ny,nx): sigma(1) = 0, sigma(nz) = 1
        ! reshape routines ensure ascending z-axis (nx,ny,nz) with sigma(nz)=1 [surface]

        ! Prepare the aa z-axis (used by ux, uy, and by uz when z_uz is absent).
        call prepare_sigma_axis(z,trc%par%is_sigma,has_srf,sigma_srf_val,z1,rev_z)

        ! Prepare the uz vertical axis. With z_uz present, uz is sampled on its
        ! native interface levels (zeta_ac); absent, it shares the aa z-axis.
        if (present(z_uz)) then
            call prepare_sigma_axis(z_uz,trc%par%is_sigma,has_srf,sigma_srf_val,z_uz1,rev_z_uz)
        else
            z_uz1    = z1
            rev_z_uz = rev_z
        end if

        call tracer_reshape1D_vec(x, x1,rev=.FALSE.)
        call tracer_reshape1D_vec(y, y1,rev=.FALSE.)

        ! Horizontal staggered axes for ux (acx) and uy (acy). Absent => aa.
        if (present(x_ux)) then
            call tracer_reshape1D_vec(x_ux, x_ux1,rev=.FALSE.)
        else
            x_ux1 = x1
        end if
        if (present(y_uy)) then
            call tracer_reshape1D_vec(y_uy, y_uy1,rev=.FALSE.)
        else
            y_uy1 = y1
        end if

        call tracer_reshape2D_field(idx_order,z_srf,z_srf1)
        call tracer_reshape2D_field(idx_order,H,H1)

        ! Each component is reshaped on its own vertical orientation: ux/uy on
        ! the aa z-axis, uz on the (possibly distinct) interface axis.
        call tracer_reshape3D_field(idx_order,ux,ux1,rev_z=rev_z)
        call tracer_reshape3D_field(idx_order,uy,uy1,rev_z=rev_z)
        call tracer_reshape3D_field(idx_order,uz,uz1,rev_z=rev_z_uz)
        
        if (dep_now) then
            ! Also reshape deposition fields. H1 is already reshaped, so it
            ! gives the (nx,ny) shape to fall back on for an absent field.
            call reshape2D_optional(idx_order,lon1,size(H1,1),size(H1,2),lon)
            call reshape2D_optional(idx_order,lat1,size(H1,1),size(H1,2),lat)

            ! Monthly climate grids. An absent field becomes an all-MV (nx,ny,nmon)
            ! array so the per-point deposition code has a uniform shape to sample.
            call reshape_monthly_optional(idx_order,t2m1, size(H1,1),size(H1,2),t2m)
            call reshape_monthly_optional(idx_order,pr1,  size(H1,1),size(H1,2),pr)

            ! d18O: prefer the monthly grid; else broadcast the annual grid across
            ! all months; else leave all-MV.
            if (present(d18O)) then
                call reshape_monthly_optional(idx_order,d18O1,size(H1,1),size(H1,2),d18O)
            else if (present(d18O_ann)) then
                call reshape_annual_to_monthly(idx_order,d18O1,size(H1,1),size(H1,2),d18O_ann)
            else
                call reshape_monthly_optional(idx_order,d18O1,size(H1,1),size(H1,2))
            end if

        end if

        ! Get axis sizes (par%is_profile marks a 2D domain with a ghost y-axis)
        nx    = size(x1,1)
        ny    = size(y1,1)
        nz    = size(z1,1)
        nz_uz = size(z_uz1,1)

        if (trim(trc%par%interp_method) .eq. "spline") then

            ! Allocate z-velocity field in sigma coordinates. uz lives on its
            ! own vertical axis, so use nz_uz here rather than nz.
            if (allocated(usig1)) deallocate(usig1)
            allocate(usig1(nx,ny,nz_uz))

            usig1 = 0.0
            do k = 1, nz_uz
                where (H1 .gt. 0.0) usig1(:,:,k) = uz1(:,:,k) / H1
            end do

            ! Each spline is built on its component's native axes.
            call interp_bspline3D_weights(bspline3d_ux,x_ux1,y1,z1,ux1)
            call interp_bspline3D_weights(bspline3d_uy,x1,y_uy1,z1,uy1)
            call interp_bspline3D_weights(bspline3d_uz,x1,y1,z_uz1,usig1)

            write(*,*) "spline weights calculated."
        end if

        ! Interpolate to get the right elevation and velocity at each particle.
        ! The particles are independent (each writes only its own index i), so
        ! this loop is embarrassingly parallel. Weights (par_pt) are thread-local;
        ! the velocity splines are firstprivate so each thread evaluates its own
        ! copy — the bspline evaluate mutates an internal search cache, which
        ! would otherwise race. The pragmas are inert unless built with openmp=1.
        !$omp parallel do default(shared) schedule(dynamic,64) &
        !$omp   private(i,par_pt,ux0,uy0,uz0,sigp) &
        !$omp   firstprivate(bspline3d_ux,bspline3d_uy,bspline3d_uz)
        do i = 1, trc%par%n

            if (trc%now%active(i) .eq. 2) then

                ! Temporarily store velocity of this time step (for accelaration calculation)
                ux0 = trc%now%ux(i)
                uy0 = trc%now%uy(i)
                uz0 = trc%now%uz(i)

                ! Linear interpolation used for surface position
                par_pt = interp_bilinear_weights(x1,y1,xout=trc%now%x(i),yout=trc%now%y(i))
                trc%now%H(i)     = interp_bilinear(par_pt,H1)
                trc%now%z_srf(i) = interp_bilinear(par_pt,z_srf1)
                trc%now%z(i)     = trc%now%z_srf(i) - trc%now%dpth(i)

                if (trim(trc%par%interp_method) .eq. "linear") then
                    ! Trilinear interpolation, each component on its own axes.
                    ! interp_vel_linear rebuilds the cartesian z-levels from the
                    ! component's sigma axis, so ux/uy/uz are sampled directly on
                    ! their native (acx / acy / zeta_ac) locations.
                    trc%now%ux(i) = interp_vel_linear(x_ux1,y1,z1,ux1, &
                        trc%now%x(i),trc%now%y(i),trc%now%z(i),trc%now%z_srf(i),trc%now%H(i))
                    trc%now%uy(i) = interp_vel_linear(x1,y_uy1,z1,uy1, &
                        trc%now%x(i),trc%now%y(i),trc%now%z(i),trc%now%z_srf(i),trc%now%H(i))
                    trc%now%uz(i) = interp_vel_linear(x1,y1,z_uz1,uz1, &
                        trc%now%x(i),trc%now%y(i),trc%now%z(i),trc%now%z_srf(i),trc%now%H(i))

                else
                    ! Spline interpolation, each spline built on its own axes. The
                    ! vertical coordinate is the particle's normalized sigma
                    ! (z relative to its own column base), not z/H — the latter
                    ! left the [0,1] knot range wherever the bed is not at z=0.
                    ! The horizontal clamp mirrors the linear path so the acx/acy
                    ! half-cell edge strip samples the nearest face.
                    sigp = (trc%now%z(i) - (trc%now%z_srf(i)-trc%now%H(i))) / trc%now%H(i)
                    trc%now%ux(i) = interp_vel_spline(bspline3d_ux,x_ux1,y1, &
                        trc%now%x(i),trc%now%y(i),sigp)
                    trc%now%uy(i) = interp_vel_spline(bspline3d_uy,x1,y_uy1, &
                        trc%now%x(i),trc%now%y(i),sigp)
                    trc%now%uz(i) = interp_vel_spline(bspline3d_uz,x1,y1, &
                        trc%now%x(i),trc%now%y(i),sigp) *trc%now%H(i)  ! sigma => m

                end if

                ! Update acceleration term
                trc%now%ax(i) = (trc%now%ux(i) - ux0) / trc%now%dt
                trc%now%ay(i) = (trc%now%uy(i) - uy0) / trc%now%dt
                trc%now%az(i) = (trc%now%uz(i) - uz0) / trc%now%dt

                ! Filler values of the tracer state, in the future these should
                ! equal the surface temperature and the accumulation rate at the time of
                ! deposition and be calculated otherwise
                trc%now%T(i)   = 260.0
                trc%now%thk(i) = 0.3

            end if

        end do
        !$omp end parallel do

        ! Update the tracer thickness, then destroy points that are too thin
        ! == TO DO ==
        ! Unimplemented. The `thk_min` namelist parameter that would control it
        ! was removed rather than left as a knob that silently does nothing;
        ! restore it here when this is written.

        ! Update the tracer positions
        call calc_position(trc%now%x,trc%now%y,trc%now%z,trc%now%ux,trc%now%uy,trc%now%uz, &
                           trc%now%ax,trc%now%ay,trc%now%az,trc%now%dt,trc%now%active)

!         call calc_position(trc%now%x,trc%now%y,trc%now%dpth,trc%now%ux,trc%now%uy,-trc%now%uz,trc%now%dt,trc%now%active)

        ! Depth below the surface, for active points only. An inactive point has
        ! z_srf == z == MV, whose difference would otherwise clamp to a depth of
        ! zero and read as a particle sitting at the surface.
        where (trc%now%active .ne. 0)
            trc%now%dpth = max(trc%now%z_srf - trc%now%z, 0.0)
        elsewhere
            trc%now%dpth = MV
        end where

        ! Destroy points that moved outside the valid region
        call tracer_deactivate(trc,x1,y1,maxval(H1))

        ! Clone eligible tracers into freed slots, before deposition claims them,
        ! so surviving old tracers take priority over fresh surface deposition.
        call tracer_clone(trc%par,trc%now,trc%dep,x1,y1)

        ! Activate new tracers if desired. The deposition weighting needs ux and
        ! uy co-located, so in the native staggered case destagger the surface
        ! slices to aa-nodes first (identity when a component is already aa).
        if (dep_now) then
            if (present(x_ux) .or. present(y_uy)) then
                ux_srf_aa = destagger_to_aa(x_ux1,y1,ux1(:,:,nz),x1,y1)
                uy_srf_aa = destagger_to_aa(x1,y_uy1,uy1(:,:,nz),x1,y1)
            else
                ux_srf_aa = ux1(:,:,nz)
                uy_srf_aa = uy1(:,:,nz)
            end if

            call tracer_activate(trc%par,trc%now,x1,y1,H=H1,lat=lat1, &
                                ux_srf=ux_srf_aa,uy_srf=uy_srf_aa,nmax=trc%par%n_max_dep)
        end if

        ! Finish activation for necessary points 
        do i = 1, trc%par%n 

            if (trc%now%active(i) .eq. 1) then 
                ! Point became active now, further initializations needed below

                ! Point is at the surface, so only bilinear interpolation is
                ! needed. z_srf and H (aa fields) use aa weights; each velocity
                ! component uses weights on its own horizontal axis and its own
                ! top sigma level (the surface).
                par_lin = interp_bilinear_weights(x1,y1,xout=trc%now%x(i),yout=trc%now%y(i))

                ! Apply interpolation weights to variables
                trc%now%dpth(i)  = 0.01   ! Always deposit just below the surface (eg 1 cm) to avoid zero z-velocity
                trc%now%z_srf(i) = interp_bilinear(par_lin,z_srf1)
                trc%now%z(i)     = trc%now%z_srf(i)-trc%now%dpth(i)

                trc%now%H(i)     = interp_bilinear(par_lin,H1)
                trc%now%ux(i)    = interp_bilinear(interp_bilinear_weights(x_ux1,y1,xout=trc%now%x(i),yout=trc%now%y(i)),ux1(:,:,nz))
                trc%now%uy(i)    = interp_bilinear(interp_bilinear_weights(x1,y_uy1,xout=trc%now%x(i),yout=trc%now%y(i)),uy1(:,:,nz))
                trc%now%uz(i)    = interp_bilinear(par_lin,uz1(:,:,nz_uz))
                trc%now%ax(i)    = 0.0 
                trc%now%ay(i)    = 0.0 
                trc%now%az(i)    = 0.0
                
                ! Initialize state variables
                trc%now%T(i)   = 260.0 
                trc%now%thk(i) = 0.3 

                ! Define deposition values 
                trc%dep%time(i) = trc%now%time 
                trc%dep%H(i)    = trc%now%H(i)
                trc%dep%x(i)    = trc%now%x(i)
                trc%dep%y(i)    = trc%now%y(i)
                trc%dep%z(i)    = trc%now%z(i) 

                ! An absent tagging field is recorded as MV rather than
                ! interpolated from a placeholder array.
                trc%dep%lon(i) = MV
                trc%dep%lat(i) = MV
                if (present(lon)) trc%dep%lon(i) = interp_bilinear(par_lin,lon1)
                if (present(lat)) trc%dep%lat(i) = interp_bilinear(par_lin,lat1)

                ! Monthly climate tags interpolated to the deposition point, and
                ! the annual quantities derived from them. An absent input arrived
                ! here as an all-MV monthly grid, so each block is guarded on the
                ! field's presence and records MV when it was not supplied.
                trc%dep%t2m(i,:)      = MV
                trc%dep%pr(i,:)       = MV
                trc%dep%d18O(i,:)     = MV
                trc%dep%t2m_ann(i)    = MV
                trc%dep%pr_ann(i)     = MV
                trc%dep%t2m_prann(i)  = MV
                trc%dep%d18O_ann(i)   = MV
                trc%dep%d18O_prann(i) = MV

                if (present(t2m)) then
                    do k = 1, nmon
                        trc%dep%t2m(i,k) = interp_bilinear(par_lin,t2m1(:,:,k))
                    end do
                    trc%dep%t2m_ann(i) = sum(trc%dep%t2m(i,:))/real(nmon,wp)
                end if

                if (present(pr)) then
                    do k = 1, nmon
                        trc%dep%pr(i,k) = interp_bilinear(par_lin,pr1(:,:,k))
                    end do
                    trc%dep%pr_ann(i) = sum(trc%dep%pr(i,:))/real(nmon,wp)
                end if

                if (present(d18O) .or. present(d18O_ann)) then
                    do k = 1, nmon
                        trc%dep%d18O(i,k) = interp_bilinear(par_lin,d18O1(:,:,k))
                    end do
                    trc%dep%d18O_ann(i) = sum(trc%dep%d18O(i,:))/real(nmon,wp)
                end if

                ! Precip-weighted annual means (temperature, d18O), formed only
                ! when precip is available and its annual sum is positive.
                if (present(pr)) then
                    pr_sum_pt = sum(trc%dep%pr(i,:))
                    if (pr_sum_pt .gt. 0.0_wp) then
                        if (present(t2m)) &
                            trc%dep%t2m_prann(i)  = sum(trc%dep%t2m(i,:) *trc%dep%pr(i,:))/pr_sum_pt
                        if (present(d18O) .or. present(d18O_ann)) &
                            trc%dep%d18O_prann(i) = sum(trc%dep%d18O(i,:)*trc%dep%pr(i,:))/pr_sum_pt
                    end if
                end if

                trc%now%active(i) = 2

            end if 

        end do 


        ! == TO DO == 
        ! - Attach whatever information we want to trace (age, deposition elevation and location, climate, isotopes, etc)
        ! - Potentially attach this information via a separate subroutine, using a flag to see if
        !   it was just deposited, then in main program calling eg, tracer_add_dep_variable(trc,"T",T),
        !   where the argument "T" should match a list of available variables, and T should be the variable
        !   to be stored from the main program. 
        !   Downside to above approach, is re-calculating the par_lin object every time. 

        ! Update summary statistics
        trc%par%n_active = count(trc%now%active.gt.0)

        ! Gridded (Eulerian) statistics on request. The stats object is only
        ! allocated when par%stats; calc_tracer_stats (tracer_stats) takes the
        ! tracer state as plain arrays, so this module does not depend on it.
        if (trc%par%stats .and. stats_now) then
            call calc_tracer_stats(trc%stats, trc%now%active, &
                                   trc%now%x, trc%now%y, trc%now%dpth, trc%now%H, &
                                   trc%dep%time, trc%dep%z, trc%dep%lon, trc%dep%lat, &
                                   trc%dep%t2m_ann, trc%dep%pr_ann, trc%dep%d18O_ann)
        end if

        return

    end subroutine tracer_update

    subroutine tracer_end(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc

        ! Allocate the state variables
        call tracer_deallocate(trc%now,trc%dep)
        call tracer_stats_end(trc%stats)

        write(*,*) "tracer:: tracer object deallocated."
        
        return 

    end subroutine tracer_end

    ! ================================================
    !
    ! tracer management 
    !
    ! ================================================
    
    subroutine tracer_activate(par,now,x,y,H,lat,ux_srf,uy_srf,nmax)
        ! Use this to activate individual or multiple tracers (not more than nmax)
        ! Only determine x/y position here, later interpolate z_srf and deposition
        ! information 

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(wp), intent(IN) :: x(:), y(:)
        real(wp), intent(IN) :: H(:,:), lat(:,:), ux_srf(:,:), uy_srf(:,:) 
        integer, intent(IN) :: nmax  

        integer :: ntot  
        real(wp) :: p(size(H,1),size(H,2)), p_init(size(H,1),size(H,2))
        integer :: i, j, k, ij(2)
        real(wp), allocatable :: jit(:,:)
        real(wp) :: xmin, ymin, xmax, ymax 

        ! How many points can be activated?
        ntot = min(nmax,count(now%active == 0))


        if (ntot .gt. 0) then 
            ! Proceed with activation, since points are available 

            ! Determine initial desired distribution of points on low resolution grid.
            ! par%weight is validated in tracer_par_load.
            select case(trim(par%weight))

                case("vel")
                    ! Weight by surface velocity: slow ice is the likeliest site
                    p_init = gen_distribution_vel(uv=sqrt(ux_srf**2+uy_srf**2),H=H, &
                                                  uv_max=par%U_max_dep,H_min=par%H_min_dep)

                case DEFAULT
                    ! Weight by ice thickness ("linear", "quadratic" or "rand")
                    p_init = gen_distribution_thickness(H,H_min=par%H_min_dep, &
                                                        alpha=par%alpha,dist=par%weight)

            end select

            ! A profile domain repeats every field along its ghost y-axis, so
            ! each column would otherwise be selected once per ghost node,
            ! spending n_max_dep on redundant tracers and shrinking the range of
            ! x actually reached. Deposit into the first node only, which makes
            ! deposition independent of the ghost axis length.
            if (par%is_profile) p_init(:,2:) = 0.0

!             ! Additionally adjust distribution according to latitude
!             where (lat .lt. 70.0)
!                 p_init = 0.0
!             end where

            ! == TO DO ==
            ! Cap the density of existing particles near the surface, so that
            ! p_init is suppressed where too many tracers already sit within
            ! some distance of the surface. The `dens_z_lim` / `dens_max`
            ! namelist parameters that would control it were removed rather
            ! than left as knobs that silently do nothing; restore them here
            ! when this is written.

            ! gen_distribution_vel returns an already-normalized field, or all
            ! zeros when no cell meets the deposition criteria (no ice above
            ! H_min_dep, or all ice faster than U_max_dep). Renormalizing here
            ! would divide by a zero sum, giving NaN and depositing nothing.
            if (sum(p_init) .eq. 0.0) return

            ! Generate random numbers to populate points 
            allocate(jit(2,ntot))

            if (par%noise) then
                call random_number(jit)
                jit = (jit - 0.5)
                jit(1,:) = jit(1,:)*(x(2)-x(1))

                ! The profile domain's y-axis is a ghost: every field is
                ! constant along it, so jittering there is meaningless.
                if (par%is_profile) then
                    jit(2,:) = 0.0
                else
                    jit(2,:) = jit(2,:)*(y(2)-y(1))
                end if
            else
                jit = 0.0
            end if

    !         write(*,*) "range jit: ", minval(jit), maxval(jit)
    !         write(*,*) "npts: ", count(now%active == 0)
    !         write(*,*) "ntot: ", ntot 
    !         stop 
        
            ! Calculate domain boundaries to be able to apply limits 
            xmin = minval(x) 
            xmax = maxval(x) 
            ymin = minval(y) 
            ymax = maxval(y) 

            if (maxval(p_init) .gt. 0.0) then 
                ! Activate points in locations with non-zero probability
                ! This if-statement ensures some valid points currently exist in the domain
                
                k = 0 
                p = p_init   ! Set probability distribution to initial distribution 

                do j = 1, par%n 

                    if (now%active(j)==0) then 

                        now%active(j) = 1
                        k = k + 1
                        par%id_max = par%id_max+1
                        now%id(j)  = par%id_max
                        now%parent(j)   = 0    ! surface-deposited original, not a clone
                        now%n_cloned(j) = 0

                        ij = maxloc(p,mask=p.gt.0.0)
                        now%x(j) = x(ij(1)) + jit(1,k)
                        now%y(j) = y(ij(2)) + jit(2,k)
                        
                        if (now%x(j) .lt. xmin) now%x(j) = xmin 
                        if (now%x(j) .gt. xmax) now%x(j) = xmax 
                        if (now%y(j) .lt. ymin) now%y(j) = ymin 
                        if (now%y(j) .gt. ymax) now%y(j) = ymax 
                        
                        p(ij(1),ij(2)) = 0.0 
                         
                    end if 

                    ! Stop when all points have been allocated
                    if (k .ge. ntot) exit 

                    ! If there are no more points with non-zero probability, reset probability
                    if (maxval(p) .eq. 0.0) p = p_init 

                end do 
              
            end if

            ! Summary 
    !         write(*,*) "tracer_activate:: ", count(now%active == 0), count(now%active .eq. 1), count(now%active .eq. 2)

        end if 

        return 

    end subroutine tracer_activate 

    subroutine tracer_deactivate(trc,x,y,Hmax)
        ! Use this to deactivate individual or multiple tracers
        implicit none 

        type(tracer_class),   intent(INOUT) :: trc  
        real(wp), intent(IN) :: x(:), y(:) 
        real(wp), intent(IN) :: Hmax 

        ! Deactivate points where:
        !  - Thickness of ice sheet at point's location is below threshold
        !  - Point is above maximum ice thickness Hmax (interp error)
        !  - Point is past maximum depth into the ice sheet 
        !  - Velocity of point is higher than maximum threshold 
        !  - x/y position is out of boundaries of the domain 
        where (trc%now%active .gt. 0 .and. &
              ( trc%now%H .lt. trc%par%H_min                            .or. &
                trc%now%H .gt. Hmax                                     .or. &
                trc%now%dpth/trc%now%H .ge. trc%par%depth_max           .or. &
                sqrt(trc%now%ux**2 + trc%now%uy**2) .gt. trc%par%U_max  .or. &
                trc%now%x .lt. minval(x) .or. trc%now%x .gt. maxval(x)  .or. &
                trc%now%y .lt. minval(y) .or. trc%now%y .gt. maxval(y) ) ) 

            trc%now%active    = 0

            trc%now%id        = mv
            trc%now%parent    = mv
            trc%now%n_cloned  = mv
            trc%now%x         = mv
            trc%now%y         = mv 
            trc%now%z         = mv 
            trc%now%sigma     = mv 
            trc%now%z_srf     = mv 
            trc%now%dpth      = mv 
            trc%now%ux        = mv 
            trc%now%uy        = mv 
            trc%now%uz        = mv 
            trc%now%thk       = mv 
            trc%now%T         = mv 
            trc%now%H         = mv 

            trc%dep%time      = mv 
            trc%dep%H         = mv 
            trc%dep%x         = mv 
            trc%dep%y         = mv 
            trc%dep%z         = mv 

        end where 

        return

    end subroutine tracer_deactivate

    subroutine tracer_clone(par,now,dep,x,y)
        ! Spawn offset copies of tracers that are deep, slow-moving and within a
        ! deposition-time window, to raise the chance that old deposition times
        ! survive to present day. A clone inherits the full state and deposition
        ! record of its parent, offset by a random horizontal displacement and a
        ! random upward vertical displacement. This injects diffusion into an
        ! otherwise deterministic advection.
        implicit none

        type(tracer_par_class),   intent(INOUT) :: par
        type(tracer_state_class), intent(INOUT) :: now
        type(tracer_dep_class),   intent(INOUT) :: dep
        real(wp), intent(IN) :: x(:), y(:)

        integer    :: i, j, c
        real(wp) :: xmin, xmax, ymin, ymax
        real(wp) :: off(3)

        if (.not. par%clone) return

        xmin = minval(x); xmax = maxval(x)
        ymin = minval(y); ymax = maxval(y)

        ! Only originals (parent == 0) that have not yet cloned are eligible.
        ! Clones carry parent > 0, so they are excluded here and clones-of-clones
        ! never occur. Freshly spawned clones are written into slots j and given
        ! parent > 0, so when the outer loop later reaches slot j it skips it.
        do i = 1, par%n

            if ( now%active(i)   .eq. 2 .and. &
                 now%parent(i)   .eq. 0 .and. &
                 now%n_cloned(i) .eq. 0 ) then

                ! Eligibility: deep enough, slow enough, within the dep-time window
                if ( now%H(i) .gt. 0.0                                        .and. &
                     now%dpth(i)/now%H(i) .ge. par%clone_depth_frac           .and. &
                     sqrt(now%ux(i)**2 + now%uy(i)**2) .lt. par%clone_U_max    .and. &
                     dep%time(i) .ge. par%clone_dep_time_min                   .and. &
                     dep%time(i) .le. par%clone_dep_time_max ) then

                    ! Fill up to n_clones free slots with offset copies of tracer i
                    c = 0
                    do j = 1, par%n

                        if (now%active(j) .eq. 0) then

                            call random_number(off)   ! off in [0,1)

                            ! Copy full state from the parent
                            now%active(j)   = 2
                            par%id_max      = par%id_max + 1
                            now%id(j)       = par%id_max
                            now%parent(j)   = now%id(i)
                            now%n_cloned(j) = 0

                            now%sigma(j) = now%sigma(i)
                            now%z_srf(j) = now%z_srf(i)
                            now%ux(j)    = now%ux(i)
                            now%uy(j)    = now%uy(i)
                            now%uz(j)    = now%uz(i)
                            now%ax(j)    = now%ax(i)
                            now%ay(j)    = now%ay(i)
                            now%az(j)    = now%az(i)
                            now%thk(j)   = now%thk(i)
                            now%T(j)     = now%T(i)
                            now%H(j)     = now%H(i)

                            ! Offset position: horizontal uniform in [-1,1]*offset,
                            ! vertical uniform upward in [0,1]*offset (raises z).
                            now%x(j) = now%x(i) + (2.0*off(1)-1.0)*par%clone_offset_xy
                            now%z(j) = now%z(i) + off(3)*par%clone_offset_z

                            if (now%x(j) .lt. xmin) now%x(j) = xmin
                            if (now%x(j) .gt. xmax) now%x(j) = xmax

                            ! The profile domain's y-axis is a ghost (all fields
                            ! constant along it), so no meaningful offset there.
                            if (par%is_profile) then
                                now%y(j) = now%y(i)
                            else
                                now%y(j) = now%y(i) + (2.0*off(2)-1.0)*par%clone_offset_xy
                                if (now%y(j) .lt. ymin) now%y(j) = ymin
                                if (now%y(j) .gt. ymax) now%y(j) = ymax
                            end if

                            ! Depth follows from the offset z; z_srf/H/velocities
                            ! re-interpolate at the clone's own location next step.
                            now%dpth(j) = max(now%z_srf(j) - now%z(j), 0.0)

                            ! Copy the deposition record: a clone shares its
                            ! parent's age and provenance exactly.
                            dep%time(j)      = dep%time(i)
                            dep%H(j)         = dep%H(i)
                            dep%x(j)         = dep%x(i)
                            dep%y(j)         = dep%y(i)
                            dep%z(j)         = dep%z(i)
                            dep%lon(j)        = dep%lon(i)
                            dep%lat(j)        = dep%lat(i)
                            dep%t2m(j,:)      = dep%t2m(i,:)
                            dep%pr(j,:)       = dep%pr(i,:)
                            dep%d18O(j,:)     = dep%d18O(i,:)
                            dep%t2m_ann(j)    = dep%t2m_ann(i)
                            dep%pr_ann(j)     = dep%pr_ann(i)
                            dep%t2m_prann(j)  = dep%t2m_prann(i)
                            dep%d18O_ann(j)   = dep%d18O_ann(i)
                            dep%d18O_prann(j) = dep%d18O_prann(i)

                            c = c + 1

                        end if

                        if (c .ge. par%n_clones) exit

                    end do

                    ! Record how many clones were made. If the array was full
                    ! (c == 0) the tracer stays eligible and retries next step.
                    now%n_cloned(i) = c

                end if

            end if

        end do

        return

    end subroutine tracer_clone

    ! ================================================
    !
    ! tracer physics / stats
    !
    ! ================================================
    
    elemental subroutine calc_position(x,y,z,ux,uy,uz,ax,ay,az,dt,active)

        implicit none 

        real(wp),   intent(INOUT) :: x, y, z 
        real(wp),   intent(IN)    :: ux, uy, uz 
        real(wp),   intent(IN)    :: ax, ay, az 
        real(prec_time), intent(IN)    :: dt 
        integer,         intent(IN)    :: active 

        if (active .gt. 0) then 

            x = x + ux*dt + 0.5*ax*dt**2 
            y = y + uy*dt + 0.5*ay*dt**2 
            z = z + uz*dt + 0.5*az*dt**2
            
        end if

        return

    end subroutine calc_position

    subroutine prepare_sigma_axis(z,is_sigma,has_srf,sigma_srf,z1,rev)
        ! Turn a raw vertical axis into an ascending axis (z1) plus the flag
        ! saying whether the input was reversed. Shared by the aa z-axis and the
        ! uz interface axis (zeta_ac) so both get identical sigma bounding and
        ! surface-orientation handling.

        implicit none

        real(wp), intent(IN)  :: z(:)
        logical,    intent(IN)  :: is_sigma, has_srf
        real(wp), intent(IN)  :: sigma_srf
        real(wp), intent(INOUT), allocatable :: z1(:)
        logical,    intent(OUT) :: rev

        real(wp), allocatable :: za(:)

        allocate(za(size(z)))
        za = z

        if (is_sigma) then
            ! Ensure z-axis is properly bounded
            where (abs(za) .lt. 1e-5) za = 0.0

            if (minval(za) .lt. 0.0 .or. maxval(za) .gt. 1.0) then
                write(0,*) "tracer:: error: sigma axis not bounded between zero and one."
                write(0,*) "z = ", za
                error stop
            end if

            ! Correct the sigma values if necessary, so that
            ! sigma==0 [base]; sigma==1 [surface]
            if (has_srf .and. sigma_srf .eq. 0.0) za = 1.0 - za
        end if

        ! Determine whether z-axis is initially ascending or descending, then
        ! reshape to ascending order.
        rev = (za(1) .gt. za(size(za)))
        call tracer_reshape1D_vec(za,z1,rev=rev)

        return

    end subroutine prepare_sigma_axis

    function interp_vel_linear(x,y,sigma,field,xp,yp,zp,z_srf_p,H_p) result(u)
        ! Trilinear interpolation of one velocity component onto a particle
        ! point. sigma is that component's own sigma axis, so the cartesian
        ! z-levels (zc) are rebuilt from it here — this is what lets ux/uy/uz be
        ! sampled on their native staggered axes without an aa-node detour.

        implicit none

        real(wp), intent(IN) :: x(:), y(:), sigma(:)
        real(wp), intent(IN) :: field(:,:,:)
        real(wp), intent(IN) :: xp, yp, zp, z_srf_p, H_p
        real(wp) :: u

        real(wp) :: zc(size(sigma))
        real(wp) :: xo, yo
        type(lin3_interp_par_type) :: par

        ! A staggered horizontal axis (acx/acy) is offset half a cell from the
        ! aa nodes, so it leaves a half-cell strip at one edge uncovered. Clamp
        ! the query point into the axis range there, giving the nearest edge
        ! velocity rather than MV. This is a no-op for an aa axis (particles are
        ! kept within the aa domain), so it does not change the aa path.
        xo = min(max(xp,min(x(1),x(size(x)))),max(x(1),x(size(x))))
        yo = min(max(yp,min(y(1),y(size(y)))),max(y(1),y(size(y))))

        ! Cartesian z of the component's sigma levels at this column
        ! (equivalent to (z_srf - depth) = z_srf - (1-sigma)*H).
        zc  = (z_srf_p - H_p) + sigma*H_p
        par = interp_trilinear_weights(x,y,zc,xout=xo,yout=yo,zout=zp)
        u   = interp_trilinear(par,field)

        return

    end function interp_vel_linear

    function interp_vel_spline(bspl,x,y,xp,yp,sig) result(u)
        ! Evaluate a pre-built velocity spline at a particle point. x and y are
        ! that component's horizontal axes, used only to clamp the query into
        ! range so the acx/acy half-cell edge strip samples the nearest face
        ! instead of evaluating outside the spline (a no-op for an aa axis).

        implicit none

        type(bspline_3d), intent(INOUT) :: bspl
        real(wp),       intent(IN)    :: x(:), y(:)
        real(wp),       intent(IN)    :: xp, yp, sig
        real(wp) :: u

        real(wp) :: xo, yo

        xo = min(max(xp,min(x(1),x(size(x)))),max(x(1),x(size(x))))
        yo = min(max(yp,min(y(1),y(size(y)))),max(y(1),y(size(y))))

        u = interp_bspline3D(bspl,xo,yo,sig)

        return

    end function interp_vel_spline

    function destagger_to_aa(x_src,y_src,field,x_aa,y_aa) result(field_aa)
        ! Interpolate a staggered 2D field onto the aa nodes (x_aa,y_aa). Used
        ! for the surface-velocity deposition heuristic, which needs ux and uy
        ! co-located. Out-of-range aa nodes (the half-cell overhang at the
        ! domain edge) clamp to the nearest source node rather than giving MV.
        ! Identity when x_src/y_src already equal x_aa/y_aa (aa input).

        implicit none

        real(wp), intent(IN) :: x_src(:), y_src(:)
        real(wp), intent(IN) :: field(:,:)
        real(wp), intent(IN) :: x_aa(:), y_aa(:)
        real(wp) :: field_aa(size(x_aa),size(y_aa))

        integer    :: i, j
        real(wp) :: xo, yo, xmn, xmx, ymn, ymx
        type(lin3_interp_par_type) :: par

        xmn = minval(x_src); xmx = maxval(x_src)
        ymn = minval(y_src); ymx = maxval(y_src)

        do j = 1, size(y_aa)
        do i = 1, size(x_aa)
            xo  = min(max(x_aa(i),xmn),xmx)
            yo  = min(max(y_aa(j),ymn),ymx)
            par = interp_bilinear_weights(x_src,y_src,xout=xo,yout=yo)
            field_aa(i,j) = interp_bilinear(par,field)
        end do
        end do

        return

    end function destagger_to_aa

    function gen_distribution_vel(uv,H,uv_max,H_min) result(p)

        implicit none 

        real(wp), intent(IN) :: uv(:,:), H(:,:)
        real(wp), intent(IN) :: uv_max, H_min 
        real(wp) :: p(size(H,1),size(H,2))

        ! Local variables
        real(wp) :: p_sum 

        
        p = 0.0
        ! Note: uv==0 must be included. Slow ice is the most likely deposition
        ! site (p -> 1), and an ice divide has uv==0 exactly by symmetry.
        where (uv .ge. 0.0 .and. uv .lt. uv_max .and. H .gt. H_min)
            p = 1.0 - uv/uv_max
        end where

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution_vel

    function gen_distribution_thickness(H,H_min,alpha,dist) result(p)

        implicit none 

        real(wp), intent(IN) :: H(:,:)
        real(wp), intent(IN) :: H_min, alpha 
        character(len=*), intent(IN) :: dist 
        real(wp) :: p(size(H,1),size(H,2))

        ! Local variables
        integer    :: k, ij(2)
        real(wp) :: p_sum, H_range

        ! No cell is thick enough to deposit into. Return an all-zero field;
        ! the caller treats a zero sum as "no valid deposition sites". This
        ! also guards the division by H_range below.
        H_range = maxval(H) - H_min
        if (H_range .le. 0.0) then
            p = 0.0
            return
        end if

        select case(trim(dist))

            case("linear")

                p = (alpha * max(H-H_min,0.0) / H_range)

            case("quadratic")

                p = (alpha * max(H-H_min,0.0) / H_range)**2

            case("rand")

                ! Random even distribution (all points equally likely)
                call random_number(p)
                where (H .le. H_min) p = 0.0

            case DEFAULT

                write(0,*) "gen_distribution_thickness:: error: unknown dist: "//trim(dist)
                error stop


        end select 

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution_thickness

    function gen_distribution_direction(x,y,u,v,theta_max) result(p)

        implicit none 

        real(wp), intent(IN) :: x(:), y(:), u(:,:), v(:,:) 
        real(wp), intent(IN) :: theta_max 
        real(wp) :: p(size(u,1),size(u,2))

        ! Local variables
        integer    :: k, ij(2)
        real(wp) :: p_sum 

        p = 1.0 

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution_direction

    elemental function calc_angle(x1,y1,x2,y2) result(theta)
        ! Given a vector, calculate the angle wrt unit circle 

        implicit none 

        real(wp), intent(IN) :: x1, y1, x2, y2 
        real(wp) :: theta 

        theta = atan2((y2-y1),(x2-x1))

        return 

    end function calc_angle 

    ! ================================================
    !
    ! Initialization routines 
    !
    ! ================================================

    subroutine tracer_set_seed(seed)
        ! Seed the intrinsic random number generator, which drives the
        ! deposition location jitter (par%noise) and the "rand" weight. A
        ! positive seed makes a run reproducible; seed <= 0 defers to the
        ! processor's nondeterministic seeding.

        implicit none

        integer, intent(IN) :: seed

        ! Local variables
        integer :: n, i
        integer, allocatable :: sd(:)

        if (seed .le. 0) then
            call random_seed()
            return
        end if

        call random_seed(size=n)
        allocate(sd(n))

        ! Spread one user-facing seed across the generator's whole state. An
        ! all-equal state is a poor starting point for the xorshift generator
        ! gfortran uses.
        do i = 1, n
            sd(i) = seed + 37*(i-1)
        end do

        call random_seed(put=sd)

        return

    end subroutine tracer_set_seed

    subroutine tracer_par_load(par,filename,is_sigma)

        implicit none 

        type(tracer_par_class), intent(OUT) :: par 
        character(len=*),       intent(IN)  :: filename 
        logical, intent(IN) :: is_sigma 

        ! Local variables
        integer :: n_time_iso

        call nml_read(filename,"trc","dt",            par%dt)
        call nml_read(filename,"trc","n",             par%n)
        call nml_read(filename,"trc","n_max_dep",     par%n_max_dep)
        call nml_read(filename,"trc","dt_dep",        par%dt_dep)
        call nml_read(filename,"trc","dt_write",      par%dt_write)
        call nml_read(filename,"trc","H_min",         par%H_min)
        call nml_read(filename,"trc","depth_max",     par%depth_max)
        call nml_read(filename,"trc","U_max",         par%U_max)
        call nml_read(filename,"trc","U_max_dep",     par%U_max_dep)
        call nml_read(filename,"trc","H_min_dep",     par%H_min_dep)
        call nml_read(filename,"trc","alpha",         par%alpha)
        call nml_read(filename,"trc","weight",        par%weight)
        call nml_read(filename,"trc","noise",         par%noise)
        call nml_read(filename,"trc","seed",          par%seed)
        call nml_read(filename,"trc","interp_method", par%interp_method)
        call nml_read(filename,"trc","par_trans_file",par%par_trans_file)

        ! Cloning is opt-in. Its knobs are read only when enabled, so a run that
        ! does not clone need not carry the extra namelist entries.
        call nml_read(filename,"trc","clone",         par%clone)
        if (par%clone) then
            call nml_read(filename,"trc","n_clones",           par%n_clones)
            call nml_read(filename,"trc","clone_depth_frac",   par%clone_depth_frac)
            call nml_read(filename,"trc","clone_U_max",        par%clone_U_max)
            call nml_read(filename,"trc","clone_dep_time_min", par%clone_dep_time_min)
            call nml_read(filename,"trc","clone_dep_time_max", par%clone_dep_time_max)
            call nml_read(filename,"trc","clone_offset_xy",    par%clone_offset_xy)
            call nml_read(filename,"trc","clone_offset_z",     par%clone_offset_z)
        else
            par%n_clones = 0
        end if

        ! Gridded statistics are opt-in. Their parameters are read only when the
        ! stats flag is set, so a run that does not want them (e.g. the RH2003
        ! profile) need not carry the extra namelist entries.
        call nml_read(filename,"trc","stats",         par%stats)
        if (par%stats) then
            call nml_read(filename,"trc","dt_write_stats", par%dt_write_stats)
            call nml_read(filename,"trc","n_depth",        par%n_depth)
            call nml_read(filename,"trc","dt_iso",         par%dt_iso)
            call nml_read(filename,"trc","n_time_iso",     n_time_iso)
            if (allocated(par%time_iso)) deallocate(par%time_iso)
            allocate(par%time_iso(n_time_iso))
            call nml_read(filename,"trc","time_iso",       par%time_iso)
        end if

        ! Define additional parameter values
        par%is_sigma  = is_sigma
        par%n_active  = 0

        ! A 3D domain unless the caller says otherwise; tracer2D_init sets this.
        par%is_profile = .FALSE.

        par%use_par_trans = .FALSE.
        if (trim(par%par_trans_file) .ne. "None") then
            par%use_par_trans = .TRUE.

            call tracer_par_trans_load(par%tpar,par%par_trans_file)
        end if

        ! Consistency checks
        if (trim(par%interp_method) .ne. "linear" .and. &
            trim(par%interp_method) .ne. "spline" ) then
            write(0,*) "tracer_init:: error: interp_method must be 'linear' &
            &or 'spline': "//trim(par%interp_method)
            error stop
        end if

        if (trim(par%weight) .ne. "vel"       .and. &
            trim(par%weight) .ne. "linear"    .and. &
            trim(par%weight) .ne. "quadratic" .and. &
            trim(par%weight) .ne. "rand" ) then
            write(0,*) "tracer_init:: error: weight must be 'vel', 'linear', &
            &'quadratic' or 'rand': "//trim(par%weight)
            error stop
        end if

        if (par%clone) then
            if (par%n_clones .lt. 1) then
                write(0,*) "tracer_init:: error: n_clones must be >= 1 when clone = .TRUE."
                error stop
            end if
            if (par%clone_depth_frac .le. 0.0 .or. par%clone_depth_frac .gt. 1.0) then
                write(0,*) "tracer_init:: error: clone_depth_frac must be in (0,1]: ", par%clone_depth_frac
                error stop
            end if
        end if


        return

    end subroutine tracer_par_load

    function uniform_depth(n_depth) result(depth_norm)
        ! n_depth normalized-depth levels, evenly spaced over (0,1] with the
        ! deepest at 1. The shallowest is 1/n_depth, not 0: the surface band is
        ! open at the top (a tracer shallower than the first level bins into it).

        implicit none

        integer, intent(IN) :: n_depth
        real(prec_wrt) :: depth_norm(n_depth)
        integer :: i

        do i = 1, n_depth
            depth_norm(i) = real(i,prec_wrt) / real(n_depth,prec_wrt)
        end do

        return

    end function uniform_depth
    
    subroutine tracer_par_update(par,tpar,time)
        ! Update transient parameter values for current time 

        implicit none 

        type(tracer_par_class), intent(INOUT) :: par 
        type(tracer_par_trans_class), intent(IN) :: tpar 
        real(prec_time) :: time 

        ! Local variables 
        integer :: i, n, k  

        n = size(tpar%time)

        ! Initially assume the first row is correct 
        k = 1
        
        ! Check to see if one of following rows is correct, update k 
        do i = 2, n 
            if (tpar%time(i) .gt. time) exit 
            k = k+1
        end do 

        par%H_min_dep = tpar%H_min_dep(k) 
        par%dt_dep    = tpar%dt_dep(k)
        par%n_max_dep = tpar%n_max_dep(k) 
        par%dt_write  = tpar%dt_write(k) 

        return 

    end subroutine tracer_par_update

    subroutine tracer_par_trans_load(tpar,filename)
        ! This subroutine will read a time series of
        ! several columns [time,par1,par2,...,parN] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        type(tracer_par_trans_class), intent(OUT) :: tpar 
        character(len=*), intent(IN) :: filename 

        ! Local variables 
        integer, parameter :: f = 191
        integer, parameter :: nmax = 10000

        integer :: i, stat, n 
        character(len=256) :: str, str1 
        real(4) :: x(nmax), y1(nmax), y2(nmax), y3(nmax), y4(nmax)

        ! Open file for reading 
        open(f,file=filename,status="old")

        ! Read the header in the first line: 
        read(f,*,IOSTAT=stat) str

        n = 0 

        do i = 1, nmax 
            read(f,'(a100)',IOSTAT=stat) str 

            ! Exit loop if the end-of-file is reached 
            if(IS_IOSTAT_END(stat)) exit 

            str1 = adjustl(trim(str))
!            str1=str
            if ( len(trim(str1)) .gt. 0 ) then 
                if ( .not. (str1(1:1) == "!" .or. &
                            str1(1:1) == "#") ) then 
                    n = n+1
                    read(str1,*) x(n), y1(n), y2(n), y3(n), y4(n)
                end if
            end if  
        end do 

        ! Close the file
        close(f) 

        if (n .eq. nmax) then 
            write(*,*) "tracer_par_trans_load:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data 
        call tracer_par_trans_allocate(tpar,n)

        tpar%time      =  x(1:n) 
        tpar%H_min_dep = y1(1:n) 
        tpar%dt_dep    = y2(1:n) 
        tpar%n_max_dep = y3(1:n) 
        tpar%dt_write  = y4(1:n) 

        write(*,*) "tracer_par_trans_load:: Time series read from file: "//trim(filename)
        write(*,"(a12,4a10)") "time", "H_min_dep", "dt_dep", "n_max_dep", "dt_write"
        do i = 1, n 
            write(*,"(g12.3,f10.1,f10.1,i10,f10.1)") tpar%time(i), tpar%H_min_dep(i), tpar%dt_dep(i), &
                       tpar%n_max_dep(i), tpar%dt_write(i) 
        end do

        return 

    end subroutine tracer_par_trans_load 

    subroutine tracer_par_trans_allocate(tpar,n)

        implicit none 

        type(tracer_par_trans_class), intent(INOUT) :: tpar 
        integer, intent(IN) :: n 

        ! Make sure all arrays are deallocated first 
        if (allocated(tpar%time))      deallocate(tpar%time)
        if (allocated(tpar%H_min_dep)) deallocate(tpar%H_min_dep)
        if (allocated(tpar%dt_dep))    deallocate(tpar%dt_dep)
        if (allocated(tpar%n_max_dep)) deallocate(tpar%n_max_dep)
        if (allocated(tpar%dt_write))  deallocate(tpar%dt_write)
        
        allocate(tpar%time(n))
        allocate(tpar%H_min_dep(n))
        allocate(tpar%dt_dep(n))
        allocate(tpar%n_max_dep(n))
        allocate(tpar%dt_write(n))
        
        return 

    end subroutine tracer_par_trans_allocate

    subroutine tracer_allocate(now,dep,n)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        integer, intent(IN) :: n

        ! Make object is deallocated
        call tracer_deallocate(now,dep)

        ! Allocate tracer 
        allocate(now%active(n))
        allocate(now%id(n))
        allocate(now%parent(n),now%n_cloned(n))
        allocate(now%x(n),now%y(n),now%z(n),now%sigma(n))
        allocate(now%z_srf(n),now%dpth(n))
        allocate(now%ux(n),now%uy(n),now%uz(n))
        allocate(now%ax(n),now%ay(n),now%az(n))
        allocate(now%thk(n))
        allocate(now%T(n))
        allocate(now%H(n))

        ! Allocate deposition properties 
        
        allocate(dep%time(n), dep%H(n))
        allocate(dep%x(n), dep%y(n), dep%z(n), dep%lon(n), dep%lat(n))
        allocate(dep%t2m(n,nmon), dep%pr(n,nmon), dep%d18O(n,nmon))
        allocate(dep%t2m_ann(n), dep%pr_ann(n), dep%t2m_prann(n))
        allocate(dep%d18O_ann(n), dep%d18O_prann(n))
        
        return

    end subroutine tracer_allocate

    subroutine tracer_deallocate(now,dep)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        
        ! Deallocate state objects
        if (allocated(now%active))    deallocate(now%active)
        if (allocated(now%id))        deallocate(now%id)
        if (allocated(now%parent))    deallocate(now%parent)
        if (allocated(now%n_cloned))  deallocate(now%n_cloned)
        if (allocated(now%x))         deallocate(now%x)
        if (allocated(now%y))         deallocate(now%y)
        if (allocated(now%z))         deallocate(now%z)
        if (allocated(now%sigma))     deallocate(now%sigma)
        if (allocated(now%z_srf))     deallocate(now%z_srf)
        if (allocated(now%dpth))      deallocate(now%dpth)
        if (allocated(now%ux))        deallocate(now%ux)
        if (allocated(now%uy))        deallocate(now%uy)
        if (allocated(now%uz))        deallocate(now%uz)
        if (allocated(now%ax))        deallocate(now%ax)
        if (allocated(now%ay))        deallocate(now%ay)
        if (allocated(now%az))        deallocate(now%az)
        if (allocated(now%thk))       deallocate(now%thk)
        if (allocated(now%T))         deallocate(now%T)
        if (allocated(now%H))         deallocate(now%H)

        ! Deallocate deposition objects
        if (allocated(dep%time))      deallocate(dep%time)
        if (allocated(dep%z))         deallocate(dep%z)
        if (allocated(dep%H))         deallocate(dep%H)
        if (allocated(dep%x))         deallocate(dep%x)
        if (allocated(dep%y))         deallocate(dep%y)
        if (allocated(dep%lon))       deallocate(dep%lon)
        if (allocated(dep%lat))       deallocate(dep%lat)
        if (allocated(dep%t2m))        deallocate(dep%t2m)
        if (allocated(dep%pr))         deallocate(dep%pr)
        if (allocated(dep%d18O))       deallocate(dep%d18O)
        if (allocated(dep%t2m_ann))    deallocate(dep%t2m_ann)
        if (allocated(dep%pr_ann))     deallocate(dep%pr_ann)
        if (allocated(dep%t2m_prann))  deallocate(dep%t2m_prann)
        if (allocated(dep%d18O_ann))   deallocate(dep%d18O_ann)
        if (allocated(dep%d18O_prann)) deallocate(dep%d18O_prann)
        
        return

    end subroutine tracer_deallocate

    subroutine reshape2D_optional(idx_order,var1,nx,ny,var)
        ! Reshape an optional deposition tagging field. When the caller omitted
        ! it, produce an (nx,ny) field of MV instead, so downstream code has a
        ! correctly-shaped array to work with. An absent optional dummy may be
        ! passed straight through to another optional dummy, so this is called
        ! unconditionally.

        implicit none

        character(len=*), intent(IN) :: idx_order
        real(wp), intent(INOUT), allocatable :: var1(:,:)
        integer,    intent(IN) :: nx, ny
        real(wp), intent(IN), optional :: var(:,:)

        if (present(var)) then

            call tracer_reshape2D_field(idx_order,var,var1)

        else

            if (allocated(var1)) deallocate(var1)
            allocate(var1(nx,ny))
            var1 = MV

        end if

        return

    end subroutine reshape2D_optional

    subroutine reshape_monthly_optional(idx_order,var1,nx,ny,var)
        ! Reshape an optional monthly deposition field to (nx,ny,nmon). The month
        ! axis is the trailing dimension; each month's horizontal slice is reshaped
        ! with the same idx_order convention as the 2D tagging fields. When the
        ! caller omitted the field, produce an all-MV (nx,ny,nmon) array so the
        ! per-point deposition code has a correctly-shaped array to sample.

        implicit none

        character(len=*), intent(IN) :: idx_order
        real(wp), intent(INOUT), allocatable :: var1(:,:,:)
        integer,  intent(IN) :: nx, ny
        real(wp), intent(IN), optional :: var(:,:,:)

        integer :: m
        real(wp), allocatable :: slice1(:,:)

        if (allocated(var1)) deallocate(var1)
        allocate(var1(nx,ny,nmon))

        if (present(var)) then

            if (size(var,3) .ne. nmon) then
                write(0,*) "reshape_monthly_optional:: error: month axis /= nmon: ", size(var,3), nmon
                error stop
            end if

            do m = 1, nmon
                call tracer_reshape2D_field(idx_order,var(:,:,m),slice1)
                var1(:,:,m) = slice1
            end do

        else

            var1 = MV

        end if

        return

    end subroutine reshape_monthly_optional

    subroutine reshape_annual_to_monthly(idx_order,var1,nx,ny,var)
        ! Broadcast an annual (nx,ny) field across all nmon months, producing
        ! (nx,ny,nmon). Used for d18O when only the annual value is supplied.

        implicit none

        character(len=*), intent(IN) :: idx_order
        real(wp), intent(INOUT), allocatable :: var1(:,:,:)
        integer,  intent(IN) :: nx, ny
        real(wp), intent(IN) :: var(:,:)

        integer :: m
        real(wp), allocatable :: ann1(:,:)

        call tracer_reshape2D_field(idx_order,var,ann1)

        if (allocated(var1)) deallocate(var1)
        allocate(var1(nx,ny,nmon))

        do m = 1, nmon
            var1(:,:,m) = ann1
        end do

        return

    end subroutine reshape_annual_to_monthly

    subroutine tracer_reshape1D_vec(var,var1,rev)

        implicit none 
     
        real(wp),    intent(IN) :: var(:)
        real(wp), intent(INOUT), allocatable :: var1(:)
        logical,    intent(IN) :: rev 

        integer :: i, nx

        nx = size(var,1)
        if (allocated(var1)) deallocate(var1)
        allocate(var1(nx))

        if (rev) then 
            do i = 1, nx
                var1(i) = var(nx-i+1)
            end do 
        else 
            var1 = var 
        end if 

        return 

    end subroutine tracer_reshape1D_vec

    subroutine tracer_reshape2D_field(idx_order,var,var1)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(wp),    intent(IN) :: var(:,:)
        real(wp), intent(INOUT), allocatable :: var1(:,:)
        integer :: i, j
        integer :: nx, ny

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

                nx = size(var,1)
                ny = size(var,2)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny))

                var1 = var 

            case("kji")
                ! z, y, x array order 

                nx = size(var,2)
                ny = size(var,1)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny))

                do i = 1, nx 
                do j = 1, ny 
                    var1(i,j)  = var(j,i)
                end do 
                end do 

            case DEFAULT 

                write(0,*) "tracer_reshape2D_field:: error: unrecognized array order: ",trim(idx_order)
                write(0,*) "    Possible choices are: ijk, kji"
                error stop

        end select 

        return 

    end subroutine tracer_reshape2D_field

    subroutine tracer_reshape3D_field(idx_order,var,var1,rev_z)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(wp),    intent(IN) :: var(:,:,:)
        real(wp), intent(INOUT), allocatable :: var1(:,:,:)
        logical,    intent(IN) :: rev_z   ! Reverse the z-axis? 
        integer :: i, j, k
        integer :: nx, ny, nz 

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

                nx = size(var,1)
                ny = size(var,2)
                nz = size(var,3)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny,nz))

                if (rev_z) then 
                    do i = 1, nx 
                    do j = 1, ny  
                    do k = 1, nz 
                        var1(i,j,k)  = var(i,j,nz-k+1) 
                    end do 
                    end do 
                    end do 
                else 
                    var1 = var 
                end if 

            case("kji")
                ! z, y, x array order 

                nx = size(var,3)
                ny = size(var,2)
                nz = size(var,1)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny,nz))

                if (rev_z) then 
                    do i = 1, nx 
                    do j = 1, ny 
                    do k = 1, nz 
                        var1(i,j,k)  = var(nz-k+1,j,i)
                    end do 
                    end do 
                    end do 
                else 
                    do i = 1, nx 
                    do j = 1, ny 
                    do k = 1, nz 
                        var1(i,j,k)  = var(k,j,i) 
                    end do 
                    end do 
                    end do 
                end if 

            case DEFAULT 

                write(0,*) "tracer_reshape3D_field:: error: unrecognized array order: ",trim(idx_order)
                write(0,*) "    Possible choices are: ijk, kji"
                error stop

        end select 

        return 

    end subroutine tracer_reshape3D_field

    subroutine which(x,ind,stat)
        ! Analagous to R::which function
        ! Returns indices that match condition x==.TRUE.

        implicit none 

        logical :: x(:)
        integer, allocatable :: tmp(:), ind(:)
        integer, optional :: stat  
        integer :: n, i  

        n = count(x)
        allocate(tmp(n))
        tmp = 0 

        n = 0
        do i = 1, size(x) 
            if (x(i)) then 
                n = n+1
                tmp(n) = i 
            end if
        end do 

        if (present(stat)) stat = n 

        if (allocated(ind)) deallocate(ind)

        if (n .eq. 0) then 
            allocate(ind(1))
            ind = -1 
        else
            allocate(ind(n))
            ind = tmp(1:n)
        end if 
        
        return 

    end subroutine which

end module tracer3D 


