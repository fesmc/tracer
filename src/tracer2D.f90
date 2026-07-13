
module tracer2D 
    ! Module to wrap a 2D (profile) version of the tracer model
    ! Makes calls to main tracer code by reshaping profile into
    ! 3D array with y-dimension thickness of 1. 

    use tracer_precision
    use tracer3D
    use ncio    

    implicit none

    ! Length of the ghost y-dimension. Three is the smallest axis bspline-fortran
    ! will build a spline on, so interp_method="spline" needs at least this many
    ! even though every field is constant along y.
    integer, parameter :: NY_GHOST = 3

    private
    public :: tracer2D_init
    public :: tracer2D_update
    public :: tracer2D_end

contains

    function ghost_yaxis() result(y)
        ! The ghost y-axis, spanning [0,1] with NY_GHOST evenly spaced points.

        implicit none

        real(wp) :: y(NY_GHOST)
        integer    :: j

        do j = 1, NY_GHOST
            y(j) = real(j-1,wp) / real(NY_GHOST-1,wp)
        end do

        return

    end function ghost_yaxis


    subroutine tracer2D_init(trc,filename,time,x,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(wp), intent(IN) :: x(:)
        logical,    intent(IN) :: is_sigma 
        real(prec_time) :: time 

        real(wp) :: y(NY_GHOST)

        ! Define the ghost y-dimension. Must match tracer2D_update.
        y = ghost_yaxis()

        ! Call 3D tracer_init
        call tracer_init(trc,filename,time,x,y,is_sigma)

        ! Mark the domain as a profile, so that tracer_activate does not jitter
        ! deposition locations along the ghost y-axis.
        trc%par%is_profile = .TRUE.

        return 

    end subroutine tracer2D_init


    subroutine tracer2D_update(trc,time,x,z,z_srf,H,ux,uz, &
                               lon,lat,t2m,pr,d18O,d18O_ann, &
                               dep_now,stats_now)

        implicit none

        type(tracer_class), intent(INOUT) :: trc
        real(prec_time), intent(IN) :: time
        real(wp), intent(IN) :: x(:), z(:)
        real(wp), intent(IN) :: z_srf(:), H(:)
        real(wp), intent(IN) :: ux(:,:), uz(:,:)

        ! Deposition tagging fields, all optional (see tracer_update). The
        ! monthly climate profiles are (nx,nmon); d18O may instead be given as
        ! the annual profile d18O_ann (nx).
        real(wp), intent(IN), optional :: lon(:), lat(:)
        real(wp), intent(IN), optional :: t2m(:,:), pr(:,:), d18O(:,:)
        real(wp), intent(IN), optional :: d18O_ann(:)

        logical,    intent(IN) :: dep_now, stats_now

        ! Local variables
        real(wp) :: y(NY_GHOST)
        real(wp), allocatable :: z_srf_2D(:,:), H_2D(:,:)
        real(wp), allocatable :: ux_3D(:,:,:), uy_3D(:,:,:), uz_3D(:,:,:)
        real(wp), allocatable :: lon_2D(:,:), lat_2D(:,:), d18O_ann_2D(:,:)
        real(wp), allocatable :: t2m_3D(:,:,:), pr_3D(:,:,:), d18O_3D(:,:,:)   ! (nx,ny,nmon)

        integer :: j, ny

        ny = size(y,1)

        ! Define ghost dimension and data
        allocate(z_srf_2D(size(x,1),ny))
        allocate(H_2D(size(x,1),ny))
        allocate(ux_3D(size(ux,1),ny,size(ux,2)))
        allocate(uy_3D(size(ux,1),ny,size(ux,2)))
        allocate(uz_3D(size(ux,1),ny,size(ux,2)))
        ! An omitted tagging field is left unallocated. Passing an unallocated
        ! allocatable to tracer_update's optional dummy makes it absent there,
        ! so absence propagates through this wrapper.
        if (present(lon))      allocate(lon_2D(size(x,1),ny))
        if (present(lat))      allocate(lat_2D(size(x,1),ny))
        if (present(t2m))      allocate(t2m_3D(size(x,1),ny,size(t2m,2)))
        if (present(pr))       allocate(pr_3D(size(x,1),ny,size(pr,2)))
        if (present(d18O))     allocate(d18O_3D(size(x,1),ny,size(d18O,2)))
        if (present(d18O_ann)) allocate(d18O_ann_2D(size(x,1),ny))

        ! Define the ghost y-dimension. Must match tracer2D_init.
        y = ghost_yaxis()

        ! Reshape input data, holding every field constant along the ghost y-axis
        do j = 1, size(y)

            z_srf_2D(:,j)    = z_srf
            H_2D(:,j)        = H

            ux_3D(:,j,:)     = ux
            uy_3D            = 0.0
            uz_3D(:,j,:)     = uz

            if (present(lon))      lon_2D(:,j)      = lon
            if (present(lat))      lat_2D(:,j)      = lat
            if (present(t2m))      t2m_3D(:,j,:)    = t2m
            if (present(pr))       pr_3D(:,j,:)     = pr
            if (present(d18O))     d18O_3D(:,j,:)   = d18O
            if (present(d18O_ann)) d18O_ann_2D(:,j) = d18O_ann

        end do

        ! Now update tracers using 3D call. Velocities are already on aa-nodes
        ! here, so the native staggered axes (x_ux/y_uy/z_uz) are left absent.
        call tracer_update(trc,time,x,y,z,z_srf_2D,H_2D,ux_3D,uy_3D,uz_3D, &
                            lon=lon_2D,lat=lat_2D,t2m=t2m_3D,pr=pr_3D,d18O=d18O_3D, &
                            d18O_ann=d18O_ann_2D, &
                            dep_now=dep_now,stats_now=stats_now,order="ijk")

        return

    end subroutine tracer2D_update

    subroutine tracer2D_end(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        
        ! Call normal tracer_end subroutine 
        call tracer_end(trc) 

        return 

    end subroutine tracer2D_end


end module tracer2D

