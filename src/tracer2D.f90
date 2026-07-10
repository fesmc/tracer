
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

        real(prec) :: y(NY_GHOST)
        integer    :: j

        do j = 1, NY_GHOST
            y(j) = real(j-1,prec) / real(NY_GHOST-1,prec)
        end do

        return

    end function ghost_yaxis


    subroutine tracer2D_init(trc,filename,time,x,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:)
        logical,    intent(IN) :: is_sigma 
        real(prec_time) :: time 

        real(prec) :: y(NY_GHOST)

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
                               lon,lat,t2m_ann,t2m_sum,pr_ann,pr_sum,d18O_ann, &
                               dep_now,stats_now)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), z(:)
        real(prec), intent(IN) :: z_srf(:), H(:)
        real(prec), intent(IN) :: ux(:,:), uz(:,:)

        ! Deposition tagging fields, all optional (see tracer_update).
        real(prec), intent(IN), optional :: lon(:), lat(:)
        real(prec), intent(IN), optional :: t2m_ann(:), t2m_sum(:)
        real(prec), intent(IN), optional :: pr_ann(:), pr_sum(:)
        real(prec), intent(IN), optional :: d18O_ann(:)

        logical,    intent(IN) :: dep_now, stats_now

        ! Local variables
        real(prec) :: y(NY_GHOST)
        real(prec), allocatable :: z_srf_2D(:,:), H_2D(:,:)
        real(prec), allocatable :: ux_3D(:,:,:), uy_3D(:,:,:), uz_3D(:,:,:)
        real(prec), allocatable :: lon_2D(:,:), lat_2D(:,:), t2m_ann_2D(:,:), t2m_sum_2D(:,:) 
        real(prec), allocatable :: pr_ann_2D(:,:), pr_sum_2D(:,:), d18O_ann_2D(:,:)
        
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
        if (present(t2m_ann))  allocate(t2m_ann_2D(size(x,1),ny))
        if (present(t2m_sum))  allocate(t2m_sum_2D(size(x,1),ny))
        if (present(pr_ann))   allocate(pr_ann_2D(size(x,1),ny))
        if (present(pr_sum))   allocate(pr_sum_2D(size(x,1),ny))
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
            if (present(t2m_ann))  t2m_ann_2D(:,j)  = t2m_ann
            if (present(t2m_sum))  t2m_sum_2D(:,j)  = t2m_sum
            if (present(pr_ann))   pr_ann_2D(:,j)   = pr_ann
            if (present(pr_sum))   pr_sum_2D(:,j)   = pr_sum
            if (present(d18O_ann)) d18O_ann_2D(:,j) = d18O_ann

        end do

        ! Now update tracers using 3D call 
        call tracer_update(trc,time,x,y,z,z_srf_2D,H_2D,ux_3D,uy_3D,uz_3D, &
                            lon_2D,lat_2D,t2m_ann_2D,t2m_sum_2D,pr_ann_2D,pr_sum_2D,d18O_ann_2D, &
                            dep_now,stats_now,order="ijk")

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

