
program tracertest 

    use nml 
    use ncio 
    use tracer 

    implicit none 

    type(tracer_class) :: trc1
    character(len=128) :: fldr, prefix, filename_nml 

    type profile_class 
        integer :: nx, nz
        integer, allocatable :: dims(:)
        real(wp) :: L, A, H0, G, B, M  
        real(wp), allocatable :: xc(:), sigma(:)
        real(wp), allocatable :: zs(:), zb(:), H(:), dHdx(:), ux(:,:), uz(:,:)
        real(wp), allocatable :: age(:)

        character(len=512) :: filename
    end type 

    type(profile_class) :: prof1 
    
    integer :: k, kmax, q 
    real(wp) :: time, time_start, time_end  
    logical    :: dep_now, write_now

    ! The namelist to run may be given as the first command-line argument, so
    ! that the RH2003 and RH2003_analytic cases can be run without renaming
    ! files over each other. Output lands in output/<namelist-stem>.
    call get_command_argument(1,filename_nml)
    if (len_trim(filename_nml) .eq. 0) filename_nml = "RH2003.nml"

    prefix = nml_stem(filename_nml)
    fldr   = "output/"//trim(prefix)

    ! Simulation parameters
    time_start = -160000.0
    time_end   = 0.0

    call calc_profile_RH2003(prof1,filename_nml)
    call profile_write(prof1,fldr=fldr,filename="profile_RH2003.nc")

    ! Initialize tracer and output file 
    call tracer2D_init(trc1,filename_nml,time=real(time_start,prec_time),x=prof1%xc,is_sigma=.TRUE.)

    prof1%filename = gen_filename(prof1,trc1%par%dt,trc1%par%interp_method)
    call tracer2D_write_init(trc1,fldr,prof1%filename)

    do k = 1, int((time_end-time_start)/trc1%par%dt)+1

        time = time_start + trc1%par%dt*(k-1)

        dep_now  = .FALSE.
        if (mod(time,trc1%par%dt_dep) .eq. 0.0) dep_now = .TRUE.

        ! This idealized profile carries no climate forcing, so the optional
        ! deposition tagging fields are omitted; they are recorded as MV.
        call tracer2D_update(trc1,time=real(time,prec_time), &
                             x=prof1%xc,z=prof1%sigma,z_srf=prof1%H,H=prof1%H, &
                             ux=prof1%ux,uz=prof1%uz, &
                             dep_now=dep_now,stats_now=.FALSE.)

        write_now = .FALSE. 
        if (mod(time,trc1%par%dt_write) .eq. 0.0) write_now = .TRUE. 

        if (write_now) then 
            call tracer2D_write(trc1,real(time,prec_time),fldr,prof1%filename)
        end if 

        if (mod(time,5000.0) .eq. 0.0) then 
            write(*,*) "time = ", time, trc1%par%dt, trc1%par%n_active
        end if 
    end do 

!     ! Write stats 
!     call tracer2D_write_stats(trc1,time,fldr,filename_stats)

contains

    function nml_stem(filename) result(stem)
        ! Strip any directory prefix and the trailing ".nml" from a namelist
        ! path, e.g. "cfg/RH2003_analytic.nml" -> "RH2003_analytic".

        implicit none

        character(len=*), intent(IN) :: filename
        character(len=128) :: stem

        integer :: i, n

        stem = adjustl(filename)

        i = index(stem,"/",back=.TRUE.)
        if (i .gt. 0) stem = stem(i+1:)

        n = len_trim(stem)
        if (n .gt. 4) then
            if (stem(n-3:n) .eq. ".nml") stem = stem(1:n-4)
        end if

        return

    end function nml_stem

    subroutine calc_profile_RH2003(prof,filename)
        ! Define a 2D profile (x-z) following 
        ! Rybak and Huybrechts (2003, Annals of Glaciology)

        implicit none 

        type(profile_class), intent(INOUT) :: prof 
        character(len=*), intent(IN) :: filename

        ! Local parameters 
!         integer, parameter :: nx = 51 
!         integer, parameter :: nz = 101 
        integer,    parameter :: ng  = 3          ! exponent
        real(wp), parameter :: rho = 910.0     ! kg/m^3
        real(wp), parameter :: gg  = 9.81        ! m/s
        real(wp), parameter :: A   = 10.0**(-16) ! Pa^3/a
!         real(wp), parameter :: M = 0.1         ! m/a
!         real(wp), parameter :: L = 10.0**6     ! m 

        ! Loaded parameters
        integer :: nx, nz 
        real(wp) :: L, G, B 
        real(wp) :: ux_fac 

        ! Local variables 
        integer :: i, j 
        real(wp) :: H0, dHdx, H 
        real(wp) :: M
        
        ! Load parameters 
        call nml_read(filename,"rh_par","nx",nx)
        call nml_read(filename,"rh_par","nz",nz)
        call nml_read(filename,"rh_par","L", L)
        call nml_read(filename,"rh_par","G", G)
        call nml_read(filename,"rh_par","B", B)
        call nml_read(filename,"rh_par","ux_fac", ux_fac)
        
        M = G - B  

        prof%nx = nx    ! Right-hand side symmetrical of domain
        prof%nz = nz
        
        allocate(prof%xc(prof%nx))
        allocate(prof%sigma(prof%nz))
        allocate(prof%zs(prof%nx))
        allocate(prof%zb(prof%nx))
        allocate(prof%H(prof%nx))
        allocate(prof%dHdx(prof%nx)) 

        allocate(prof%ux(prof%nx,prof%nz))
        allocate(prof%uz(prof%nx,prof%nz))

        allocate(prof%age(prof%nz)) 

        ! Define x-dimension (0-1000 km) 
        do i = 1, prof%nx 
            prof%xc(i) = 0.0 + (i-1)/real(nx-1) * 1000.0 
        end do

        prof%xc = prof%xc*1e3   ! [km] => [m] 

        ! Define z-dimension (0-1 sigma)
        do i = 1, prof%nz 
            prof%sigma(i) = 0.0 + (i-1)/real(prof%nz-1) * 1.0 
        end do 
        prof%sigma(1) = 1e-8   ! To avoid singularities 

        ! Calculate H0 (should be H0=3598.4 m)
        H0 = (20.0*M/A)**(1.0/(2.0*(ng+1)))*(1/(rho*gg))**(ng/(2.0*(ng+1)))*L**(1.0/2.0)
        
        ! Calculate H(x) 
        prof%H = H0 * (1.0-(abs(prof%xc)/L)**((ng+1.0)/ng))**(real(ng)/(2.0*(ng+1.0)))

        prof%ux   = 0.0 
        prof%uz   = 0.0 
        prof%age  = 0.0 
        prof%dHdx = 0.0

        ! Calculate slope 
        do i = 2, prof%nx-1 
            prof%dHdx(i) = (prof%H(i+1)-prof%H(i-1))/(prof%xc(i+1)-prof%xc(i-1))
        end do 

        ! Calculate velocities 
        do i = 1, prof%nx 

            H    = prof%H(i)
            dHdx = prof%dHdx(i)

            do j = 1, prof%nz

                prof%ux(i,j) = -(2.0*A)/(ng+1.0)*(rho*gg)**ng * dHdx**(ng-1.0) &
                  * dHdx * (H**(ng+1.0)-(H-prof%sigma(j)*H)**(ng+1.0))

                prof%ux(i,j) = prof%ux(i,j)*ux_fac   ! Sets ux=0 for analytic solution at summit

                prof%uz(i,j) = prof%sigma(j)*(-G+B+prof%ux(i,j)*dHdx) - B

            end do 
        end do 


        ! Calculate analytical age at the divide
        prof%age = -(H0/G)*log(prof%sigma)
        
!         ! Write summary 
!         write(*,"(a,500f8.2)") "xc: ",    prof%xc*1e-3 
!         write(*,"(a,500f8.2)") "sigma: ", prof%sigma 
!         write(*,"(a,f10.2)") "H0 = ", H0 
!         write(*,"(a,2f10.2)") "range(ux): ", minval(prof%ux), maxval(prof%ux)
!         write(*,"(a,2f10.2)") "range(uz): ", minval(prof%uz), maxval(prof%uz) 
        
        ! Additionally store constants and parameters 
        prof%L  = L 
        prof%nx = nx 
        prof%nz = nz 
        prof%H0 = H0 
        prof%A  = A 
        prof%G  = G 
        prof%B  = B 
        prof%M  = M 
        
        
        return 

    end subroutine calc_profile_RH2003

    function gen_filename(prof,dt,interp_method) result(filename)

        implicit none

        type(profile_class), intent(IN) :: prof
        real(prec_wrt),      intent(IN) :: dt
        character(len=*),    intent(IN) :: interp_method
        character(len=512) :: filename

        character(len=5) :: str_nx, str_nz, str_prec, str_dt

        write(str_nx,"(i5)") prof%nx
        str_nx = adjustl(str_nx)
        write(str_nz,"(i5)") prof%nz
        str_nz = adjustl(str_nz)
        write(str_prec,*) "sp"
        if (kind(1.d0)==wp) write(str_prec,*) "dp"
        str_prec = adjustl(str_prec)

        write(str_dt,"(f5.1)") dt
        str_dt = adjustl(str_dt)

        ! Generate filename based on parameter values. interp_method is part of
        ! the stem so that a linear and a spline run of the same configuration
        ! do not overwrite each other.
        write(filename,"(a,a,a1,a,a1,f4.2,a1,f4.2,a1,a,a1,a,a1,a,a3)") &
            "RH2003_", trim(str_nx), "_", trim(str_nz),"_", prof%G, "_", prof%B,"_", &
            trim(str_prec),"_",trim(interp_method),"_",trim(str_dt),".nc"

        write(*,*) "Filename: ", trim(filename)

        return

    end function gen_filename

    subroutine profile_write(prof,fldr,filename)

        implicit none 

        type(profile_class), intent(IN) :: prof 
        character(len=*), intent(IN)    :: fldr, filename 

        ! Local variables 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",x=prof%xc*1e-3,units="kilometers")
        call nc_write_dim(path_out,"sigma",x=prof%sigma,units="1")
        call nc_write_dim(path_out,"par",x=1)

        ! Write parameters 
        call nc_write(path_out,"nx",  prof%nx,  dim1="par")
        call nc_write(path_out,"nz",  prof%nz,  dim1="par")
        call nc_write(path_out,"L",   prof%L,   dim1="par")
        call nc_write(path_out,"H0",  prof%H0,  dim1="par")
        call nc_write(path_out,"A",   prof%A,   dim1="par")
        call nc_write(path_out,"G",   prof%G,   dim1="par")
        call nc_write(path_out,"B",   prof%B,   dim1="par")
        call nc_write(path_out,"M",   prof%M,   dim1="par")
        
        call nc_write(path_out,"H",   prof%H,   dim1="xc",missing_value=MV)
        call nc_write(path_out,"dHdx",prof%dHdx,dim1="xc",missing_value=MV)
        call nc_write(path_out,"ux",prof%ux,dim1="xc",dim2="sigma",missing_value=MV)
        call nc_write(path_out,"uz",prof%uz,dim1="xc",dim2="sigma",missing_value=MV)
        call nc_write(path_out,"umag",sqrt(prof%ux**2+prof%uz**2),dim1="xc",dim2="sigma",missing_value=MV)
        call nc_write(path_out,"age",prof%age,dim1="sigma",missing_value=MV)
        

        return 

    end subroutine profile_write 

end program tracertest 

