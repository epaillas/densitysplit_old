program density_profiles
    implicit none
    
    real*8 :: rgrid_x, rgrid_y, rgrid_z, vol, rhomed
    real*8 :: boxsize_x, boxsize_y, boxsize_z
    real*8 :: disx, disy, disz, dis
    real*8 :: xvc, yvc, zvc
    real*8 :: dmax, dmin, rfilter
    real*8 :: pi = 4.*atan(1.)
    real*8 :: qperp, qpara
    
    integer*8 :: ng, nc
    integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: indx, indy, indz, nrows, ncols
    integer*8 :: ipx, ipy, ipz, ndif_x, ndif_y, ndif_z
    integer*8 :: ngrid_x, ngrid_y, ngrid_z
    integer*4 :: argstat1, argstat2
    
    integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
    integer*8, dimension(:), allocatable :: ll
    
    real*8, dimension(3) :: r
    real*8, allocatable, dimension(:,:)  :: tracers, centres
    real*8, dimension(:), allocatable :: DD, delta
  
    logical :: has_velocity = .false.
    
    character(20), external :: str
    character(len=500) :: input_tracers, input_centres, output_den
    character(len=10) :: dmax_char, dmin_char
    character(len=10) :: boxchar_x, boxchar_y, boxchar_z, rfilter_char
    character(len=10) :: ngridchar_x, ngridchar_y, ngridchar_z
    character(len=10) :: qperp_char, qpara_char
    
    if (iargc() .lt. 12) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) input_data'
        write(*,*) '2) input_centres'
        write(*,*) '3) output_den'
        write(*,*) '4) boxsize_x'
        write(*,*) '5) boxsize_y'
        write(*,*) '6) boxsize_z'
        write(*,*) '7) dmin'
        write(*,*) '8) dmax'
        write(*,*) '9) rfilter'
        write(*,*) '10) ngrid_x'
        write(*,*) '11) ngrid_y'
        write(*,*) '12) ngrid_z'
        write(*,*) '13) qperp (optional)'
        write(*,*) '14) qpara (optional)'
        write(*,*) ''
        stop
      end if
      
    call get_command_argument(number=1, value=input_tracers)
    call get_command_argument(number=2, value=input_centres)
    call get_command_argument(number=3, value=output_den)
    call get_command_argument(number=4, value=boxchar_x)
    call get_command_argument(number=5, value=boxchar_y)
    call get_command_argument(number=6, value=boxchar_z)
    call get_command_argument(number=7, value=dmin_char)
    call get_command_argument(number=8, value=dmax_char)
    call get_command_argument(number=9, value=rfilter_char)
    call get_command_argument(number=10, value=ngridchar_x)
    call get_command_argument(number=11, value=ngridchar_y)
    call get_command_argument(number=12, value=ngridchar_z)
    call get_command_argument(number=13, value=qperp_char, status=argstat1)
    call get_command_argument(number=14, value=qpara_char, status=argstat2)
    
    read(boxchar_x, *) boxsize_x
    read(boxchar_y, *) boxsize_y
    read(boxchar_z, *) boxsize_z
    read(dmin_char, *) dmin
    read(dmax_char, *) dmax
    read(rfilter_char, *) rfilter
    read(ngridchar_x, *) ngrid_x
    read(ngridchar_y, *) ngrid_y
    read(ngridchar_z, *) ngrid_z

    if (argstat1 == 0 .and. argstat2 == 0) then
      read(qperp_char, *) qperp
      read(qpara_char, *) qpara
    else
      qperp = 1.0
      qpara = 1.0
    end if
    
    write(*,*) '-----------------------'
    write(*,*) 'Running tophat_filter.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'input_tracers: ', trim(input_tracers)
    write(*, *) 'input_centres: ', trim(input_centres)
    write(*, *) 'boxsize_x: ', trim(boxchar_x)
    write(*, *) 'boxsize_y: ', trim(boxchar_y)
    write(*, *) 'boxsize_z: ', trim(boxchar_z)
    write(*, *) 'output_den: ', trim(output_den)
    write(*, *) 'dmin: ', trim(dmin_char), ' Mpc'
    write(*, *) 'dmax: ', trim(dmax_char), ' Mpc'
    write(*, *) 'rfilter: ', trim(rfilter_char), 'Mpc'
    write(*, *) 'ngrid_x: ', trim(ngridchar_x)
    write(*, *) 'ngrid_y: ', trim(ngridchar_y)
    write(*, *) 'ngrid_z: ', trim(ngridchar_z)
    write(*,*) ''
  
    open(10, file=input_tracers, status='old', form='unformatted')
    read(10) nrows
    read(10) ncols
    allocate(tracers(ncols, nrows))
    read(10) tracers
    close(10)
    ng = nrows
    if (ncols .eq. 6) then
      has_velocity = .true.
      write(*,*) 'Tracer file has velocity information.'
    end if
    write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
    write(*,*) 'tracers(min), tracers(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))
  
    open(11, file=input_centres, status='old', form='unformatted')
    read(11) nrows
    read(11) ncols
    allocate(centres(ncols, nrows))
    read(11) centres
    close(11)
    nc = nrows
    write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)
    write(*,*) 'centres(min), centres(max) = ', minval(centres(1,:)), maxval(centres(1,:))

    ! Account for potential geometrical distortions
    boxsize_x = boxsize_x / qperp
    boxsize_y = boxsize_y / qperp
    boxsize_z = boxsize_z / qpara

    tracers(1,:) = tracers(1,:) / qperp
    tracers(2,:) = tracers(2,:) / qperp
    tracers(3,:) = tracers(3,:) / qpara

    centres(1,:) = centres(1,:) / qperp
    centres(2,:) = centres(2,:) / qperp
    centres(3,:) = centres(3,:) / qpara


    if (qperp .ne. 1.0 .or. qpara .ne. 1.0) then
      write(*,*) 'Positions have been shifted due to geometrical distortions'
      write(*,*) 'qperp, qpara: ', qperp, qpara
      write(*,*) 'boxsize_x: ', boxsize_x
      write(*,*) 'boxsize_y: ', boxsize_y
      write(*,*) 'boxsize_z: ', boxsize_z
      write(*,*) 'tracers(min), tracers(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))
      write(*,*) 'centres(min), centres(max) = ', minval(centres(1,:)), maxval(centres(1,:))
    end if
  
    allocate(DD(nc))
    allocate(delta(nc))
    
    ! Mean density inside the box
    rhomed = ng / (boxsize_x * boxsize_y * boxsize_z)
    
    ! Construct linked list for tracers
    write(*,*) ''
    write(*,*) 'Constructing linked list...'
    allocate(lirst(ngrid_x, ngrid_y, ngrid_z))
    allocate(nlirst(ngrid_x, ngrid_y, ngrid_z))
    allocate(ll(ng))
    rgrid_x = (boxsize_x) / real(ngrid_x)
    rgrid_y = (boxsize_y) / real(ngrid_y)
    rgrid_z = (boxsize_z) / real(ngrid_z)
    
    lirst = 0
    ll = 0
    
    do i = 1, ng
      indx = int((tracers(1, i)) / rgrid_x + 1.)
      indy = int((tracers(2, i)) / rgrid_y + 1.)
      indz = int((tracers(3, i)) / rgrid_z + 1.)
    
      if(indx.gt.0.and.indx.le.ngrid_x.and.indy.gt.0.and.indy.le.ngrid_y.and.&
      indz.gt.0.and.indz.le.ngrid_z)lirst(indx,indy,indz)=i
    
      if(indx.gt.0.and.indx.le.ngrid_x.and.indy.gt.0.and.indy.le.ngrid_y.and.&
      indz.gt.0.and.indz.le.ngrid_z)nlirst(indx,indy,indz) = &
      nlirst(indx, indy, indz) + 1
    end do
    
    do i = 1, ng
      indx = int((tracers(1, i))/ rgrid_x + 1.)
      indy = int((tracers(2, i))/ rgrid_y + 1.)
      indz = int((tracers(3, i))/ rgrid_z + 1.)
      if(indx.gt.0.and.indx.le.ngrid_x.and.indy.gt.0.and.indy.le.ngrid_y.and.&
      &indz.gt.0.and.indz.le.ngrid_z) then
        ll(lirst(indx,indy,indz)) = i
        lirst(indx,indy,indz) = i
      endif
    end do
    
    write(*,*) 'Linked list successfully constructed'
    write(*,*) ''
    write(*,*) 'Starting loop over centres...'
    
    DD = 0
    delta = 0
    
    do i = 1, nc
      xvc = centres(1, i)
      yvc = centres(2, i)
      zvc = centres(3, i)
  
      ipx = int((xvc) / rgrid_x + 1.)
      ipy = int((yvc) / rgrid_y + 1.)
      ipz = int((zvc) / rgrid_z + 1.)
  
      ndif_x = int(dmax / rgrid_x + 1.)
      ndif_y = int(dmax / rgrid_y + 1.)
      ndif_z = int(dmax / rgrid_z + 1.)
    
      do ix = ipx - ndif_x, ipx + ndif_x
        do iy = ipy - ndif_y, ipy + ndif_y
          do iz = ipz - ndif_z, ipz + ndif_z
    
            ix2 = ix
            iy2 = iy
            iz2 = iz
    
            if (ix2 .gt. ngrid_x) ix2 = ix2 - ngrid_x
            if (ix2 .lt. 1) ix2 = ix2 + ngrid_x
            if (iy2 .gt. ngrid_y) iy2 = iy2 - ngrid_y
            if (iy2 .lt. 1) iy2 = iy2 + ngrid_y
            if (iz2 .gt. ngrid_z) iz2 = iz2 - ngrid_z
            if (iz2 .lt. 1) iz2 = iz2 + ngrid_z
    
            ii = lirst(ix2,iy2,iz2)
            if(ii.ne.0) then
              do
                ii = ll(ii)
                disx = tracers(1, ii) - xvc
                disy = tracers(2, ii) - yvc
                disz = tracers(3, ii) - zvc
  
                if (disx .lt. -boxsize_x/2) disx = disx + boxsize_x
                if (disx .gt. boxsize_x/2) disx = disx - boxsize_x
                if (disy .lt. -boxsize_y/2) disy = disy + boxsize_y
                if (disy .gt. boxsize_y/2) disy = disy - boxsize_y
                if (disz .lt. -boxsize_z/2) disz = disz + boxsize_z
                if (disz .gt. boxsize_z/2) disz = disz - boxsize_z
    
                r = (/ disx, disy, disz /)
                dis = norm2(r)
  
                if (dis .gt. dmin .and. dis .lt. dmax) then
                  DD(i) = DD(i) + 1
                end if
  
    
                if(ii.eq.lirst(ix2,iy2,iz2)) exit
    
              end do
            end if
          end do
        end do
      end do
  
  
    vol = 4./3 * pi * (dmax ** 3 - dmin ** 3)
    delta(i) = DD(i) / (vol * rhomed) - 1
  
    end do
    
    write(*,*) ''
    write(*,*) 'Calculation finished. Writing output...'
    
    open(12, file=output_den, status='replace', form='unformatted')
  
    write(12) nc
    write(12) delta
  
    end program density_profiles
    