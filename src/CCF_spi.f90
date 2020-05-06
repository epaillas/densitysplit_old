program density_profiles
    implicit none
    
    real*4 :: rgrid, boxsize, diff_vol, rhomed
    real*4 :: disx, disy, disz, para, perp
    real*4 :: xvc, yvc, zvc
    real*4 :: perpwidth, perpmax, perpmin
    real*4 :: parawidth, paramin, paramax
    real*4 :: pi = 4.*atan(1.)
    
    integer*4 :: ng, nc, nperpbin, perpind, nparabin, paraind
    integer*4 :: i, ii, jj, ix, iy, iz, ix2, iy2, iz2
    integer*4 :: indx, indy, indz, nrows, ncols
    integer*4 :: ipx, ipy, ipz, ndif
    integer*4 :: ngrid
    
    integer*4, dimension(:, :, :), allocatable :: lirst, nlirst
    integer*4, dimension(:), allocatable :: ll
    
    real*4, dimension(3) :: r, com
    real*4, allocatable, dimension(:,:)  :: tracers, centres
    real*4, dimension(:, :, :), allocatable :: DD, cum_DD, delta, cum_delta
    real*4, dimension(:), allocatable :: perpbin, perpbin_edges, parabin, parabin_edges
  
    logical :: has_velocity = .false.
    
    character(20), external :: str
    character(len=500) :: input_tracers, input_centres, output_den
    character(len=10) :: perpmax_char, perpmin_char, nperpbin_char, ngrid_char, box_char
    
    if (iargc() .ne. 8) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) input_data'
        write(*,*) '2) input_centres'
        write(*,*) '3) output_den'
        write(*,*) '4) boxsize'
        write(*,*) '5) perpmin'
        write(*,*) '6) perpmax'
        write(*,*) '7) nperpbin'
        write(*,*) '8) ngrid'
        write(*,*) ''
        stop
      end if
      
    call getarg(1, input_tracers)
    call getarg(2, input_centres)
    call getarg(3, output_den)
    call getarg(4, box_char)
    call getarg(5, perpmin_char)
    call getarg(6, perpmax_char)
    call getarg(7, nperpbin_char)
    call getarg(8, ngrid_char)
    
    read(box_char, *) boxsize
    read(perpmin_char, *) perpmin
    read(perpmax_char, *) perpmax
    read(nperpbin_char, *) nperpbin
    read(ngrid_char, *) ngrid
    
    write(*,*) '-----------------------'
    write(*,*) 'Running CCF_spi.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'input_tracers: ', trim(input_tracers)
    write(*, *) 'input_centres: ', trim(input_centres)
    write(*, *) 'boxsize: ', trim(box_char)
    write(*, *) 'output_den: ', trim(output_den)
    write(*, *) 'perpmin: ', trim(perpmin_char), ' Mpc'
    write(*, *) 'perpmax: ', trim(perpmax_char), ' Mpc'
    write(*, *) 'nperpbin: ', trim(nperpbin_char)
    write(*, *) 'ngrid: ', trim(ngrid_char)
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
    write(*,*) 'pos(min), pos(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))
  
    open(11, file=input_centres, status='old', form='unformatted')
    read(11) nrows
    read(11) ncols
    allocate(centres(ncols, nrows))
    read(11) centres
    close(11)
    nc = nrows
    write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)
  
    nparabin = nperpbin
    allocate(perpbin(nperpbin))
    allocate(parabin(nparabin))
    allocate(perpbin_edges(nperpbin + 1))
    allocate(parabin_edges(nparabin + 1))
    allocate(DD(nc, nperpbin, nparabin))
    allocate(cum_DD(nc, nperpbin, nparabin))
    allocate(delta(nc, nperpbin, nparabin))
    allocate(cum_delta(nc, nperpbin, nparabin))
    
    
    perpwidth = (perpmax - perpmin) / nperpbin
    do i = 1, nperpbin + 1
      perpbin_edges(i) = perpmin+(i-1)*perpwidth
    end do
    do i = 1, nperpbin
      perpbin(i) = perpbin_edges(i+1)-perpwidth/2.
    end do
  
    paramin = perpmin
    paramax = perpmax
  
    parawidth = (paramax - paramin) / nparabin
    do i = 1, nparabin + 1
      parabin_edges(i) = paramin+(i-1)*parawidth
    end do
    do i = 1, nparabin
      parabin(i) = parabin_edges(i+1)-parawidth/2.
    end do
    
    ! Mean density inside the box
    rhomed = ng / (boxsize ** 3)
    
    ! Construct linked list for tracers
    write(*,*) ''
    write(*,*) 'Constructing linked list...'
    allocate(lirst(ngrid, ngrid, ngrid))
    allocate(nlirst(ngrid, ngrid, ngrid))
    allocate(ll(ng))
    rgrid = (boxsize) / real(ngrid)
    
    lirst = 0
    ll = 0
    
    do i = 1, ng
      indx = int((tracers(1, i)) / rgrid + 1.)
      indy = int((tracers(2, i)) / rgrid + 1.)
      indz = int((tracers(3, i)) / rgrid + 1.)
    
      if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
      indz.gt.0.and.indz.le.ngrid)lirst(indx,indy,indz)=i
    
      if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
      indz.gt.0.and.indz.le.ngrid)nlirst(indx,indy,indz) = &
      nlirst(indx, indy, indz) + 1
    end do
    
    do i = 1, ng
      indx = int((tracers(1, i))/ rgrid + 1.)
      indy = int((tracers(2, i))/ rgrid + 1.)
      indz = int((tracers(3, i))/ rgrid + 1.)
      if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
      &indz.gt.0.and.indz.le.ngrid) then
        ll(lirst(indx,indy,indz)) = i
        lirst(indx,indy,indz) = i
      endif
    end do
    
    write(*,*) 'Linked list successfully constructed'
    write(*,*) ''
    write(*,*) 'Starting loop over centres...'
    
    DD = 0
    cum_DD = 0
    delta = 0
    cum_delta = 0
    
    do i = 1, nc
      xvc = centres(1, i)
      yvc = centres(2, i)
      zvc = centres(3, i)
  
      ipx = int((xvc) / rgrid + 1.)
      ipy = int((yvc) / rgrid + 1.)
      ipz = int((zvc) / rgrid + 1.)
  
      ndif = int(perpmax / rgrid + 1.)
    
      do ix = ipx - ndif, ipx + ndif
        do iy = ipy - ndif, ipy + ndif
          do iz = ipz - ndif, ipz + ndif
    
            ix2 = ix
            iy2 = iy
            iz2 = iz
    
            if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
            if (ix2 .lt. 1) ix2 = ix2 + ngrid
            if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
            if (iy2 .lt. 1) iy2 = iy2 + ngrid
            if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
            if (iz2 .lt. 1) iz2 = iz2 + ngrid
    
            ii = lirst(ix2,iy2,iz2)
            if(ii.ne.0) then
              do
                ii = ll(ii)
                disx = tracers(1, ii) - xvc
                disy = tracers(2, ii) - yvc
                disz = tracers(3, ii) - zvc
  
                if (disx .lt. -boxsize/2) disx = disx + boxsize
                if (disx .gt. boxsize/2) disx = disx - boxsize
                if (disy .lt. -boxsize/2) disy = disy + boxsize
                if (disy .gt. boxsize/2) disy = disy - boxsize
                if (disz .lt. -boxsize/2) disz = disz + boxsize
                if (disz .gt. boxsize/2) disz = disz - boxsize
    
                r = (/ disx, disy, disz /)
                com = (/ 0, 0, 1 /)
                para = abs(dot_product(r, com))
                perp = sqrt(norm2(r) ** 2 - para ** 2)
  
                if (perp .gt. perpmin .and. perp .lt. perpmax .and. &
                    para .gt. paramin .and. para .lt. paramax ) then
                  perpind = int((perp - perpmin) / perpwidth + 1)
                  paraind = int((para - paramin) / parawidth + 1)
                  DD(i, perpind, paraind) = DD(i, perpind, paraind) + 1
                end if
  
    
                if(ii.eq.lirst(ix2,iy2,iz2)) exit
    
              end do
            end if
          end do
        end do
      end do
  
      do ii = 1, nperpbin
        do jj = 1, nparabin
          diff_vol = 2 * pi * (paramax - paramin) * (perpmax**2 - perpmin**2)
          delta(i, ii, jj) = DD(i, ii, jj) / (diff_vol * rhomed) - 1
        end do
      end do
    end do
    
    write(*,*) ''
    write(*,*) 'Calculation finished. Writing output...'
    
    open(12, file=output_den, status='replace', form='unformatted')
  
    write(12) nc
    write(12) size(perpbin)
    write(12) size(parabin)
    write(12) perpbin
    write(12) parabin
    write(12) delta
  
    end program density_profiles
    