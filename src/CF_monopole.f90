program CF_monopole
    implicit none
    
    real*8 :: rgrid, boxsize, vol, int_vol, rhomed
    real*8 :: disx, disy, disz, dis, vr, vlos
    real*8 :: xvc, yvc, zvc
    real*8 :: velx, vely, velz
    real*8 :: rwidth, dmax, dmin
    real*8 :: pi = 4.*atan(1.)
    
    integer*8 :: ng, nc, nrbin, rind
    integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: indx, indy, indz, nrows, ncols
    integer*8 :: ipx, ipy, ipz, ndif
    integer*8 :: ngrid
    
    integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
    integer*8, dimension(:), allocatable :: ll
    
    real*8, dimension(3) :: r, vel, com
    real*8, allocatable, dimension(:,:)  :: tracers
    real*8, dimension(:), allocatable :: DD, int_DD, delta, int_delta
    real*8, dimension(:), allocatable :: VV_r, VV_los, VV2_los, mean_vr, std_vlos
    real*8, dimension(:), allocatable :: rbin, rbin_edges
  
    logical :: has_velocity = .false.
    
    character(20), external :: str
    character(len=500) :: input_tracers, output_den
    character(len=10) :: dmax_char, dmin_char, nrbin_char, ngrid_char, box_char
    
    if (iargc() .ne. 7) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) input_tracers'
        write(*,*) '2) output_den'
        write(*,*) '3) boxsize'
        write(*,*) '4) dmin'
        write(*,*) '5) dmax'
        write(*,*) '6) nrbin'
        write(*,*) '7) ngrid'
        write(*,*) ''
        stop
      end if
      
    call getarg(1, input_tracers)
    call getarg(2, output_den)
    call getarg(3, box_char)
    call getarg(4, dmin_char)
    call getarg(5, dmax_char)
    call getarg(6, nrbin_char)
    call getarg(7, ngrid_char)
    
    read(box_char, *) boxsize
    read(dmin_char, *) dmin
    read(dmax_char, *) dmax
    read(nrbin_char, *) nrbin
    read(ngrid_char, *) ngrid
    
    write(*,*) '-----------------------'
    write(*,*) 'Running density_profiles.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'input_tracers: ', trim(input_tracers)
    write(*, *) 'boxsize: ', trim(box_char)
    write(*, *) 'output_den: ', trim(output_den)
    write(*, *) 'dmin: ', trim(dmin_char), ' Mpc'
    write(*, *) 'dmax: ', trim(dmax_char), ' Mpc'
    write(*, *) 'nrbin: ', trim(nrbin_char)
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
  
    allocate(rbin(nrbin))
    allocate(rbin_edges(nrbin + 1))
    allocate(DD(nrbin))
    allocate(int_DD(nrbin))
    allocate(delta(nrbin))
    allocate(int_delta(nrbin))
    ! if (has_velocity) then
    !   allocate(VV_r(nrbin))
    !   allocate(VV_los(nrbin))
    !   allocate(VV2_los(nrbin))
    !   allocate(mean_vr(nrbin))
    !   allocate(std_vlos(nrbin))
    ! end if
    
    
    rwidth = (dmax - dmin) / nrbin
    do i = 1, nrbin + 1
      rbin_edges(i) = dmin+(i-1)*rwidth
    end do
    do i = 1, nrbin
      rbin(i) = rbin_edges(i+1)-rwidth/2.
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
    int_DD = 0
    delta = 0
    int_DD = 0
    ! if (has_velocity) then
    !   VV_r = 0
    !   VV_los = 0
    !   VV2_los = 0
    !   mean_vr = 0
    !   std_vlos = 0
    ! end if
    
    do i = 1, ng
      xvc = tracers(1, i)
      yvc = tracers(2, i)
      zvc = tracers(3, i)
  
      ipx = int((xvc) / rgrid + 1.)
      ipy = int((yvc) / rgrid + 1.)
      ipz = int((zvc) / rgrid + 1.)
  
      ndif = int(dmax / rgrid + 1.)
    
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
                if (ii .eq. i) cycle

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
                dis = norm2(r)
  
                ! if (has_velocity) then
                !   velx = tracers(4, ii)
                !   vely = tracers(5, ii)
                !   velz = tracers(6, ii)
                !   vel = (/ velx, vely, velz /)
                !   com = (/ 0, 0, 1 /)
                !   vr = dot_product(vel, r) / norm2(r)
                !   vlos = dot_product(vel, com)
                ! end if
  
                if (dis .gt. dmin .and. dis .lt. dmax) then
                  rind = int((dis - dmin) / rwidth + 1)
                  DD(rind) = DD(rind) + 1
  
                !   if (has_velocity) then
                !     VV_r(rind) = VV_r(rind) + vr
                !     VV_los(rind) = VV_los(rind) + vlos
                !     VV2_los(rind) = VV2_los(rind) + vlos**2
                !   end if
                end if
    
                if(ii.eq.lirst(ix2,iy2,iz2)) exit
    
              end do
            end if
          end do
        end do
      end do
    end do
  
    int_DD(1) = DD(1)
    do i = 2, nrbin
      int_DD(i) =  int_DD(i - 1) + DD(i)
    end do
  
    do i = 1, nrbin
      vol = 4./3 * pi * (rbin_edges(i + 1) ** 3 - rbin_edges(i) ** 3)
      int_vol = 4./3 * pi * rbin_edges(i + 1) ** 3
      delta(i) = DD(i) / (2 * vol * rhomed * nc) - 1
      int_delta(i) = int_DD(i) / (2 * int_vol * rhomed * nc) - 1
  
    !   if (has_velocity) then
    !     mean_vr(i) = VV_r(i) / DD(i)
    !     std_vlos(i) = sqrt((VV2_los(i) - (VV_los(i) ** 2 / DD(i))) / (DD(i) - 1))
    !   end if
  
    end do
    
    write(*,*) ''
    write(*,*) 'Calculation finished. Writing output...'
    
    open(12, file=output_den, status='replace')
    do i = 1, nrbin
      write(12, fmt='(3f15.5)') rbin(i), delta(i), int_delta(i)
    end do
  
    end program CF_monopole
    