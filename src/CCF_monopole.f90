program density_profiles
  implicit none
  
  real*4 :: rgrid, boxsize, diff_vol, cum_vol, rhomed
  real*4 :: disx, disy, disz, dis, vr, vlos
  real*4 :: xvc, yvc, zvc
  real*4 :: velx, vely, velz
  real*4 :: rwidth, dmax, dmin
  real*4 :: pi = 4.*atan(1.)
  
  integer*4 :: ng, nc, nrbin, rind
  integer*4 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*4 :: indx, indy, indz, nrows, ncols
  integer*4 :: ipx, ipy, ipz, ndif
  integer*4 :: ngrid
  
  integer*4, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*4, dimension(:), allocatable :: ll
  
  real*4, dimension(3) :: r, vel, com
  real*4, allocatable, dimension(:,:)  :: tracers, centres
  real*4, dimension(:, :), allocatable :: DD, cum_DD, delta, cum_delta
  real*4, dimension(:, :), allocatable :: VV_r, VV_los, VV2_los, mean_vr, std_vlos
  real*4, dimension(:), allocatable :: rbin, rbin_edges

  logical :: has_velocity = .false.
  
  character(20), external :: str
  character(len=500) :: input_tracers, input_centres, output_den
  character(len=10) :: dmax_char, dmin_char, nrbin_char, ngrid_char, box_char
  
  if (iargc() .ne. 8) then
      write(*,*) 'Some arguments are missing.'
      write(*,*) '1) input_data'
      write(*,*) '2) input_centres'
      write(*,*) '3) output_den'
      write(*,*) '4) boxsize'
      write(*,*) '5) dmin'
      write(*,*) '6) dmax'
      write(*,*) '7) nrbin'
      write(*,*) '8) ngrid'
      write(*,*) ''
      stop
    end if
    
  call getarg(1, input_tracers)
  call getarg(2, input_centres)
  call getarg(3, output_den)
  call getarg(4, box_char)
  call getarg(5, dmin_char)
  call getarg(6, dmax_char)
  call getarg(7, nrbin_char)
  call getarg(8, ngrid_char)
  
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
  write(*, *) 'input_centres: ', trim(input_centres)
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

  open(11, file=input_centres, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(centres(ncols, nrows))
  read(11) centres
  close(11)
  nc = nrows
  write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)

  allocate(rbin(nrbin))
  allocate(rbin_edges(nrbin + 1))
  allocate(DD(nc, nrbin))
  allocate(cum_DD(nc, nrbin))
  allocate(delta(nc, nrbin))
  allocate(cum_delta(nc, nrbin))
  if (has_velocity) then
    allocate(VV_r(nc, nrbin))
    allocate(VV_los(nc, nrbin))
    allocate(VV2_los(nc, nrbin))
    allocate(mean_vr(nc, nrbin))
    allocate(std_vlos(nc, nrbin))
  end if
  
  
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
  cum_DD = 0
  delta = 0
  cum_delta = 0
  if (has_velocity) then
    VV_r = 0
    VV_los = 0
    VV2_los = 0
    mean_vr = 0
    std_vlos = 0
  end if
  
  do i = 1, nc
    xvc = centres(1, i)
    yvc = centres(2, i)
    zvc = centres(3, i)

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

              if (has_velocity) then
                velx = tracers(4, ii)
                vely = tracers(5, ii)
                velz = tracers(6, ii)
                vel = (/ velx, vely, velz /)
                com = (/ 0, 0, 1 /)
                vr = dot_product(vel, r) / norm2(r)
                vlos = dot_product(vel, com)
              end if

              if (dis .gt. dmin .and. dis .lt. dmax) then
                rind = int((dis - dmin) / rwidth + 1)
                DD(i, rind) = DD(i, rind) + 1

                if (has_velocity) then
                  VV_r(i, rind) = VV_r(i, rind) + vr
                  VV_los(i, rind) = VV_los(i, rind) + vlos
                  VV2_los(i, rind) = VV2_los(i, rind) + vlos**2
                end if
              end if

  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do


    cum_DD(i, 1) = DD(i, 1)
    do ii = 2, nrbin
      cum_DD(i, ii) = cum_DD(i, ii - 1) + DD(i, ii)
    end do
  
    do ii = 1, nrbin

      diff_vol = 4./3 * pi * (rbin_edges(ii+1)**3 - rbin_edges(ii)**3)
      cum_vol = 4./3 * pi * rbin_edges(ii+1)**3

      delta(i, ii) = DD(i, ii) / (diff_vol * rhomed) - 1
      cum_delta(i, ii) = cum_DD(i, ii) / (cum_vol * rhomed) - 1

      ! if (has_velocity) then
      !   if (DD(i, ii) .gt. 1) then
      !     mean_vr(i, ii) = VV_r(i, ii) / DD(i, ii)
      !     std_vlos(i, ii) = sqrt((VV2_los(i, ii) - (VV_los(i, ii) ** 2 / DD(i, ii))) / (DD(i, ii) - 1))
      !   else if (DD(i, ii) .eq. 1) then
      !     mean_vr(i, ii) = VV_r(i, ii) / DD(i, ii)
      !     std_vlos(i, ii) = 0
      !   else
      !     mean_vr(i, ii) = 0
      !     std_vlos(i, ii) = 0
      !   end if
      ! end if

      if (has_velocity) then
        mean_vr(i, ii) = VV_r(i, ii) / DD(i, ii)
        std_vlos(i, ii) = sqrt((VV2_los(i, ii) - (VV_los(i, ii) ** 2 / DD(i, ii))) / (DD(i, ii) - 1))
      end if

    end do
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_den, status='replace', form='unformatted')

  write(12) nc
  write(12) size(rbin)
  write(12) rbin
  write(12) DD
  write(12) delta
  write(12) cum_delta
  if (has_velocity) then
    write(12) mean_vr
    write(12) std_vlos
  end if

  end program density_profiles
  