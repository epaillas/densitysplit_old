program density_profiles
  implicit none
  
  integer, parameter:: dp=kind(0.d0)
  
  real(dp) :: rgrid, boxsize, diff_vol, cum_vol, rhomed
  real(dp) :: disx, disy, disz, dis, vr
  real(dp) :: xvc, yvc, zvc
  real(dp) :: velx, vely, velz, comx, comy, comz
  real(dp) :: rwidth, dmax, dmin
  real(dp) :: pi = 4.*atan(1.)
  
  integer*4 :: ng, nc, nrbin, rind
  integer*4 :: id, iargc
  integer*4 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*4 :: indx, indy, indz, nrows, ncols
  integer*4 :: ipx, ipy, ipz, ndif
  integer*4 :: ngrid
  
  integer*4, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*4, dimension(:), allocatable :: ll
  
  real(dp), dimension(3) :: r, v
  real(dp), allocatable, dimension(:,:)  :: tracers, centres
  real(dp), dimension(:, :), allocatable :: DD, RR, cum_DD, delta, cum_delta
  real(dp), dimension(:, :), allocatable :: radial_vel, radial_vel2, mean_vr, std_vr
  real(dp), dimension(:), allocatable :: rbin, rbin_edges
  
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
  write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)

  open(11, file=input_centres, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(centres(ncols, nrows))
  read(11) centres
  close(11)
  nc = nrows
  write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)

  print*, maxval(centres(1, :))
  print*, minval(centres(1, :))

  allocate(rbin(nrbin))
  allocate(rbin_edges(nrbin+1))
  allocate(DD(nc, nrbin))
  allocate(cum_DD(nc, nrbin))
  allocate(RR(nc, nrbin))
  allocate(delta(nc, nrbin))
  allocate(cum_delta(nc, nrbin))
  allocate(radial_vel(nc, nrbin))
  allocate(radial_vel2(nc, nrbin))
  allocate(mean_vr(nc, nrbin))
  allocate(std_vr(nc, nrbin))
  
  
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
  radial_vel = 0
  radial_vel2 = 0
  mean_vr = 0
  std_vr = 0
  
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

              velx = tracers(4, ii)
              vely = tracers(5, ii)
              velz = tracers(6, ii)

              comx = (xvc + tracers(1, ii)) / 2
              comy = (yvc + tracers(2, ii)) / 2
              comz = (zvc + tracers(3, ii)) / 2
  
              if (comx .lt. -boxsize/2) comx = comx + boxsize
              if (comx .gt. boxsize/2) comx = comx - boxsize
              if (comy .lt. -boxsize/2) comy = comy + boxsize
              if (comy .gt. boxsize/2) comy = comy - boxsize
              if (comz .lt. -boxsize/2) comz = comz + boxsize
              if (comz .gt. boxsize/2) comz = comz - boxsize
  
              if (disx .lt. -boxsize/2) disx = disx + boxsize
              if (disx .gt. boxsize/2) disx = disx - boxsize
              if (disy .lt. -boxsize/2) disy = disy + boxsize
              if (disy .gt. boxsize/2) disy = disy - boxsize
              if (disz .lt. -boxsize/2) disz = disz + boxsize
              if (disz .gt. boxsize/2) disz = disz - boxsize
  
              r = (/ disx, disy, disz /)
              v = (/ velx, vely, velz /)

              dis = norm2(r)
              vr = dot_product(v, r) / norm2(r)

              if (dis .gt. dmin .and. dis .lt. dmax) then
                rind = int((dis - dmin) / rwidth + 1)
                DD(i, rind) = DD(i, rind) + 1
                radial_vel(i, rind) = radial_vel(i, rind) + vr
                radial_vel2(i, rind) = radial_vel2(i, rind) + vr ** 2
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
      if (DD(i, ii) .ne. 0) then
        mean_vr(i, ii) = radial_vel(i, ii) / DD(i, ii)
        std_vr(i, ii) = sqrt((radial_vel2(i, ii) - (radial_vel(i, ii) ** 2 / DD(i, ii))) / (DD(i, ii) - 1))
      else
        mean_vr(i, ii) = 0
        std_vr(i, ii) = 0
      end if

    end do
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_den, status='replace', form='unformatted')

  write(12) nc
  write(12) size(rbin, dim=1)
  write(12) rbin
  write(12) DD
  write(12) delta
  write(12) cum_delta
  write(12) mean_vr

  end program density_profiles
  