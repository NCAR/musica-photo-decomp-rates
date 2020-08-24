program test
  use phot_kind_mod, only: rk => kind_phot
  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run
  use tuv_photolysis, only: tuv_photolysis_readnl, tuv_photolysis_init
  use wavelength_grid,only: n_wavelen=>nwave
  use molec_ox_xsect, only: molec_ox_xsect_init
  use molec_ox_xsect, only: molec_ox_xsect_run

  implicit none

  integer :: errflg
  character(len=444) :: errmsg
  character(len=*), parameter :: nml_file = '/terminator-data1/home/fvitt/photodev/no-atmos-test/test/integration/data/test_nml'

  integer, parameter :: nlevels = 20

  real(rk) :: dz, hz, press_top
  real(rk) :: zenith
  real(rk) :: albedo
  real(rk) :: alt(nlevels)
  real(rk) :: press(nlevels)
  real(rk) :: temp(nlevels)
  real(rk) :: o2vmr(nlevels)
  real(rk) :: o3vmr(nlevels)
  real(rk) :: so2vmr(nlevels)
  real(rk) :: no2vmr(nlevels)
  real(rk) :: cldfrac(nlevels)
  real(rk) :: cldwat(nlevels)
  
  real(rk), allocatable :: srb_o2_xs(:,:)
  real(rk), allocatable :: dto2(:,:)
  real(rk), allocatable :: radfld(:,:)

  integer :: k, wn
  logical :: test_pass

  write(*,*) 'Test radiation transfer with no absorbers and no reflection'

  o2vmr = 0._rk
  o3vmr = 0._rk
  so2vmr = 0._rk
  no2vmr = 0._rk
  cldfrac = 0._rk
  cldwat = 0._rk

  temp = 300._rk
  albedo = 0._rk
  zenith = 0._rk
  
  dz = 10000._rk
  hz = 8000._rk
  alt(nlevels) = dz*0.5_rk
  press(nlevels) = 1.e5_rk * exp(-alt(nlevels)/hz)

  do k = nlevels-1,1, -1
     alt(k) = alt(k+1)+dz
     press(k) = 1.e5_rk * exp(-alt(k)/hz)
  end do

  press_top = 1.e5_rk * exp(-(alt(1)+.5_rk*dz)/hz)
    
  call tuv_photolysis_readnl(nml_file, errmsg, errflg)
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call molec_ox_xsect_init( errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call tuv_radiation_transfer_init( rk, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  allocate(srb_o2_xs(n_wavelen,nlevels), dto2(nlevels,n_wavelen))
  allocate(radfld(n_wavelen,nlevels))

  call molec_ox_xsect_run( nlevels, zenith, alt, temp, press, press_top, o2vmr, dto2, srb_o2_xs, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if
  
  call tuv_radiation_transfer_run( nlevels, n_wavelen, zenith, albedo, press, press_top, alt, temp, o3vmr, &
                                   so2vmr, no2vmr, cldfrac, cldwat, dto2, radfld, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  test_pass = .true.
  ! radfld is normalized so without any absorption or any reflection the actinic fluxes should be 1
  do wn = 1,n_wavelen
     write(*,*) radfld(wn,1), abs(radfld(wn,1)-1._rk)
     if ( abs(radfld(wn,1)-1._rk) > 1.e-3_rk ) then
        test_pass = .false.
     endif
  enddo

  if( .not. test_pass ) then
     write(*,*) "Test failed!"
     stop 3
  end if

  write(*,*) 'Test completed successfully'

end program test
