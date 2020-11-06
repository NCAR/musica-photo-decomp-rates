program test_jrates
  use phot_kind_mod, only: rk => kind_phot
  use jrates_mod, only: jrates_init, jrates_calc
  use wavelength_grid, only: wavelength_grid_init, nwave
  use params_mod, only: input_data_root
  use etf_mod, only: read_etf, etf=>photon_flux

  use wavelength_grid, only: wavelength_grid_init, nwave
  use environ_conditions_mod, only: environ_conditions_create, environ_conditions

  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run
  use molec_ox_xsect, only: molec_ox_xsect_init
  use molec_ox_xsect, only: molec_ox_xsect_run

  use params_mod, only: input_data_root
  use etf_mod, only: read_etf, etf=>photon_flux

  implicit none

  character(len=200) :: errmsg
  integer            :: errflg
  integer, parameter :: nrxns = 15
  character(len=80)  :: rxn_eqns(nrxns)
  character(len=8), parameter :: jnames(nrxns) = (/ &
       'jno2    ', &
       'jo3a    ', &
       'jo3b    ', &
       'jch4a   ', &
       'jch4b   ', &
       'jco2    ', &
       'jcof2   ', &
       'jcofcl  ', &
       'jhbr    ', &
       'jhf     ', &
       'jsf6    ', &
       'jglyald ', &
       'jhyac   ', &
       'jmacr_a ', &
       'jmacr_b ' &
       /)

  type(environ_conditions),pointer :: colEnvConds => null()
  integer :: nlevels
  real(rk) :: zenith
  real(rk) :: albedo
  real(rk) :: esfact
  real(rk), allocatable :: alt(:)
  real(rk), allocatable :: press_mid(:)
  real(rk), allocatable :: press_int(:)
  real(rk), allocatable :: temp(:)
  real(rk), allocatable :: dens(:)

  real(rk), allocatable :: o2vmrcol(:)
  real(rk), allocatable :: o3vmrcol(:)
  real(rk), allocatable :: so2vmrcol(:)
  real(rk), allocatable :: no2vmrcol(:)

  real(rk), allocatable :: cldfrac(:)
  real(rk), allocatable :: cldwat(:)
  real(rk), allocatable :: srb_o2_xs(:,:)
  real(rk), allocatable :: dto2(:,:)
  real(rk), allocatable :: radfld(:,:)
  real(rk), allocatable :: rate_consts(:,:)

  real(rk), parameter :: kboltz = 1.38064852e-16_rk ! boltzmann constant (erg/K)
  integer :: j,k
  character(len=300) :: etf_wave_grid_filepath
  character(len=80) :: wavelength_grid_file

  character(len=128) :: env_conds_file='NONE'
  real :: env_conds_lat=45.
  real :: env_conds_lon=180.
  character(len=*), parameter :: nml_file = 'data/test_nml'
  integer :: unitn

  namelist /drv_opts/ env_conds_file, env_conds_lat, env_conds_lon
  namelist /drv_opts/ input_data_root, wavelength_grid_file

  write(*,*) 'BEGIN TEST...'
  open(newunit=unitn, file=trim(nml_file), status='old')
  read(unit=unitn, nml=drv_opts)
  close(unitn)

  etf_wave_grid_filepath = trim(input_data_root)//'/'//trim(wavelength_grid_file)

  ! init wavel length grid
  call wavelength_grid_init( etf_wave_grid_filepath, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  call jrates_init( jnames, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  ! load env conds
  colEnvConds => environ_conditions_create( env_conds_file, lat=env_conds_lat, lon=env_conds_lon )
  nlevels = colEnvConds%nlayers()

  !rad xfer ...

  call molec_ox_xsect_init( errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  call tuv_radiation_transfer_init( rk, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  allocate(srb_o2_xs(nwave,nlevels), dto2(nlevels,nwave))

  allocate(radfld(nwave,nlevels))

  allocate(alt(nlevels))
  allocate(press_mid(nlevels))
  allocate(press_int(nlevels+1))
  allocate(temp(nlevels))
  allocate(dens(nlevels))

  allocate(o2vmrcol(nlevels))
  allocate(o3vmrcol(nlevels))
  allocate(so2vmrcol(nlevels))
  allocate(no2vmrcol(nlevels))

  allocate(cldfrac(nlevels))
  allocate(cldwat(nlevels))


  zenith = colEnvConds%getsrf('SZA')

  albedo = colEnvConds%getsrf('ASDIR')
  press_mid(:nlevels) = colEnvConds%press_mid(nlevels)
  press_int(:nlevels+1) = colEnvConds%press_int(nlevels+1)
  alt(:nlevels) = colEnvConds%getcol('Z3',nlevels) ! meters
  temp(:nlevels) = colEnvConds%getcol('T',nlevels)
  dens(:) = 10._rk*press_mid(:)/(kboltz*temp(:)) ! # molecules / cm3 in each layer

  o2vmrcol(:nlevels) = colEnvConds%getcol('O2',nlevels)
  o3vmrcol(:nlevels) = colEnvConds%getcol('O3',nlevels)
  so2vmrcol(:nlevels) = colEnvConds%getcol('SO2',nlevels)
  no2vmrcol(:nlevels) = colEnvConds%getcol('NO2',nlevels)

  cldfrac = 0._rk
  cldwat = 0._rk

  call molec_ox_xsect_run( nlevels, zenith, alt, temp, press_mid, press_int(1), o2vmrcol, dto2, srb_o2_xs, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  call tuv_radiation_transfer_run( nlevels, nwave, zenith, albedo, press_mid, press_int(1), alt, temp, o3vmrcol, &
                                   so2vmrcol, no2vmrcol, cldfrac, cldwat, dto2, radfld, errmsg, errflg )
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  call read_etf( etf_wave_grid_filepath, errmsg, errflg)
  if (errflg/=0) then
     write(*,*) 'ERROR: ',errflg,trim(errmsg)
     stop errflg
  endif

  esfact = 1._rk

  allocate(rate_consts(nlevels,nrxns))
  rate_consts = -huge(1._rk)

  ! invert radfld to be top-down vert coord
  radfld(1:nwave,1:nlevels) = radfld(1:nwave,nlevels:1:-1)
  call jrates_calc( radfld, etf, esfact, temp, dens, rate_consts, rxn_eqns=rxn_eqns )

  write(*,*) " Env conds file: "//trim(env_conds_file)
  write(*,*) " Env conds lon : ",env_conds_lon
  write(*,*) " Env conds lat : ",env_conds_lat
  write(*,*) ' zen angle: ',zenith

  write(10,*) " Env conds file: "//trim(env_conds_file)
  write(10,*) " Env conds lon : ",env_conds_lon
  write(10,*) " Env conds lat : ",env_conds_lat
  write(10,*) ' zen angle: ',zenith

  do j = 1,nrxns

     write( *,*) trim(jnames(j))//'   '//trim(rxn_eqns(j))
     write(10,*) trim(jnames(j))//'   '//trim(rxn_eqns(j))
     write( *,'("  rate = ",e12.4," /sec")' ) rate_consts(nlevels,j)
     do k=1,nlevels
         write(10,'("  rate = ",e24.16," /sec")' ) rate_consts(k,j)
     end do

     write(10,*) ' '
     write(*,*)
  end do

  write(*,*) 'END TEST'

end program test_jrates
