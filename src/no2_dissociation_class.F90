!-----------------------------------------------------------------------------!
!  photolysis:                                                                !
!         NO2 + hv -> NO + O(3P)                                              !
!  Cross section from JPL94 (can also have Davidson et al.)                   !
!  Quantum yield from Gardiner, Sperry, and Calvert                           !
!-----------------------------------------------------------------------------!
module no2_dissociation_class
  use phot_kind_mod, only: rk => kind_phot
  use dissociation_class, only: dissociation_t

  implicit none

  type, extends(dissociation_t) :: no2_dissociation_t
     real(rk), allocatable :: ydel(:), yg(:)
   contains
     procedure :: initialize
     procedure :: update_rates
  end type no2_dissociation_t

  real(rk), parameter :: initval = -huge(1._rk)

contains

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine initialize(self, rxn_spec, errmsg, errflg)
    use phot_specs_mod,  only: phot_specs_t
    use wavelength_grid, only: nwave, wl
    use photo_utils,     only: base_read, add_pnts_inter2
    use params_mod,      only: deltax
    use module_xsections,only: rdxs_init

    class(no2_dissociation_t), intent(inout) :: self
    type(phot_specs_t),  intent(in) :: rxn_spec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: kdata = 250
    integer ::  n,  j
    real(rk) :: x1(kdata), x2(kdata)
    real(rk) :: y1(kdata), y2(kdata)
    real(rk) :: yg1(kdata),yg2(kdata)

    integer :: nsav, nw
    real(rk) :: xsav(kdata)

    call rdxs_init( nwave, wl, errmsg, errflg )

    allocate(self%yg(nwave))
    allocate(self%ydel(nwave))

    nw = nwave+1

    self%nchannels = rxn_spec%nchannels

    allocate(self%equations(self%nchannels))
    allocate(self%tag_names(self%nchannels))
    do j = 1,self%nchannels
       self%equations(j) = rxn_spec%equations(j)
       self%tag_names(j) = rxn_spec%tag_names(j)
    end do

    n = rxn_spec%nread(1)
    nsav = rxn_spec%nread(1)
    CALL base_read( filespec=rxn_spec%filename(1) , errmsg=errmsg, errflg=errflg, &
         skip_cnt=rxn_spec%nskip(1),rd_cnt=n,x=x1,y=y1,y1=y2 )
    xsav(1:n) = x1(1:n)
    CALL add_pnts_inter2(x1,y1,yg1,kdata,n, &
         nw,wl,self%equations(1),deltax,(/y1(1),0._rk/), errmsg, errflg)
    n = nsav
    x1(1:n) = xsav(1:n)
    CALL add_pnts_inter2(x1,y2,yg2,kdata,n, &
         nw,wl,self%equations(1),deltax,(/y2(1),0._rk/), errmsg, errflg)

    self%ydel(1:nwave) = yg1(1:nwave) - yg2(1:nwave)
    self%yg(1:nwave) = yg1(1:nwave)

    do j=1,self%nchannels
       write(*,*) ' Read in data for '//trim(self%equations(j))//' -- tag: '//trim(self%tag_names(j))
    end do

  end subroutine initialize

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine update_rates(self, temp, dens, fluxes, etf, esfact)
    use wavelength_grid, only: wl, nwave

    use module_xsections, only : no2xs_jpl06a

    class(no2_dissociation_t), intent(inout) :: self
    real(rk), intent(in) :: temp(:)
    real(rk), intent(in) :: dens(:)
    real(rk), intent(in) :: fluxes(:,:)
    real(rk), intent(in) :: etf(:)
    real(rk), intent(in) :: esfact

    real(rk) :: t(size(temp))
    real(rk) :: no2xs(size(temp),nwave), qy(nwave), sq(size(temp),nwave)
    real(rk) :: xsect(nwave)
    integer :: nz, iw, iz

    nz = size(temp)

    if (.not.associated(self%rateconsts)) then
       allocate(self%rateconsts(nz,self%nchannels))
    end if

    self%rateconsts = initval

    call no2xs_jpl06a(nz,temp,nwave+1,wl, no2xs)

    t(1:nz) = .02_rk*(temp(1:nz) - 298._rk)
    do iw = 1, nwave
       qy(1:nz) = self%yg(iw) + self%ydel(iw)*t(1:nz)
       sq(1:nz,iw) = no2xs(1:nz,iw)*max( qy(1:nz),0._rk )
    enddo

    do iz = 1, nz
       xsect(1:nwave) = sq(iz,1:nwave)*etf(1:nwave)*esfact
       self%rateconsts(iz,1) = dot_product( fluxes(1:nwave,iz),xsect(1:nwave) )
    end do

  end subroutine update_rates

end module no2_dissociation_class
