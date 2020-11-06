!------------------------------------------------------------------------------
! generic dissociation class
!------------------------------------------------------------------------------
module dissociation_class
  use phot_kind_mod, only: rk => kind_phot

  implicit none

  type :: dissociation_t
     integer :: nchannels
     character(len=80), pointer :: equations(:)=>null()
     character(len=16), pointer :: tag_names(:)=>null()
     real(rk), pointer :: xsqy(:,:)=>null()
     real(rk), pointer :: rateconsts(:,:)=>null()
   contains
     procedure :: initialize => diss_initialize
     procedure :: update_rates => diss_update_rates
  end type dissociation_t

  real(rk), parameter :: initval = -huge(1._rk)

contains

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine diss_initialize(self, rxn_spec, errmsg, errflg)
    use phot_specs_mod,  only: phot_specs_t
    use wavelength_grid, only: nw=>nwave, wl
    use photo_utils,     only: base_read, add_pnts_inter2
    use params_mod,      only: deltax

    class(dissociation_t), intent(inout) :: self
    type(phot_specs_t),  intent(in) :: rxn_spec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: j

    integer, PARAMETER :: kdata=1000
    real(rk) :: yg(kdata), y1(kdata), x1(kdata)
    integer :: nread

    self%nchannels = rxn_spec%nchannels

    allocate(self%equations(self%nchannels))
    allocate(self%tag_names(self%nchannels))
    do j = 1,self%nchannels
       self%equations(j) = rxn_spec%equations(j)
       self%tag_names(j) = rxn_spec%tag_names(j)
    end do

    allocate(self%xsqy(nw,self%nchannels))
    self%xsqy(:,:) = initval

    nread = rxn_spec%nread(1)

    if( rxn_spec%nskip(1)>0 ) then
       CALL base_read( filespec=rxn_spec%filename(1), errmsg=errmsg, errflg=errflg, &
                       skip_cnt=rxn_spec%nskip(1), rd_cnt=nread, x=x1, y=y1 )
    else
       CALL base_read( filespec=rxn_spec%filename(1), errmsg=errmsg, errflg=errflg, &
                       rd_cnt=nread, x=x1, y=y1 )
    endif
    y1(1:nread) = y1(1:nread) * rxn_spec%xfac(1)

    CALL add_pnts_inter2(x1,y1,yg,kdata,nread, nw+1, wl,self%equations(1),deltax,(/0._rk,0._rk/), &
         errmsg, errflg)

    do j=1,self%nchannels
       self%xsqy(:nw,j) = yg(:nw)*rxn_spec%quantum_yields(j)
       write(*,*) ' Read in data for '//trim(self%equations(j))//' -- tag: '//trim(self%tag_names(j))
    end do

  end subroutine diss_initialize

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine diss_update_rates(self, temp, dens, fluxes, etf, esfact)
    use wavelength_grid, only: nwave

    class(dissociation_t), intent(inout) :: self
    real(rk), intent(in) :: temp(:)
    real(rk), intent(in) :: dens(:)
    real(rk), intent(in) :: fluxes(:,:)
    real(rk), intent(in) :: etf(:)
    real(rk), intent(in) :: esfact

    integer :: j, nlev
    real(rk) :: xsect(nwave)

    nlev=size(fluxes,dim=2)

    if (.not.associated(self%rateconsts)) then
       allocate(self%rateconsts(nlev,self%nchannels))
    end if

    do j = 1,self%nchannels
       xsect(1:nwave) = self%xsqy(1:nwave,j)*etf(1:nwave)*esfact
       self%rateconsts(:,j) = matmul( transpose( fluxes ), xsect(1:nwave) )
    end do

  end subroutine diss_update_rates

end module dissociation_class
