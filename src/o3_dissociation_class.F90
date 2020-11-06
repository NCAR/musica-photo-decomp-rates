!-----------------------------------------------------------------------------!
!   o3 photolysis reactions:                                                  !
!              (a) O3 + hv -> O2 + O(1D)         -- jo3_a                     !
!              (b) O3 + hv -> O2 + O(3P)         -- jo3_b                     !
!   cross sections:                                                           !
!               120nm - 185 nm = Ackerman, 1971                               !
!               185nm - 827 nm = JPL06 (293-298K)                             !
!               196nm - 342 nm = JPL06 (218 K)                                !
!   quantum yield:  JPL 06 recommendation                                     !
!-----------------------------------------------------------------------------!
module o3_dissociation_class
  use phot_kind_mod, only: rk => kind_phot
  use dissociation_class, only: dissociation_t

  implicit none

  type, extends(dissociation_t) :: o3_dissociation_t
     real(rk), allocatable :: yg298(:), yg218(:), xso3(:)
   contains
     procedure :: initialize
     procedure :: update_rates
  end type o3_dissociation_t

  real(rk), parameter :: initval = -huge(1._rk)

contains

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine initialize(self, rxn_spec, errmsg, errflg)
    use phot_specs_mod,  only: phot_specs_t
    use wavelength_grid, only: nwave, wl, wc
    use photo_utils,     only: base_read, add_pnts_inter2
    use params_mod,      only: deltax

    class(o3_dissociation_t), intent(inout) :: self
    type(phot_specs_t),  intent(in) :: rxn_spec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: kdata = 250
    integer :: i, iw, iz, n, idum, kk, j
    real(rk) :: x1(kdata), x2(kdata), x3(kdata)
    real(rk) :: y1(kdata), y2(kdata), y3(kdata)
    real(rk) :: yg(nwave)

    integer :: nw

    self%nchannels = rxn_spec%nchannels

    allocate(self%equations(self%nchannels))
    allocate(self%tag_names(self%nchannels))
    do j = 1,self%nchannels
       self%equations(j) = rxn_spec%equations(j)
       self%tag_names(j) = rxn_spec%tag_names(j)
    end do


    nw = nwave+1
    allocate( self%yg298(nwave), self%yg218(nwave), self%xso3(nwave) )

    !-----------------------------------------------------------
    !     Ref: Atkinson, ultraviolet solar radiation related to
    !          mesopheric procesess, 149-159, in fiocco, g. (ed.),
    !          mesospheric models and related exp., d. reidel,
    !          dordrecht, 1971.
    !
    !          120 nm through 200 nm
    !-----------------------------------------------------------
    x1 = initval
    y1 = initval
    n = rxn_spec%nread(1)
    call base_read( filespec=rxn_spec%filename(1), &
         errmsg=errmsg, errflg=errflg, skip_cnt=rxn_spec%nskip(1), &
         rd_cnt=n, x=x1,y=y1 )
    call add_pnts_inter2(x1,y1,yg, kdata, n, &
         nw,wl,'O3 + hv ->',deltax,(/0._rk,0._rk/), errmsg, errflg)

    do iw = 1, nwave
       self%xso3(iw) = yg(iw)
    enddo

    !-------------------------------------------------------
    !     ... REF: JPL06 218K
    !         from 196.078 to 342.5 nm
    !-------------------------------------------------------
    x1 = initval
    y1 = initval
    y2 = initval
    n = rxn_spec%nread(2)
    call base_read( filespec=rxn_spec%filename(2), &
         errmsg=errmsg, errflg=errflg, &
         skip_cnt=rxn_spec%nskip(2), &
         rd_cnt=n, &
         x=x1,y=y1, y1=y2)
    x2(:n) = 0.5_rk*(x1(:n) + y1(:n))
    call add_pnts_inter2(x2,y2,self%yg218, kdata, n, &
         nw,wl, 'O3 + hv ->',deltax,(/0._rk,0._rk/), errmsg, errflg)

    !-------------------------------------------------------
    !     ... REF: JPL06 293-298K
    !         from 185.185 to 827.500 nm
    !-------------------------------------------------------
    x1 = initval
    y1 = initval
    y3 = initval
    n = rxn_spec%nread(3)
    call base_read( filespec=rxn_spec%filename(3), &
         errmsg=errmsg, errflg=errflg, &
         skip_cnt=rxn_spec%nskip(3), &
         rd_cnt=n, &
         x=x1,y=y1, y1=y3)
    x3(:n) = 0.5_rk*(x1(:n) + y1(:n))
    call add_pnts_inter2(x3,y3,self%yg298, kdata, n, &
         nw,wl,'O3 + hv ->',deltax,(/0._rk,0._rk/), errmsg, errflg)

    do iw = 1, nwave
       if (wc(iw) .ge. 184.0_rk) then
          self%xso3(iw) = self%yg298(iw)
       endif
    enddo

    do j=1,self%nchannels
       write(*,*) ' Read in data for '//trim(self%equations(j))//' -- tag: '//trim(self%tag_names(j))
    end do

  end subroutine initialize

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  subroutine update_rates(self, temp, dens, fluxes, etf, esfact)
    use wavelength_grid, only: wc, nwave

    class(o3_dissociation_t), intent(inout) :: self
    real(rk), intent(in) :: temp(:)
    real(rk), intent(in) :: dens(:)
    real(rk), intent(in) :: fluxes(:,:)
    real(rk), intent(in) :: etf(:)
    real(rk), intent(in) :: esfact

    integer :: j, iw, iz
    integer :: nz

    real(rk), parameter :: A(3) = (/0.8036_rk, 8.9061_rk, 0.1192_rk /)
    real(rk), parameter :: X(3) = (/ 304.225_rk, 314.957_rk, 310.737_rk/)
    real(rk), parameter :: om(3) = (/ 5.576_rk, 6.601_rk, 2.187_rk/)

    real(rk) :: QY_O1D(size(temp),nwave) ! nz, nw

    real(rk) :: qy(size(temp),nwave) ! nz, nw
    real(rk) :: T, kt, q1, q2
    real(rk) :: so3(size(temp),nwave)

    real(rk) :: xsqy(nwave)
    real(rk) :: xsect(nwave)


    nz = size(temp)

    if (.not.associated(self%rateconsts)) then
       allocate(self%rateconsts(nz,self%nchannels))
    end if

    self%rateconsts = initval

    !-------------------------------------------------------
    !     ... for hartley and huggins bands, use
    !         temperature-dependent values from
    !         JPL06
    !-------------------------------------------------------
    !-------------------------------------------------------
    !     ... Cross Sections and Quantum Yields
    !-------------------------------------------------------
    do iw = 1, nwave
       do iz = 1, nz

          so3(iz,iw) = self%xso3(iw)

          if ((wc(iw) .ge. 196.078_rk) .and. (wc(iw) .le. 342.5_rk)) then

             if (temp(iz) .lt. 218._rk) then
                so3(iz,iw) = self%yg218(iw)
             endif
             if ((temp(iz) .ge. 218._rk) .and. (temp(iz) .le. 298._rk)) then
                so3(iz,iw) = self%yg218(iw)+(self%yg298(iw)-self%yg218(iw))/(298._rk-218._rk)* &
                     (temp(iz)-218._rk)
             endif
             if (temp(iz) .gt. 298._rk) then
                so3(iz,iw) = self%yg298(iw)
             endif
          endif

       enddo
    enddo


    !------------------------------------------------------
    !     ... QY JPL06
    !         Valid from 306-328 nm
    !                    200-320 K
    !------------------------------------------------------
    do iz = 1, nz
       T = max(min(320.0_rk, temp(iz)),200._rk)

       do iw = 1, nwave

          kt = 0.695_rk * T
          q1 = 1._rk
          q2 = exp(-825.518_rk/kT)

          IF(wc(iw) .LE. 305._rk) THEN
             QY_O1D (iz,iw) = 0.90_rk
          ELSEIF((wc(iw) .GT. 305) .AND. (wc(iw) .LE. 328._rk)) THEN

             QY_O1D(iz,iw)  = 0.0765_rk + &
                  A(1)*                (q1/(q1+q2))*EXP(-((X(1)-wc(iw))/om(1))**4)+ &
                  A(2)*(T/300._rk)**2 *(q2/(q1+q2))*EXP(-((X(2)-wc(iw))/om(2))**2)+ &
                  A(3)*(T/300._rk)**1.5_rk         *EXP(-((X(3)-wc(iw))/om(3))**2)

          ELSEIF(wc(iw) .GT. 328._rk .AND. wc(iw) .LE. 345._rk) THEN
             QY_O1D(iz,iw) = 0.08_rk
          ELSEIF(wc(iw) .GT. 340._rk) THEN
             QY_O1D(iz,iw) = 0._rk
          ENDIF

          QY_O1D(iz,iw) = min(qy_O1D(iz,iw),1.0_rk)
       enddo

    enddo

    xsqy(:) = 0._rk

    chan_loop: do j = 1,self%nchannels
       !------------------------------------------------------
       !     ... derive the cross section*qy
       !------------------------------------------------------

       if (j == 1) then
          qy(:nz,:nwave) = qy_O1D(:nz,:nwave)
       else
          qy(:nz,:nwave) = 1.0_rk - qy_O1D(:nz,:nwave)
       end if

       do iz = 1, nz
          do iw = 1, nwave
             xsqy(iw) = qy(iz,iw)*so3(iz,iw)
          enddo

          xsect(1:nwave) = xsqy(1:nwave)*etf(1:nwave)*esfact
          self%rateconsts(iz,j) = dot_product( fluxes(1:nwave,iz),xsect(1:nwave) )
       end do

    enddo chan_loop

  end subroutine update_rates

end module o3_dissociation_class
