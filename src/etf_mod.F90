module etf_mod
  use phot_kind_mod, only: rk => kind_phot
  use wavelength_grid,only: nwave, wl, wc
  use params_mod, only : hc

  implicit none

  real(rk), allocatable, protected :: photon_flux(:) ! /cm2/sec

contains

  subroutine read_etf(filepath, errmsg, errflg )
    !---------------------------------------------------------------------
    !	... read in ETF
    !---------------------------------------------------------------------
    use netcdf

    !---------------------------------------------------------------------
    !	... dummy args
    !---------------------------------------------------------------------
    character(len=*), intent(in)  :: filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !---------------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------------
    integer :: astat, ret
    integer :: m
    integer :: ncid, dimid, varid
    character(len=64) :: varname

    real(rk) :: etfl(nwave)
    real(rk) :: dw(nwave)

    !---------------------------------------------------------------------
    !	... open file
    !---------------------------------------------------------------------
    ret = nf90_open( trim(filepath), nf90_noclobber, ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_etf: failed to open file ' // trim(filepath)
       return
    end if

    !---------------------------------------------------------------------
    !	... allocate module arrays
    !---------------------------------------------------------------------
    allocate( photon_flux(nwave), stat=astat )
    if( astat /= 0 ) then
       errmsg = 'read_etf: failed to allocate'
       errflg = astat
       return
    endif

    dw(:nwave) = wl(2:nwave+1) - wl(1:nwave)

    !---------------------------------------------------------------------
    !	... read arrays
    !---------------------------------------------------------------------
    ret = nf90_inq_varid( ncid, 'etf', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_etf: failed to get etfl variable id'
       return
    end if

    if (varid>0) then
       ret = nf90_get_var( ncid, varid, etfl )
       if( ret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'read_etf: failed to read etfl variable'
          return
       end if
       photon_flux(:nwave) = dw(:nwave)*etfl(:nwave)*1.e-13_rk*wc(:nwave)/hc !  watts/m2/nm --> /cm2/sec
    end if

    !---------------------------------------------------------------------
    !	... close file
    !---------------------------------------------------------------------
    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'read_etf: failed to close file ' // trim(filepath)
       return
    end if

  end subroutine read_etf

end module etf_mod
