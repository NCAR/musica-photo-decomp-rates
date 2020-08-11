!> \file
!> Tests for the share/numer_mod.F90

!> Test module for the numer_mod functions
program test_share_numer

  implicit none

  call test_inter1( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests for the inter1 function
  subroutine test_inter1( )

    use numer_mod, only : inter1
    use phot_kind_mod, only: rk => kind_phot

    logical :: pass = .true.

    integer, parameter :: n = 10
    integer, parameter :: ng = 5
    integer :: i
    real(rk) :: x(n), y(n)
    real(rk) :: xg(ng), yg(ng)
    real(rk), parameter :: epsilon = 1.e-6_rk

    y = (/ 1._rk, 2._rk, 3._rk, 4._rk, 5._rk, 6._rk, 7._rk, 8._rk, 9._rk, 10._rk /)
    x = (/ 1._rk, 2._rk, 3._rk, 4._rk, 5._rk, 6._rk, 7._rk, 8._rk, 9._rk, 10._rk /)
    xg= (/    1.5_rk,       3.5_rk,       5.5_rk,       7.5_rk,        9.5_rk    /)
    
    ! Here you could set up arrays, parameters, etc to test the inter1( )
    ! function and evaluate the results. If the result look good, do nothing.
    ! If the results look bad abort with an error code

    call inter1( ng, xg, yg, n, x, y )

    check: do i = 1,ng
       if ( yg(i)<xg(i)-epsilon .or. yg(i)>xg(i)+epsilon ) then
          pass=.false.
          exit check
       endif
    end do check

    if( .not. pass ) then
      write(*,*) 'yg: ',yg
      write(*,*) "Something is wrong!"
      stop 3
    end if

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_share_numer
