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

    use numer_mod,                     only : inter1

    logical :: results = .true.

    ! Here you could set up arrays, parameters, etc to test the inter1( )
    ! function and evaluate the results. If the result look good, do nothing.
    ! If the results look bad abort with an error code

    ! call inter1( ng, nx, yg, n, x, y )

    if( .not. results ) then
      write(*,*) "Something's wrong!"
      stop 3
    end if

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_share_numer
