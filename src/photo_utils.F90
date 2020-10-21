module photo_utils

  use phot_kind_mod, only: rk => kind_phot
  use numer_mod, only: addpnt, inter2, inter4

  implicit none

contains

  SUBROUTINE add_pnts_inter2(xin,yin,yout,kdata,n,nw,wl,jlabel,deltax,yends, errmsg, errflg)

    integer, intent(in) :: kdata
    integer, intent(in) :: n
    integer, intent(in) :: nw
    real(rk), intent(in)    :: deltax
    real(rk), intent(in)    :: wl(nw)
    real(rk), intent(in)    :: xin(kdata)
    real(rk), intent(in)    :: yin(kdata)
    real(rk), intent(in)    :: yends(2)
    real(rk), intent(inout) :: yout(nw-1)
    character(len=*), intent(in) :: jlabel
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m
    real(rk)    :: xwrk(kdata), ywrk(kdata)

    errmsg = ' '
    errflg = 0

    m = n 
    xwrk(1:n) = xin(1:n)
    ywrk(1:n) = yin(1:n)
    CALL addpnt(xwrk,ywrk,kdata,m,xin(1)*(1._rk-deltax),yends(1),errmsg, errflg)
    CALL addpnt(xwrk,ywrk,kdata,m,              0._rk,yends(1),errmsg, errflg)
    CALL addpnt(xwrk,ywrk,kdata,m,xin(n)*(1._rk+deltax),yends(2),errmsg, errflg)
    CALL addpnt(xwrk,ywrk,kdata,m,          1.e+38_rk,yends(2),errmsg, errflg)

    CALL inter2(nw,wl,yout,m,xwrk,ywrk,errmsg, errflg)

  END SUBROUTINE add_pnts_inter2

  SUBROUTINE base_read( filespec, errmsg, errflg, skip_cnt, rd_cnt,x, y, y1, y2, y3, y4, y5 )

    integer, optional, intent(in) :: skip_cnt
    integer, intent(inout)        :: rd_cnt
    real(rk), intent(inout)           :: x(:), y(:)
    real(rk), optional, intent(inout) :: y1(:), y2(:), y3(:)
    real(rk), optional, intent(inout) :: y4(:), y5(:)
    character(len=*), intent(in)  :: filespec
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i, idum, unitn
    integer :: y_to_rd
    integer :: ios

    errmsg = ' '
    errflg = 0

    y_to_rd = 1
    if( present(y5) ) y_to_rd = y_to_rd + 1
    if( present(y4) ) y_to_rd = y_to_rd + 1
    if( present(y3) ) y_to_rd = y_to_rd + 1
    if( present(y2) ) y_to_rd = y_to_rd + 1
    if( present(y1) ) y_to_rd = y_to_rd + 1

    OPEN(newunit=unitn,FILE=trim(filespec),STATUS='old',IOSTAT=ios)
    IF( ios /= 0 ) then
       write(errmsg,'(''base_read: failed to open '',a)') trim(filespec)
       errflg = ios
       return
    ENDIF

    if( present(skip_cnt) ) then
       DO i = 1, skip_cnt
          READ(unitn,*,IOSTAT=ios)
          IF( ios /= 0 ) exit
       END DO
    else
       READ(unitn,*,IOSTAT=ios) idum,rd_cnt
       IF( ios == 0 ) then
          DO i = 1, idum-2
             READ(unitn,*,IOSTAT=ios)
             IF( ios /= 0 ) exit
          ENDDO
       ENDIF
    endif

    IF( ios /= 0 ) then
       write(errmsg,'(''base_read: failed to read '',a)') trim(filespec)
       errflg = ios
       return
    ENDIF

    select case( y_to_rd )
    case( 1 )
       DO i = 1, rd_cnt
          READ(unitn,*,IOSTAT=ios) x(i), y(i)
          IF( ios /= 0 ) exit
       END DO
    case( 2 )
       DO i = 1, rd_cnt
          READ(unitn,*,IOSTAT=ios) x(i), y(i), y1(i)
          IF( ios /= 0 ) exit
       END DO
    case( 3 )
       DO i = 1, rd_cnt
          READ(unitn,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i)
          IF( ios /= 0 ) exit
       END DO
    case( 4 )
       DO i = 1, rd_cnt
          READ(unitn,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i), y3(i)
          IF( ios /= 0 ) exit
       END DO
    case( 5 )
       DO i = 1, rd_cnt
          READ(unitn,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i), y3(i),y4(i)
          IF( ios /= 0 ) exit
       END DO
    case( 6 )
       DO i = 1, rd_cnt
          READ(unitn,*,IOSTAT=ios) x(i), y(i), y1(i), y2(i), y3(i),y4(i),y5(i)
          IF( ios /= 0 ) exit
       END DO
    end select

    CLOSE (unitn)

    IF( ios /= 0 ) then
       write(errmsg,'(''base_read: failed to read '',a)') trim(filespec)
       errflg = ios
       return
    ENDIF

  END SUBROUTINE base_read

end module photo_utils
