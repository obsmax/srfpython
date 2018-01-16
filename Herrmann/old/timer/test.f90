program test


    REAL*8 del1
    del1 = -2.0
    write (*,*) 1.0D+00, DSIGN(1.0D+00,del1)


end program test


!    real :: start, finish
!    integer :: i, j
!    call cpu_time(start)
!        ! put code to test here
!        do i=1,100000000
!            j = j + 0.
!        end do
!    call cpu_time(finish)
!    print '("Elapsed Time = ",f6.3," seconds.")',finish-start