program test
    real :: start, finish
    integer :: i, j
    call cpu_time(start)
        ! put code to test here
        do i=1,100000000
            j = j + 0.
        end do
    call cpu_time(finish)
    print '("Elapsed Time = ",f6.3," seconds.")',finish-start
end program test