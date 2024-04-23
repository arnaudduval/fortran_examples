program conversion
    implicit none


    real, dimension(2, 3) :: a
    real, dimension(2, 3) :: b
    complex, dimension(2, 3) :: c
    real, dimension(2, 3) :: d
    real, dimension(2, 3) :: e

    integer :: i, j

    a = reshape((/ 1., 2., 3., 4., 5., 6. /), shape(a))
    b = reshape((/ 12., 11., 10., 9., 8., 7. /), shape(b))

    c = cmplx(a, b)

    do i = 1, 2
        write(*,*) (c(i, j), j = 1,3)
    enddo
    write(*,*)

    d = realpart(c)
    e = imagpart(c)

    do i = 1, 2
        write(*,*) (d(i, j), j = 1,3)
    enddo
    write(*,*)
    do i = 1, 2
        write(*,*) (e(i, j), j = 1,3)
    enddo





end program conversion