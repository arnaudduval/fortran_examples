real function CN3(x, y, dx, dy, E, nu)
    !! Compute Love coefficient Fz <-> uz for a general point x, y
    implicit none

    real, intent(in) :: x, y, dx, dy, E, nu

    !! Local variables
    double precision, parameter :: PI=4.D0*datan(1.D0)
    real :: xa1, xa2, yb1, yb2
    real :: cN31, cN32, cN33, cN34, cN35, cN36, cN37, cN38

    xa1 = x + dx/2
    xa2 = x - dx/2.
    yb1 = y + dy/2.
    yb2 = y - dy/2.

    cN31=xa1*log(-yb1+sqrt(yb1*yb1+xa1*xa1))
    cN32=yb1*log(-xa1+sqrt(yb1*yb1+xa1*xa1))
    cN33=xa2*log(-yb1+sqrt(yb1*yb1+xa2*xa2))
    cN34=yb2*log(-xa1+sqrt(yb2*yb2+xa1*xa1))
    cN35=xa1*log(-yb2+sqrt(yb2*yb2+xa1*xa1))
    cN36=yb1*log(-xa2+sqrt(yb1*yb1+xa2*xa2))
    cN37=xa2*log(-yb2+sqrt(yb2*yb2+xa2*xa2))
    cN38=yb2*log(-xa2+sqrt(yb2*yb2+xa2*xa2))

    cN3 = -1.*((1. - nu**2)/(real(PI)*E))*(cN31+cN32-cN33-cN34-cN35-cN36+cN37+cN38)
end function CN3

subroutine wrap_around(g, g_wrap, nx, ny)
    !! perform symmetric wraparound of coefficients
    implicit none

    integer, intent(in) :: nx, ny
    real, dimension(nx, ny), intent(in) :: g
    real, dimension(2*nx, 2*ny), intent(out) :: g_wrap

    integer :: i, j

    g_wrap(:,:) = 0.

    do i=1, nx
        do j=1, ny
            g_wrap(i, j) = g(i, j)
            if (j+ny+1.le.2*ny) then
                !Second section, extension along y
                g_wrap(i, j+ny+1) = g(i, ny+1-j)
            endif
            if (i+nx+1.le.2*nx) then
                !Third section, extension along x
                g_wrap(i+nx+1, j) = g(nx+1-i, j)
            endif
            if (i+nx+1.le.2*nx.and.j+ny+1.le.2*ny) then
                !Fourth section, extension along x and y
                g_wrap(i+nx+1, j+ny+1) = g(nx+1-i, ny+1-j)
            endif
        enddo
    enddo


end subroutine wrap_around


program test_fftw
    use iso_c_binding
    implicit none

    include 'fftw3.f03'

    integer :: nx, ny
    real :: dx, dy
    real :: nu, E, Estar
    real, allocatable, dimension(:, :) :: g_love, g_love_wrap
    real :: cN3

    integer :: i, j

    !! Specific for fftw
    complex, allocatable, dimension(:, :) :: in, out
    type(c_ptr) :: plan


    !! Initialization
    E = 210000.
    nu = 0.3
    Estar = E/(1. - nu**2.)
    dx = 0.01
    dy = 0.01
    nx = 49
    ny = 49


    allocate(g_love(nx, ny), g_love_wrap(2*nx, 2*ny))
    allocate(in(2*nx, 2*ny), out(2*nx, 2*ny))

    do i=1, nx
        do j=1, ny
            g_love(i, j) = CN3((i-1)*dx, (j-1)*dy, dx, dy, E, nu)
        enddo
    enddo

    call wrap_around(g_love, g_love_wrap, nx, ny)

    !! Implicit conversion to complex type
    in(:,:) = g_love_wrap(:,:)

    write(*,*) minval(real(in)), maxval(real(in))
    write(*,*) minval(aimag(in)), maxval(aimag(in))


    !! WARNING: attention should be care to storage convention (row major/col major)
    !! It is not a problem for this purely symmetric case
    plan = fftwf_plan_dft_2d(2*ny, 2*nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    call fftwf_execute_dft(plan, in, out)
    call fftwf_destroy_plan(plan)



    write(*,*) minval(real(out)), maxval(real(out))
    write(*,*) minval(aimag(out)), maxval(aimag(out))


end program test_fftw