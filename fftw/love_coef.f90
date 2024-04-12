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

real function Hertz_pressure_distribution(r, Ph, a)
    !! Compute Hertz pressure distribution at a given radius
    implicit none

    real, intent(in) :: r   !! Radius
    real, intent(in) :: Ph  !! Hertz pressure
    real, intent(in) :: a   !! Hertz radius

    if (r >= a) then
        Hertz_pressure_distribution = 0.
    else
        Hertz_pressure_distribution = Ph*(1. - (r**2.) / (a**2.))**0.5
    endif
end function Hertz_pressure_distribution

real function Hertz_pressure_distribution_xy(x, y, Ph, a)
    !! Compute Hertz pressure distribution at given coordinates x/y
    implicit none

    real, intent(in) :: x, y
    real, intent(in) :: Ph
    real, intent(in) :: a

    real :: Hertz_pressure_distribution

    Hertz_pressure_distribution_xy = Hertz_pressure_distribution((x**2. + y**2.)**0.5, Ph, a)
end function Hertz_pressure_distribution_xy

real function Hertz_pressure(load, a)
    !! Compute Hertz pressure value for a shpere/plane contact
    implicit none

    real, parameter :: pi = 3.141592653
    real, intent(in) :: load
    real, intent(in) :: a   !! Hertz radius

    Hertz_pressure = 3.*load/(2.*pi*a**2.)
end function Hertz_pressure

real function Hertz_radius(load, radius, E, nu)
    !! Compute Hertz radius fo a shpare plane contact
    implicit none

    real, intent(in) :: load
    real, intent(in) :: radius  !! Sphere radius
    real, intent(in) :: E, nu  !! Linear elasticity parameters

    real :: Estar

    Estar = 0.5*E/(1. - nu**2.)
    Hertz_radius = (3.*load*radius/(4.*Estar))**(1./3.)
end function Hertz_radius



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

subroutine legacy_fft2d_forward(input_r, input_i, output_r, output_i, n1, n2)
    !! Compute forward FFT2D with legacy source code
    implicit none

    integer, intent(in) :: n1, n2
    real, dimension(n1, n2), intent(in) :: input_r, input_i
    real, dimension(n1, n2), intent(out) :: output_r, output_i

    integer :: n3

    !! legacy routine fft computes in place: a copy must be made
    output_r(:, :) = input_r(:, :)
    output_i(:, :) = input_i(:, :)

    n3 = n1*n2
    call fft(output_r, output_i, n3, n1, n1, 1)
    call fft(output_r, output_i, n3, n2, n3, 1)
end subroutine legacy_fft2d_forward

subroutine legacy_fft2d_backward(input_r, input_i, output_r, output_i, n1, n2)
    !! Compute forward FFT2D with legacy source code
    implicit none

    integer, intent(in) :: n1, n2
    real, dimension(n1, n2), intent(in) :: input_r, input_i
    real, dimension(n1, n2), intent(out) :: output_r, output_i

    integer :: n3

    !! legacy routine fft computes in place: a copy must be made
    output_r(:, :) = input_r(:, :)
    output_i(:, :) = input_i(:, :)

    n3 = n1*n2
    call fft(output_r, output_i, n3, n1, n1, -1)
    call fft(output_r, output_i, n3, n2, n3, -1)
end subroutine legacy_fft2d_backward



program test_fftw
    use iso_c_binding
    implicit none

    include 'fftw3.f03'

    integer :: nx, ny
    real :: dx, dy
    real :: nu, E
    real :: load, radius
    real, allocatable, dimension(:, :) :: g_love, g_love_wrap
    complex, allocatable, dimension(:,:) :: gf_love_wrap
    real, allocatable, dimension(:, :) :: pres, pres_zp
    complex, allocatable, dimension(:,:) :: presf_zp
    complex, allocatable, dimension(:, :) :: uz3, uz3f

    real :: cN3

    integer :: i, j
    real :: Hertz_radius, Hertz_pressure, Hertz_pressure_distribution_xy
    real :: a, Ph

    !! Specific for fftw
    complex, allocatable, dimension(:, :) :: in, out
    type(c_ptr) :: plan

    !! Specific for lecacy fft
    real, allocatable, dimension(:, :) :: out_leg_r, out_leg_i, zero


    !! Initialization
    E = 210000.
    nu = 0.3
    dx = 0.01
    dy = 0.01
    nx = 11
    ny = 11
    load = 1.
    radius = 1000.


    allocate(g_love(nx, ny), g_love_wrap(2*nx, 2*ny), gf_love_wrap(2*nx, 2*ny))
    allocate(pres(nx, ny), pres_zp(2*nx, 2*ny), presf_zp(2*nx, 2*ny))
    allocate(uz3(2*nx, 2*ny), uz3f(2*nx, 2*ny))
    allocate(in(2*nx, 2*ny), out(2*nx, 2*ny))
    allocate(out_leg_r(2*nx, 2*ny), out_leg_i(2*nx, 2*ny), zero(2*nx, 2*ny))
    zero(:, :) = 0.0

    !! Compute influence coefficients
    do i=1, nx
        do j=1, ny
            g_love(i, j) = CN3((i-1)*dx, (j-1)*dy, dx, dy, E, nu)
        enddo
    enddo

    call wrap_around(g_love, g_love_wrap, nx, ny)

    !! Implicit conversion to complex type
    in(:,:) = g_love_wrap(:,:)

    write(*,*) "Influence coefficients in physical space"
    write(*,*) "----------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(in)), maxval(real(in))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(in)), maxval(aimag(in))
    write(*,*)


    !! WARNING: attention should be care to storage convention (row major/col major)
    !! It is not a problem for this purely symmetric case
    plan = fftwf_plan_dft_2d(2*ny, 2*nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    call fftwf_execute_dft(plan, in, out)
    call fftwf_destroy_plan(plan)
    gf_love_wrap(:,:) = out(:,:)

    write(*,*) "Influence coefficients in Fourier space"
    write(*,*) "---------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(gf_love_wrap)), maxval(real(gf_love_wrap))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(gf_love_wrap)), maxval(aimag(gf_love_wrap))
    write(*,*)

    write(*,*) "Influence coefficients in physical space (are they changed ?)"
    write(*,*) "-------------------------------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(in)), maxval(real(in))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(in)), maxval(aimag(in))
    write(*,*)

    call legacy_fft2d_forward(g_love_wrap, zero, out_leg_r, out_leg_i, 2*nx, 2*ny)

    write(*,*) "Legagcy influence coefficients in Fourier space"
    write(*,*) "-----------------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(gf_love_wrap)), maxval(real(gf_love_wrap))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(gf_love_wrap)), maxval(aimag(gf_love_wrap))
    write(*,*)

    !! Compute input pressure
    a = Hertz_radius(load, radius, E, nu)
    Ph = Hertz_pressure(load, a)

    do i = 1, nx
        do j = 1, ny
            pres(i, j) = Hertz_pressure_distribution_xy(dx*(i - nx/2 - 1), dy*(i - ny/2 - 1), Ph, a)
        enddo
    enddo

    pres_zp(:, :) = 0.
    pres_zp(:nx, :ny) = pres(:,:)

    in(:,:) = pres_zp(:,:)
    plan = fftwf_plan_dft_2d(2*ny, 2*nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    call fftwf_execute_dft(plan, in, out)
    call fftwf_destroy_plan(plan)
    presf_zp(:,:) = out(:,:)
    write(*,*) "Pressure in Fourier space"
    write(*,*) "---------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(presf_zp)), maxval(real(presf_zp))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(presf_zp)), maxval(aimag(presf_zp))
    write(*,*)

    !! Compute normal displacement
    uz3f = presf_zp * gf_love_wrap
    write(*,*) "Normal displacement in Fourier space"
    write(*,*) "------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(uz3f)), maxval(real(uz3f))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(uz3f)), maxval(aimag(uz3f))
    write(*,*)

    plan = fftwf_plan_dft_2d(2*ny, 2*nx, uz3f, uz3, FFTW_BACKWARD, FFTW_ESTIMATE)
    call fftwf_execute_dft(plan, uz3f, uz3)
    call fftwf_destroy_plan(plan)

    write(*,*) "Normal displacement in physical space"
    write(*,*) "-------------------------------------"
    write(*,*) "Real part, min/max: ", minval(real(uz3/(4*nx*ny))), maxval(real(uz3/(4*nx*ny)))
    write(*,*) "Imaginary part, min/max: ", minval(aimag(uz3/(4*nx*ny))), maxval(aimag(uz3/(4*nx*ny)))
    write(*,*)


end program test_fftw
