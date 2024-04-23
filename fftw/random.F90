!! TODO factorize those functions in a separate file

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

subroutine dft(in_r, in_i, out_r, out_i, nx, ny)
    !! Compute slow Fourier transform
    implicit none

    integer, intent(in) :: nx, ny
    real, dimension(nx, ny), intent(in) :: in_r, in_i
    real, dimension(nx, ny), intent(out) :: out_r, out_i

    integer :: k, l, m, n
    complex :: temp
    real :: pi = 4.0 * atan(1.0)
    complex :: arg

    out_r = 0.0
    out_i = 0.0

    do k =0, nx-1
        do l=0, ny-1

            do n=0, nx-1
                do m=0, ny-1
                    arg = -2.*pi*cmplx(0.0, real(n)*real(k)/real(nx) + real(m)*real(l)/real(ny))
                    ! write(*,*) k, l, n, m, cmplx(in_r(n+1, m+1), in_i(n+1, m+1)), cexp(arg)
                    temp = cmplx(in_r(n+1, m+1), in_i(n+1, m+1))*cexp(arg)
                    out_r(k+1, l+1) = out_r(k+1, l+1) + real(temp)
                    out_i(k+1, l+1) = out_i(k+1, l+1) + aimag(temp)
                enddo
            enddo

        enddo
    enddo


end subroutine dft



program random_fft
    use iso_c_binding
    implicit none

    include 'fftw3.f03'

    integer nx, ny

    real, allocatable, dimension(:, :) :: in_i, in_r, out_i, out_r, out_id, out_rd
    complex, allocatable, dimension(:, :) :: in, out
    integer :: i, j

    type(c_ptr) :: plan

    nx = 4
    ny = 7

    allocate(in_r(nx, ny), in_i(nx, ny), out_r(nx, ny), out_i(nx, ny), out_rd(nx, ny), out_id(nx, ny))
    allocate(in(nx, ny), out(nx, ny))
    call random_init(.true., .true.)

    call random_number(in_r)
    call random_number(in_i)
    !!in_i(:,:) = 0.0

    do i = 1, nx
        do j = 1, ny
            in(i, j) = cmplx(in_r(i, j), in_i(i, j))
        enddo
    enddo


    do i = 1, nx
        write(*,*) in_r(i, :)
    enddo
    write(*,*)
    do i = 1, nx
        write(*,*) in_i(i, :)
    enddo
    write(*,*)
    do i = 1, nx
        write(*,*) in(i, :)
    enddo

    !! FFT with FFTW
    plan = fftwf_plan_dft_2d(ny, nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    call fftwf_execute_dft(plan, in, out)
    call fftwf_destroy_plan(plan)

    !! FFT with legacy code
    call legacy_fft2d_forward(in_r, in_i, out_r, out_i, nx, ny)

    !! DFT slow
    call dft(in_r, in_i, out_rd, out_id, nx, ny)

    write(*,*)
    write(*,*) "Error Legacy/FFTW:", norm2(out_rd-real(out))/norm2(out_rd), &
    norm2(out_id-aimag(out))/norm2(out_id)
    write(*,*)

    write(*, *) "FFT ISAAC"
    do i = 1, nx
        write(*,*) out_r(i, :)
    enddo
    write(*,*)
    do i = 1, nx
        write(*,*) out_i(i, :)
    enddo
    write(*,*)
    write(*,*) "DFT"
    do i = 1, nx
        write(*,*) out_rd(i, :)
    enddo
    write(*,*)
    do i = 1, nx
        write(*,*) out_id(i, :)
    enddo
    write(*,*)
    write(*,*) "FFTW"
    do i = 1, nx
        write(*,*) out(i, :)
    enddo

end program random_fft
