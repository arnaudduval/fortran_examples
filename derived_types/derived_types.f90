module types
    type nurbspatch
        integer :: dimension
        integer, dimension(:), allocatable :: Jpqr         !! degrees
        double precision, dimension(:), allocatable :: Ukv1, Ukv2, Ukv3

        contains

        procedure :: initPatch

    end type nurbspatch

    contains
        subroutine initPatch(self, dim, Jpqr, Ukv, kv_sizes)
            implicit none

            class(nurbspatch) :: self
            integer :: dim
            integer, dimension(dim) :: Jpqr
            integer, dimension(dim) :: kv_sizes
            double precision, dimension(:) :: Ukv

            self%dimension = dim

            if (allocated(self%Jpqr)) deallocate(self%Jpqr)
            allocate(self%Jpqr(dim))
            self%Jpqr(:) = Jpqr(:)

            if (allocated(self%Ukv1)) deallocate(self%Ukv1)
            allocate(self%Ukv1(kv_sizes(1)))
            self%Ukv1(:) = Ukv(:kv_sizes(1))

            if (dim > 1) then
                if (allocated(self%Ukv2)) deallocate(self%Ukv2)
                allocate(self%Ukv2(kv_sizes(2)))
                self%Ukv2(:) = Ukv(kv_sizes(1)+1:kv_sizes(1)+kv_sizes(2))
            endif

            if (dim > 2) then
                if (allocated(self%Ukv3)) deallocate(self%Ukv3)
                allocate(self%Ukv3(kv_sizes(3)))
                self%Ukv3(:) = Ukv(kv_sizes(1)+kv_sizes(2)+1:   &
                                   kv_sizes(1)+kv_sizes(2)+kv_sizes(3))
            endif

        end subroutine initPatch
end module types



program test
    use types, only: nurbspatch

    implicit none

    type(nurbspatch) :: patch

    call patch%initPatch(3, (/ 2, 2, 2 /),       &
                    &   (/ 0.0D0, 0.0D0, 0.0D0, 0.25D0, 0.5D0, 0.75D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    &      0.5D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0 /), &
                    &   (/ 9, 7, 6 /) )

    write(*,'(10F5.2)') patch%Ukv1
    write(*,'(10F5.2)') patch%Ukv2
    write(*,'(10F5.2)') patch%Ukv3

    call manip(patch)

end program test


subroutine manip(patch)
    use types, only: nurbspatch

    implicit none

    type(nurbspatch) :: patch

    write(*,*) 'In subroutine'
    write(*,'(10F5.2)') patch%Ukv1
    write(*,'(10F5.2)') patch%Ukv2
    write(*,'(10F5.2)') patch%Ukv3


end subroutine manip