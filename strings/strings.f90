subroutine string_as_argument(mystring)
    implicit none

    character(len=*), intent(in) :: mystring

    write(*,*) mystring


end subroutine string_as_argument


program test
    implicit none

    character(len=256) :: mystring

    mystring = 'Hello world!'

    call string_as_argument(mystring)

end program test