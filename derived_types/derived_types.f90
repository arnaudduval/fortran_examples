module types
    type color
        character(len=25) :: name
        integer, dimension(3) :: code_rvb
    end type color

    type pen
        character(len=25) :: name
        integer :: thickness
    end type pen
end module types



program test
    use types, only: color, pen

    implicit none

    type(color), dimension(3) :: tabrvb
    type(pen) :: crayon

    tabrvb =    &
        (/ color("Red", (/ 255,0,0 /)),                 &
           color("Green", (/ 0,255,0 /)),               &
           color("Blue", (/ 0,0,255 /)) /)

    crayon = pen("HB5", 12)


    write(*,*) tabrvb(1)
    write(*,*) crayon

    call printcolor(tabrvb(1))
    call print_color_and_pen(tabrvb(1), crayon)

end program test


subroutine printcolor(col)
    use types, only: color

    implicit none

    type(color) :: col

    write(*,*) col%name, col%code_rvb
end subroutine printcolor

subroutine print_color_and_pen(col, stylo)
    use types, only: color, pen

    implicit none

    type(color) :: col
    type(pen) :: stylo

    write(*,*) col%name, col%code_rvb
    write(*,*) stylo%name, stylo%thickness
end subroutine print_color_and_pen