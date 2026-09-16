subroutine refine_grid()
    integer :: ii,jj,hh,iifine,jjfine
    integer :: I_new, J_new

    real, intent(in) :: h
    integer, intent(inout) :: I,J
    integer, intent(inout) :: out
    integer,allocatable,intent(inout) :: material(:,:)

    I_new=I*h
    J_new=J*h

    do ii=1,I_new
        do hh=1,h
            iifine=(ii-1)*h+hh
            dxnew(iifine)=dx(ii)/h
        enddo
    enddo




































