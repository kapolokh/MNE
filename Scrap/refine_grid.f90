subroutine refine_grid(I,J,h,out,material,dx,dy)
    implicit none

    integer :: ii,jj,hh,iifine,jjfine,hi,hj
    integer :: I_new, J_new
    real,allocatable :: dxnew(:),dynew(:)
    real,allocatable :: mat_new(:,:)
   
    integer, intent(inout) :: I,J
    integer, intent(in) :: h
    integer, intent(inout) :: out
    integer,allocatable,intent(inout) :: material(:,:)
    real, allocatable,intent(inout) :: dx(:),dy(:)
 
  
    I_new=I*h
    J_new=J*h
    allocate(dxnew(I_new))
    allocate(dynew(J_new))
    allocate(mat_new(I_new,J_new))


    do ii=1,I
        do hh=1,h
            iifine=(ii-1)*h+hh
            dxnew(iifine)=dx(ii)/h
        enddo
    enddo

    do jj=1,J
        do hh=1,h
            jjfine=(jj-1)*h+hh
            dynew(jjfine)=dy(jj)/h
        enddo
    enddo

    do ii=1,I
        do jj=1,J
            do hi=1,h
                do hj=1,h
                    iifine=(ii-1)*h+hi
                    jjfine=(jj-1)*h+hj
                    mat_new(iifine,jjfine)=material(ii,jj)
                enddo
            enddo
        enddo
    enddo

    I=I_new
    J=J_new
    dx=dxnew
    dy=dynew
    material=mat_new

end subroutine refine_grid 



































