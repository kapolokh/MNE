subroutine meshmap()
implicit none
!Input
integer, intent(in) :: I,J,material(I,J)
real(8), intent(in) :: h
!Output
integer, allocatable, intent(out) :: cell_id(:,:), cell_row(:), cell_col(:)
!Local


!Isolate triangles from rectangular material input
nt=0
do ii=1,I
    do jj=1,J
        if (material(ii,jj)>0) then
            nt=nt+1
        endif
    enddo
enddo

!Allocate all arrays based on # of triangles in input
allocate(cell_id(I,J),cell_row(I),cell_col(J))


!Assign indexes to triangles: start from bottom left and end at top right
!We map the triangle index onto the row/column x and y.
t=0
do ii=1,I
    do jj=1,J
        if (material(ii,jj)==0) cycle                !skip the nonexistent material marked by 0 in input
        t=t+1
        cell_id(ii,jj)=t                             !Mapping triangle id to physical location
        cell_row(t)=ii
        cell_col(t)=jj
        cell_mat(t)=material(ii,jj)
    enddo
enddo


!Based on input, infer triangle orientation (upward/downward)
do t=1,nt
    ii=cell_row(t)
    jj=cell_col(t)
    if (2*((jj+I-ii)/2)/=jj+I-ii) then
        !Upward triangle if(ii+I-jj) is odd, since integers /2 will produce some rounding it will not equal orig
        !Assign coordinates to the vertices (refer to map below)
        updown(t)=1                  !Triangle
        ix(1,t)=jj-1                 !                    j-1   j   j+1
        iy(1,t)=ii+1                 !            i             3
        ix(2,t)=jj+1                 !                        /   \
        iy(2,t)=ii+1                 !                       /     \
        ix(3,t)=jj                   !                      /       \
        iy(3,t)=ii                   !            i+1      1---------2
    else
        !Downward triangle





































































































end subroutine meshmap