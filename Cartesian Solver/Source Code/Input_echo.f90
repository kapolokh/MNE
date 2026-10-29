subroutine Input_echo(out,I,J,K,G,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,SigmaA,nuSigmaF,Sigma1_2,material)
        implicit none
        integer, intent(in) :: out,I,J,K,G
        integer, intent(in) :: BCL,BCR,BCB,BCT
        real, intent(in) :: mu(:),eta(:),w(:),dx(:),dy(:),SigmaT(:,:),SigmaA(:,:),nuSigmaF(:,:),Sigma1_2(:,:)
        integer,intent(in) :: material(:,:)
        integer :: n,ii,jj,m,gg
        !We print out some of the necessary input variables to output file
        !Write the Discrete Ordinates Line
        write(out,'(/,a)') "Discrete Ordinates/Octant:"
        write(out,'(a3,1x,a10,1x,a10,1x,a10)') "n", "mu", "eta", "w"
        do n=1,K 
            write(out,'(i3,3(1x,f10.5))') n, mu(n), eta(n), w(n)
        end do

        !Write the BC Line
        write(out,'(/,a)') "Boundary Conditions:"
        write(out,'(a6,1x,a6,2x,a6,1x,a6)')  "Left", "Right", "Bottom", "Top"
        write(out,'(i6,3(1x,i6))') BCL,BCR,BCB,BCT

        !Write the cell data line
        write(out,'(/,a)') "Computational Cell Data:"
        write(out,'(a3,1x,a3,1x,a8,1x,a5,1x,a10,1x,a10,1x,a10,1x,a10,1x,a10,1x,a10,1x,a10)') &
        "i", "j", "Material", "Group", "dx", "dy", "SigmaT", "SigmaA", "nuSigmaF", "Sigma1_2"
        do jj = 1,J
            do ii = 1,I
                m = material(ii,jj)
                    do gg = 1,G
                        write(out,'(i3,1x,i3,1x,i8,1x,i5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5)') &
                        ii,jj,m,gg,dx(ii),dy(jj),SigmaT(m,gg),SigmaA(m,gg),nuSigmaF(m,gg),Sigma1_2(m,gg)
                    end do
            end do
        end do

    end subroutine Input_echo