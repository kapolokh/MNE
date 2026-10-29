subroutine Power_iteration(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,Sigma1_2,nuSigmaF, &
                            material,toli,maxinner,tolo,maxouter,phi,keff,BCB,BCR,BCL,BCT)
    implicit none
    integer, intent(in) :: out,I,J,K,M,G,maxinner,maxouter,BCB,BCR,BCL,BCT
    real, intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:)
    real, intent(in) :: SigmaT(:,:),SigmaA(:,:),Sigma1_2(:,:),nuSigmaF(:,:)
    integer, intent(in) :: material(:,:)
    real, intent(in) :: toli,tolo
    real, intent(inout) :: phi(:,:,:)
    real, intent(out) :: keff
    real, allocatable :: SigmaS_self(:,:), q(:,:), Fsrc(:,:), chi(:)
    real :: Fsum_old, Fsum_new, diffk
    integer :: ii,jj,mm,gg,outer

    allocate(q(I,J),Fsrc(I,J),SigmaS_self(M,G),chi(G))

    !Self-scatter only (total out-scatter minus what leaves to the next group)
    do mm=1,M
        do gg=1,G
            SigmaS_self(mm,gg) = SigmaT(mm,gg) - SigmaA(mm,gg) - Sigma1_2(mm,gg)
        end do
    end do

    !Fission spectrum - fast fission only, no thermal fission neutrons born
    chi = 0.0
    chi(1) = 1.0

    phi = 1.0
    keff = 1.0

    do outer = 1,maxouter

        !Fission source from the current flux guess: F(i,j) = sum_g' nuSigmaF(g')*phi(g')
        Fsrc = 0.0
        do jj=1,J
            do ii=1,I
                mm = material(ii,jj)
                do gg=1,G
                    Fsrc(ii,jj) = Fsrc(ii,jj) + nuSigmaF(mm,gg)*phi(ii,jj,gg)
                end do
            end do
        end do
        Fsum_old = sum(Fsrc)

        !Solve groups in sequence (no upscatter, so group g only needs group g-1's
        !already-updated flux, plus the fission source built above)
        do gg=1,G
            do jj=1,J
                do ii=1,I
                    mm = material(ii,jj)
                    q(ii,jj) = chi(gg)*Fsrc(ii,jj)/keff
                    if (gg > 1) then
                        q(ii,jj) = q(ii,jj) + Sigma1_2(mm,gg-1)*phi(ii,jj,gg-1)
                    end if
                end do
            end do
            call inner(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT(:,gg),SigmaS_self(:,gg),material,q, &
           maxinner,toli,phi(:,:,gg),.false.,BCL,BCR,BCB,BCT)
        end do

        !Updated fission source and k
        Fsrc = 0.0
        do jj=1,J
            do ii=1,I
                mm = material(ii,jj)
                do gg=1,G
                    Fsrc(ii,jj) = Fsrc(ii,jj) + nuSigmaF(mm,gg)*phi(ii,jj,gg)
                end do
            end do
        end do
        Fsum_new = sum(Fsrc)

        keff = keff * (Fsum_new/Fsum_old)
        diffk = abs(Fsum_new/Fsum_old - 1.0)

        write(out,'(a,i5,a,f12.6,a,es12.5)') "Outer iteration ", outer, "  k = ", keff, "  diff = ", diffk

        if (diffk <= tolo) then
            write(out,'(/,a)') "Power iteration converged successfully."
            write(out,'(a,i6)') "Outer iterations consumed:", outer
            exit
        end if

        if (outer == maxouter) then
            write(out,'(/,a)') "Power iteration terminated unsuccessfully."
        end if

    end do

end subroutine Power_iteration