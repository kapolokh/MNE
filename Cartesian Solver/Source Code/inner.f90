subroutine inner(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,SigmaS,material,source,maxinner,toli,phi,verbose, &
                  BCL,BCR,BCB,BCT)
    implicit none
    !Input:
    !I,J - number of spatial cells
    !K - # of angles per octant
    !M - # of materials
    !dx,dy - width of the cell
    !mu,eta - ordinates
    !w - weight
    !SigmaT - total X/S
    !SigmaS - self-scatter X/S for this group
    !material - material map
    !source - fixed source in each material cell
    !toli - convergence criterion
    !BCL,BCR,BCB,BCT - boundary flags, 0 = vacuum, 1 = reflective

    !Output:
    !phi - scalar flux after inner iterations

    integer,intent(in) :: I,J,K,M,out
    real, intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:)
    real, intent(in) :: SigmaT(:),SigmaS(:),source(:,:)
    integer, intent(in) :: material(:,:)
    integer, intent(in) :: maxinner
    real, intent(in) :: toli
    logical, intent(in) :: verbose
    integer, intent(in) :: BCL,BCR,BCB,BCT
    real, intent(out) :: phi(:,:)
    real,allocatable :: phi_old(:,:), q(:,:)
    real,allocatable :: leftflux(:,:,:), rightflux(:,:,:), bottomflux(:,:,:), topflux(:,:,:)
    real :: diff,maxdiff
    integer :: ii,jj,mm,it

    allocate (phi_old(I,J) , q(I,J))
    allocate (leftflux(4,K,J), rightflux(4,K,J), bottomflux(4,K,I), topflux(4,K,I))
    leftflux = 0.0
    rightflux = 0.0
    bottomflux = 0.0
    topflux = 0.0

    phi_old = 0.0d0
    phi = 0.0d0

    do it = 1,maxinner
        do jj = 1,J
            do ii = 1,I
                mm = material(ii,jj)
                q(ii,jj) = source(ii,jj) + SigmaS(mm)*phi_old(ii,jj)
            end do
        end do

        call sweep(I,J,K,M,dx,dy,mu,eta,w,SigmaT,material,q,phi,BCL,BCR,BCB,BCT, &
                   leftflux,rightflux,bottomflux,topflux)

        maxdiff = 0.0d0

        do jj = 1,J
            do ii = 1,I
                if (phi_old(ii,jj) /= 0.0d0) then
                    diff = abs(phi(ii,jj)/phi_old(ii,jj)-1.0d0)
                else
                    diff = abs(phi(ii,jj) - phi_old(ii,jj))
                end if
                if (diff > maxdiff) then
                    maxdiff = diff
                end if
            end do
        end do

        if (maxdiff <= toli) then
            if (verbose) then
                write(out,'(/,a)') "Inner iterations converged successfully."
                write(out,'(a,i6)') "Iterations consumed:", it
                write(out,'(a,1x,es12.5)') "Convergence criterion achieved:", maxdiff
            end if
            exit
        end if

        if (it == maxinner) then
            if (verbose) then
                write(out,'(/,a)') "Inner iterations terminated unsuccessfully."
                write(out,'(a,i6)') "Iterations consumed:", it
                write(out,'(a,1x,es12.5)') "Convergence criterion achieved:", maxdiff
            end if
            exit
        end if

        phi_old = phi
    end do

    if (verbose) then
        write(out,'(/,a)') "Discrete Ordinates Method Solution"
        write(out,'(a3,1x,a3,1x,a26)') "i","j","Cell-Averaged Scalar Flux"
        do jj = 1,J
            do ii = 1,I
                write(out,'(i3,1x,i3,1x,es16.8)') ii,jj,phi(ii,jj)
            end do
        end do
    end if

end subroutine inner