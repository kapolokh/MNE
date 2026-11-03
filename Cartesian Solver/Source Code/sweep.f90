subroutine sweep(I,J,K,M,dx,dy,mu,eta,w,SigmaT,material,source,phi,BCL,BCR,BCB,BCT, &
                  leftflux,rightflux,bottomflux,topflux)
    implicit none
    integer, intent(in) :: I,J,K,M
    real, intent(in) :: dx(:), dy(:), mu(:), eta(:), w(:), SigmaT(:), source(:,:)
    integer, intent(in) :: material(:,:)
    integer, intent(in) :: BCL,BCR,BCB,BCT
    real, intent(out) :: phi(:,:)
    real, intent(inout) :: leftflux(:,:,:), rightflux(:,:,:), bottomflux(:,:,:), topflux(:,:,:)
    !leftflux/rightflux dims: (quad,n,jj)   bottomflux/topflux dims: (quad,n,ii)
    real,allocatable :: xbound(:,:), ybound(:,:)
    real :: mun, etan
    real :: psi
    real :: psi_x_in
    real :: psi_y_in
    real :: psi_x_out
    real :: psi_y_out
    integer :: mat
    integer :: quad, qx, qy
    integer :: j_start, j_end, j_step,i_start, i_end, i_step
    integer :: n,ii,jj

    allocate(xbound(0:I,J),ybound(I,0:J))
    phi = 0.0

    !Adjust ordinate signs depending on direction of the sweep
    !Do sweep in 4 directions to get all particles travelling in all directions in mesh.
    !qx/qy identify which quadrant reflects into this one across the x/y boundary
    do quad = 1,4
        do n = 1,K
            select case(quad)
            case(1)
                mun = mu(n)
                etan = eta(n)
                qx = 2
                qy = 3
            case(2)
                mun = -mu(n)
                etan = eta(n)
                qx = 1
                qy = 4
            case(3)
                mun = mu(n)
                etan = -eta(n)
                qx = 4
                qy = 1
            case(4)
                mun = -mu(n)
                etan = -eta(n)
                qx = 3
                qy = 2
            end select

            !Define starting and ending positions for the sweep
            if (mun > 0.0) then
                i_start = 1
                i_end = I
                i_step = 1
            else
                i_start = I
                i_end = 1
                i_step = -1
            end if

            if (etan > 0.0) then
                j_start = 1
                j_end = J
                j_step = 1
            else
                j_start = J
                j_end = 1
                j_step = -1
            end if

            xbound = 0.0
            ybound = 0.0

            !Reflective inflow at the domain edge: use the reflected quadrant's
            !exiting flux, saved from its most recent sweep, only where the
            !relevant BC flag is 1. Otherwise stays 0.0 (vacuum).
            if (mun > 0.0 .and. BCL == 1) then
                do jj=1,J
                    xbound(0,jj) = leftflux(qx,n,jj)
                end do
            else if (mun < 0.0 .and. BCR == 1) then
                do jj=1,J
                    xbound(I,jj) = rightflux(qx,n,jj)
                end do
            end if

            if (etan > 0.0 .and. BCB == 1) then
                do ii=1,I
                    ybound(ii,0) = bottomflux(qy,n,ii)
                end do
            else if (etan < 0.0 .and. BCT == 1) then
                do ii=1,I
                    ybound(ii,J) = topflux(qy,n,ii)
                end do
            end if

            do jj = j_start,j_end,j_step
                do ii = i_start,i_end, i_step
                    !Define incoming fluxes into cells
                    if (mun > 0.0) then
                        psi_x_in = xbound(ii-1,jj)
                    else
                        psi_x_in = xbound(ii,jj)
                    end if
                    if (etan > 0.0) then
                        psi_y_in = ybound(ii,jj-1)
                    else
                        psi_y_in = ybound(ii,jj)
                    end if

                    mat = material(ii,jj)

                    call ddsolve(mun,etan,dx(ii),dy(jj),SigmaT(mat),psi_x_in,psi_y_in,psi_x_out,psi_y_out,psi,source(ii,jj))

                    !Each sweep adds to value of scalar flux
                    phi(ii,jj) = phi(ii,jj) + w(n)*psi

                    !Define outgoing fluxes out of the cells
                    if (mun > 0.0) then
                        xbound(ii,jj) = psi_x_out
                    else
                        xbound(ii-1,jj) = psi_x_out
                    end if
                    if (etan > 0.0) then
                        ybound(ii,jj) = psi_y_out
                    else
                        ybound(ii,jj-1) = psi_y_out
                    end if
                end do
            end do

            !Save this quadrant's exiting boundary flux for the next sweep() call
            if (mun > 0.0) then
                rightflux(quad,n,:) = xbound(I,:)
            else
                leftflux(quad,n,:) = xbound(0,:)
            end if
            if (etan > 0.0) then
                topflux(quad,n,:) = ybound(:,J)
            else
                bottomflux(quad,n,:) = ybound(:,0)
            end if

        end do
    end do
end subroutine sweep