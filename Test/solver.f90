program solver_n
    implicit none
    integer :: out
    integer :: I,J,K,M,G,BCL,BCR,BCB,BCT
    real, allocatable :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:,:),SigmaA(:,:),SigmaS(:,:),nuSigmaF(:,:),Sigma1_2(:,:),phi(:,:,:)
    integer,allocatable :: material(:,:)
    integer :: count0, count1, rate
    real :: keff
    real :: diff
    real :: tol, tolk
    integer :: maxiter,maxouter
    integer :: ii,jj,gg


    call system_clock(count=count0)
    call system_clock(count_rate = rate)

    open(newunit=out, file="Test", status="replace", action = "write")
    call Version_data(out)
    call Input_data(out,I,J,K,M,G,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w,SigmaT,SigmaA,SigmaS,nuSigmaF,Sigma1_2,material,tol,maxiter)
    call Input_check(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,SigmaS,BCL,BCR,BCB,BCT,material,tol,maxiter)
    call Input_echo(out,I,J,K,G,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,SigmaS,SigmaA,nuSigmaF,Sigma1_2,material)

    allocate(phi(I,J,G))

    call Power_iteration(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,Sigma1_2,nuSigmaF,material,tol,maxiter,tolk,maxouter,phi,keff)

    write(out,'(/,a)') "Discrete Ordinates Method Solution"
    do gg=1,G
        write(out,'(/,a,i0)') "Group ", gg
        write(out,'(a3,1x,a3,1x,a26)') "i","j","Cell-Averaged Scalar Flux"
        do jj = 1,J
            do ii = 1,I
                write(out,'(i3,1x,i3,1x,es16.8)') ii,jj,phi(ii,jj,gg)
            end do
        end do
    end do
write(out,'(/,a,f12.6)') "k-effective = ", keff



    call system_clock(count=count1)
    diff = real(count1-count0)/real(rate)
    write(out,'(/,a,f10.4)') "Run Time (s):", diff

    close(out)

    contains
    subroutine Version_data(out)
        implicit none
        integer, intent(in) :: out
        character(len=32) :: execution_date, execution_time
        call date_and_time(date=execution_date,time=execution_time)
        write(out,'(a)') "Transport Solver"
        write(out,'(a)') "Version 4.0"
        write(out,'(a)') "Kirill Polokhalo"
        write(out,'(a)') "Date:", execution_date
        write(out,'(a)') "Time:", execution_time
        write(out,'(a)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    end subroutine Version_data


    subroutine Input_data(out,I,J,K,M,G,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w,SigmaT,SigmaA,SigmaS,nuSigmaF,Sigma1_2,material,tol,maxiter)
        implicit none
        integer, intent(in) :: out
        integer, intent(out) :: I,J,K,M,G
        integer,intent(out) :: BCL,BCR,BCB,BCT
        real, allocatable, intent(out) :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:,:), & 
        SigmaA(:,:),SigmaS(:,:),nuSigmaF(:,:),Sigma1_2(:,:)
        integer,allocatable,intent(out) :: material(:,:)
        real, intent(out) :: tol
        integer, intent(out) :: maxiter
        integer :: n,ii,jj,gg

        !Reading size of mesh and preallocating the matrix 
        read(*,*) I,J
        allocate(dx(I), dy(J)) 
        !Creating cells in matrix in y and x direction
        read(*,*) (dx(ii), ii=1,I)
        read(*,*) (dy(jj), jj=1,J)
        
        !Reading # of angles per octant and allocating vectors that will be of this size
        read(*,*) K 
        allocate(mu(K),eta(K),w(K))
        do n=1,K                        !We might have more than 1 angle, so we make it dependable on K. So far we only have one of each mu etc, but this should work for any #, hence the loop
            read(*,*) mu(n),eta(n),w(n)
        enddo

        !Reading # of materials & groups
        read(*,*) G
        read(*,*) M 
        allocate(SigmaT(M,G),SigmaA(M,G),nuSigmaF(M,G),Sigma1_2(M,G))   !We created the array for each X/S of length of materials, and then we fill it out in the order that our input is positioned in (T followed by S)
        do n=1,M 
            do gg=1,G
                read(*,*) SigmaT(n,gg),SigmaA(n,gg),nuSigmaF(n,gg),Sigma1_2(n,gg)
            end do
        end do

        !Create/Calculate in-group scattering cross-section
        allocate(SigmaS(M,G))
        do n=1,M 
            do gg=1,G
                SigmaS(n,gg)= SigmaT(n,gg)-SigmaA(n,gg)
            end do
        end do

        !Reading BCs
        read(*,*) BCL,BCR,BCB,BCT
        
        !Create material matrix
        allocate(material(I,J))
        do jj=1,J
            read(*,*) (material(ii,jj), ii=1,I)    !Could not do a nested loop, therefore did this trick since read function only reads horizontally all characters on the line.
        end do

        !Read convergence criteria
        read(*,*) tol,maxiter
        read(*,*) tolk,maxouter

        write(out,'(a)') "Input has been read, variables assigned, matrices created"

    end subroutine Input_data

    subroutine Input_check(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,SigmaS,BCL,BCR,BCB,BCT,material,tol,maxiter)
        implicit none
        integer :: n,ii,jj,gg
        real,intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:,:),SigmaA(:,:),SigmaS(:,:)
        integer,intent(in) :: BCL,BCR,BCB,BCT, material(:,:)
        integer, intent(in) :: out,I,J,K,M,G
        real,intent(in) :: tol
        integer,intent(in) :: maxiter
        
        !Input Dimension Check
        if (J <= 0 .or. I <= 0) then
            write(out,'(a,i0,a,f12.6)') "ERROR: I and J must be >= 1"
            stop 1
        end if
        if (K < 1) then 
             write(out,'(a,i0,a,f12.6)') "ERROR: # of angles per octant has to be >= 1"
            stop 1
        end if
        if (M < 1) then 
             write(out,'(a,i0,a,f12.6)') "ERROR: # of materials has to be >= 1"
            stop 1
        end if
        if (G < 1) then 
             write(out,'(a,i0,a,f12.6)') "ERROR: # of energy groups has to be >= 1"
            stop 1
        end if
        do jj=1,J 
            if(dy(jj) <= 0.0) then
            write(out,'(a,i0,a,f12.6)') "ERROR: dx cell size has to be > 0"
            stop 1
            end if
        end do
        do ii=1,I 
            if(dx(ii) <= 0.0) then
            write(out,'(a,i0,a,f12.6)') "ERROR: dy cell size has to be > 0"
            stop 1
            end if
        end do
        
        !Angle Physicality Check
        do n=1,K
            if (w(n) <= 0.0 .or. abs(mu(n)) > 1.0 .or. abs(eta(n)) > 1.0 ) then
            write(out,'(a,i0,a,f12.6)') "ERROR: angle component value not physical"
            stop 2
            end if
            if(eta(n)**2+mu(n)**2 > 1.0) then
            write(out,'(a,i0,a,f12.6)') "ERROR: eta^2 + mu^2 >1 "
            stop 2
            end if
        end do 

        !X/S Check
        do n=1,M 
            do gg=1,G
                if(SigmaT(n,gg) < 0.0) then
                write(out,'(a,i0,a,f12.6)') "ERROR: Total X/S less than 0"
                stop 3
                end if
                if(SigmaA(n,gg) < 0.0) then
                write(out,'(a,i0,a,f12.6)') "ERROR: Absorption X/S less than 0"
                stop 3
                end if
                if(SigmaS(n,gg) < 0.0) then
                write(out,'(a,i0,a,f12.6)') "ERROR: Scattering X/S less than 0"
                stop 3
                end if
                if (SigmaS(n,gg) > SigmaT(n,gg)) then
                write(out,'(a,i0,a,f12.6)') "ERROR: Scattering X/S greater than Total X/S"
                stop 3
                end if
            end do
        end do

        !BC Check
        if (BCL /= 0 .and. BCL /=1 ) then
        write(out,'(a,i0,a,f12.6)') "ERROR: BCL neither 0 or 1"
        stop 4
        end if
        if (BCR /= 0 .and. BCR /=1 ) then
        write(out,'(a,i0,a,f12.6)') "ERROR: BCR neither 0 or 1"
        stop 4
        end if
        if (BCB /= 0 .and. BCB /=1 ) then
        write(out,'(a,i0,a,f12.6)') "ERROR: BCB neither 0 or 1"
        stop 4
        end if
        if (BCT /= 0 .and. BCT /=1 ) then
        write(out,'(a,i0,a,f12.6)') "ERROR: BCT neither 0 or 1"
        stop 4
        end if

        do jj = 1,J
            do ii = 1,I
                if (material(ii,jj) < 1 .or. material(ii,jj)> M) then
                    write(out,'(a,i0,a,f12.6)') "ERROR: Material Type non-existant "
                    stop 5
                end if
            end do
        end do
        
        if (tol <= 0.0d0) then
            write(out,'(a)') "Error: convergence criterion must be greater than 0"
            stop 7
        end if

        if (maxiter < 1 ) then
            write(out,'(a)') "Error: maximum # of iterations has to be grater than 1"
            stop 8
        end if

    end subroutine Input_check

    subroutine Input_echo(out,I,J,K,G,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,SigmaS,SigmaA,nuSigmaF,Sigma1_2,material)
        implicit none
        integer, intent(in) :: out,I,J,K,G
        integer, intent(in) :: BCL,BCR,BCB,BCT
        real, intent(in) :: mu(:),eta(:),w(:),dx(:),dy(:),SigmaT(:,:),SigmaS(:,:),SigmaA(:,:),nuSigmaF(:,:),Sigma1_2(:,:)
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
        "i", "j", "Material", "Group", "dx", "dy", "SigmaT", "SigmaA", "SigmaS", "nuSigmaF", "Sigma1_2"
        do jj = 1,J
            do ii = 1,I
                m = material(ii,jj)
                    do gg = 1,G
                        write(out,'(i3,1x,i3,1x,i8,1x,i5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5)') &
                        ii,jj,m,gg,dx(ii),dy(jj),SigmaT(m,gg),SigmaA(m,gg),SigmaS(m,gg),nuSigmaF(m,gg),Sigma1_2(m,gg)
                    end do
            end do
        end do

    end subroutine Input_echo


    subroutine Power_iteration(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,Sigma1_2,nuSigmaF, &
                            material,tol,maxiter,tolk,maxouter,phi,keff)
    implicit none
    integer, intent(in) :: out,I,J,K,M,G,maxiter,maxouter
    real, intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:)
    real, intent(in) :: SigmaT(:,:),SigmaA(:,:),Sigma1_2(:,:),nuSigmaF(:,:)
    integer, intent(in) :: material(:,:)
    real, intent(in) :: tol,tolk
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
           maxiter,tol,phi(:,:,gg),.false.,BCL,BCR,BCB,BCT)
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

        if (diffk <= tolk) then
            write(out,'(/,a)') "Power iteration converged successfully."
            write(out,'(a,i6)') "Outer iterations consumed:", outer
            exit
        end if

        if (outer == maxouter) then
            write(out,'(/,a)') "Power iteration terminated unsuccessfully."
        end if

    end do

end subroutine Power_iteration


   subroutine inner(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,SigmaS,material,source,maxiter,tol,phi,verbose, &
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
    !tol - convergence criterion
    !BCL,BCR,BCB,BCT - boundary flags, 0 = vacuum, 1 = reflective

    !Output:
    !phi - scalar flux after inner iterations

    integer,intent(in) :: I,J,K,M,out
    real, intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:)
    real, intent(in) :: SigmaT(:),SigmaS(:),source(:,:)
    integer, intent(in) :: material(:,:)
    integer, intent(in) :: maxiter
    real, intent(in) :: tol
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

    do it = 1,maxiter
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

        if (maxdiff <= tol) then
            if (verbose) then
                write(out,'(/,a)') "Inner iterations converged successfully."
                write(out,'(a,i6)') "Iterations consumed:", it
                write(out,'(a,1x,es12.5)') "Convergence criterion achieved:", maxdiff
            end if
            exit
        end if

        if (it == maxiter) then
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

    subroutine ddsolve(mu,eta,dx,dy,Sigma_T,psi_x_in,psi_y_in,psi_x_out,psi_y_out,psi,q)
    !   Inputs:
    !   dx,dy - represent the cell width in x and y dimension
    !   mu, eta - direction cosines for discrete ordinate
    !   Sigma_T - total macroscopic cross section
    !   q - source value for this cell
    !   psi_x_in - incoming angular flux coming into the cell on x-face of cell
    !   psi_y_in - incoming angular flux coming into the cell on y-face of cell

    !   Outputs:
    !   psi - cell averaged flux
    !   psi_x_out - outgoing angular flux coming from cell on x-face of cell
    !   psi_y_out - outgoing angular flux coming from cell on y-face of cell
    
    implicit none
    real, intent(in) :: mu,eta,dx,dy,Sigma_T,psi_x_in,psi_y_in,q
    real, intent(out) :: psi_x_out, psi_y_out, psi

    psi = (q + (2.0*abs(mu)/dx)*psi_x_in + (2.0*abs(eta)/dy)*psi_y_in) &
        / (Sigma_T+(2.0*abs(mu)/dx)+(2.0*abs(eta)/dy))
    psi_x_out = 2.0*psi-psi_x_in
    psi_y_out = 2.0*psi-psi_y_in

end subroutine ddsolve 


end program solver_n