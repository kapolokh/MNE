program solver_n
    implicit none
    integer :: out
    integer :: I,J,K,M,BCL,BCR,BCB,BCT
    real, allocatable :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:),Sigma1_2(:),source(:,:)
    integer,allocatable :: material(:,:)
    integer :: count0, count1, rate
    real :: diff
    real :: tol
    integer :: maxiter

    call system_clock(count=count0)
    call system_clock(count_rate = rate)

    open(newunit=out, file="Sn_bwr_Output", status="replace", action = "write")
    call Version_data(out)
    call Input_data(out,I,J,K,M,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w,SigmaT,Sigma1_2,material,source,tol,maxiter)
    call Input_check(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,Sigma1_2,BCL,BCR,BCB,BCT,material,source,tol,maxiter)
    call Input_echo(out,I,J,K,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,Sigma1_2,source,material)
    call Transport_solver(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,Sigma1_2,material,source,tol,maxiter)
    
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


    subroutine Input_data(out,I,J,K,M,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w,SigmaT,Sigma1_2,material,source,tol,maxiter)
        implicit none
        integer, intent(in) :: out
        integer, intent(out) :: I,J,K,M
        integer,intent(out) :: BCL,BCR,BCB,BCT
        real, allocatable, intent(out) :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:),Sigma1_2(:)
        integer,allocatable,intent(out) :: material(:,:)
        real,allocatable,intent(out) :: source(:,:)   !Source is a decimal value
        real, intent(out) :: tol
        integer, intent(out) :: maxiter
        integer :: n,ii,jj

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

        !Reading # of materials
        read(*,*) M 
        allocate(SigmaT(M),Sigma1_2(M))   !We created the array for each X/S of length of materials, and then we fill it out in the order that our input is positioned in (T followed by S)
        do n=1,M 
            read(*,*) SigmaT(n),Sigma1_2(n)
        end do

        !Reading BCs
        read(*,*) BCL,BCR,BCB,BCT
        

        !Create material matrix
        allocate(material(I,J),source(I,J))
        do jj=1,J
            read(*,*) (material(ii,jj), ii=1,I)    !Could not do a nested loop, therefore did this trick since read function only reads horizontally all characters on the line.
        end do

        !Creating source matrix
        do jj=1,J
            read(*,*) (source(ii,jj), ii=1,I)        
        end do

        !Read convergence criteria
        read(*,*) tol,maxiter

        write(out,'(a)') "Input has been read, variables assigned, matrices created"

    end subroutine Input_data

    subroutine Input_check(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,Sigma1_2,BCL,BCR,BCB,BCT,material,source,tol,maxiter)
        implicit none
        integer :: n,ii,jj
        real,intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:),Sigma1_2(:), source(:,:)
        integer,intent(in) :: BCL,BCR,BCB,BCT, material(:,:)
        integer, intent(in) :: out,I,J,K,M
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
            if(SigmaT(n) < 0.0) then
            write(out,'(a,i0,a,f12.6)') "ERROR: Total X/S less than 0"
            stop 3
            end if
            if(Sigma1_2(n) < 0.0) then
            write(out,'(a,i0,a,f12.6)') "ERROR: Scattering X/S less than 0"
            stop 3
            end if
            if (Sigma1_2(n) > SigmaT(n)) then
            write(out,'(a,i0,a,f12.6)') "ERROR: Scattering X/S greater than Total X/S"
            stop 3
            end if
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
                if (source(ii,jj) < 0.0) then
                    write(out,'(a,i0,a,f12.6)') "ERROR: Source value < 0 "
                    stop 6
                end if
            end do
        end do
        
        if (tol <= 0.0d0) then
            write(out,'(a)') "Error: convergence criterion must be greater than 0"
            stop 7
        end if

        if (maxiter < 1 ) then
            write(out,'(a)') "Error: maximum # of iterations ahs to be grater than 1"
            stop 8
        end if

    end subroutine Input_check

    subroutine Input_echo(out,I,J,K,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,Sigma1_2,source,material)
        implicit none
        integer, intent(in) :: out,I,J,K
        integer, intent(in) :: BCL,BCR,BCB,BCT
        real, intent(in) :: mu(:),eta(:),w(:),dx(:),dy(:),source(:,:),SigmaT(:),Sigma1_2(:)
        integer,intent(in) :: material(:,:)
        integer :: n,ii,jj,m
        !We print out some of the necessary input variavles to output file
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
        write(out,'(a3,1x,a3,1x,a8,1x,a10,1x,a10,1x,a10,1x,a10,1x,a10)') &
         "i", "j", "Material", "dx", "dy","SigmaT","Sigma1_2","Source"
        do jj = 1,J 
            do ii = 1,I 
                m = material(ii,jj)
                write(out,'(i3,1x,i3,1x,i8,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5,1x,f10.5)') &
                 ii,jj,m,dx(ii),dy(jj),SigmaT(m),Sigma1_2(m),source(ii,jj)
            end do
        end do

    end subroutine Input_echo

    subroutine Transport_solver(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,Sigma1_2,material,source,tol,maxiter)
        implicit none
        integer,intent(in) :: out
        integer, intent(in) :: I,J,K,M,maxiter
        real, intent(in) :: dx(:), dy(:), mu(:), eta(:), w(:), SigmaT(:), Sigma1_2(:), source(:,:), tol
        integer, intent(in) :: material(:,:)
        real, allocatable :: phi (:,:)
        allocate(phi(I,J))
        write(out,'(/a)') "will solve D.O. here"

        call inner(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,Sigma1_2,material,source,maxiter,tol,phi)

    end subroutine Transport_solver

    subroutine inner(out,I,J,K,M,dx,dy,mu,eta,w,SigmaT,Sigma1_2,material,source,maxiter,tol,phi)
        implicit none
        !Input:
        !I,J - number of spatial cells
        !K - # of angles per octant
        !M - # of materials
        !dx,dy - width of the cell
        !mu,eta - ordiantes
        !w - weight
        !SigmaT - total X/S
        !Sigma1_2 - scattering X/S
        !material - material map
        !source - fixed source in each material cell
        !tol - convergence criterion

        !Output:
        !phi - scalar flux after inner iterations

        integer,intent(in) :: I,J,K,M,out
        real, intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:)
        real, intent(in) :: SigmaT(:),Sigma1_2(:),source(:,:)
        integer, intent(in) :: material(:,:)
        integer, intent(in) :: maxiter
        real, intent(in) :: tol
        real, intent(out) :: phi(:,:)
        real,allocatable :: phi_old(:,:), q(:,:)
        real :: diff,maxdiff
        integer :: ii,jj,mm,it

        allocate (phi_old(I,J) , q(I,J))

        phi_old = 0.0d0
        phi = 0.0d0

        do it = 1,maxiter
            do jj = 1,J
                do ii = 1,I
                    mm = material(ii,jj)
                    q(ii,jj) = source(ii,jj) + Sigma1_2(mm)*phi_old(ii,jj)
                end do
            end do

            call sweep(I,J,K,M,dx,dy,mu,eta,w,SigmaT,material,q,phi)

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
                write(out,'(/,a)') "Inner iterations converged successfully."
                write(out,'(a,i6)') "Iterations consumed:", it
                write(out,'(a,1x,es12.5)') "Convergence criterion achieved:", maxdiff
                exit
            end if

            if (it == maxiter) then
                write(out,'(/,a)') "Inner iterations terminated unsuccessfully."
                write(out,'(a,i6)') "Iterations consumed:", it
                write(out,'(a,1x,es12.5)') "Convergence criterion achieved:", maxdiff
                exit
            end if

            phi_old = phi
        end do

        write(out,'(/,a)') "Discrete Ordinates Method Solution"
        write(out,'(a3,1x,a3,1x,a26)') "i","j","Cell-Averaged Scalar Flux"

        do jj = 1,J
            do ii = 1,I
                write(out,'(i3,1x,i3,1x,es16.8)') ii,jj,phi(ii,jj)
            end do
        end do
        
    end subroutine inner 


    subroutine sweep(I,J,K,M,dx,dy,mu,eta,w,SigmaT,material,source,phi)
        implicit none
        integer, intent(in) :: I,J,K,M
        real, intent(in) :: dx(:), dy(:), mu(:), eta(:), w(:), SigmaT(:), source(:,:)
        integer, intent(in) :: material(:,:)
        real, intent(out) :: phi(:,:)
        real,allocatable :: xbound(:,:), ybound(:,:)
        real :: mun, etan
        real :: psi
        real :: psi_x_in
        real :: psi_y_in
        real :: psi_x_out
        real :: psi_y_out
        integer :: mat 
        integer :: quad
        integer :: j_start, j_end, j_step,i_start, i_end, i_step
        integer :: n,ii,jj
        
        allocate(xbound(0:I,J),ybound(I,0:J))
        phi = 0.0

        !Adjust ordiante signs depending on direction of the sweep
        !Do sweep in 4 directions to get all particles travelling in all directions in mesh.
        do quad = 1,4
            do n = 1,K
                select case(quad)
                case(1)
                    mun = mu(n)
                    etan = eta(n)
                case(2)
                    mun = -mu(n)
                    etan = eta(n)
                case(3)
                    mun = mu(n)
                    etan = -eta(n)
                case(4)
                    mun = -mu(n)
                    etan = -eta(n)
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