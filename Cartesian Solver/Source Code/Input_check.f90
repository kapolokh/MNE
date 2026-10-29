subroutine Input_check(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,BCL,BCR,BCB,BCT,material,toli,maxinner)
        implicit none
        integer :: n,ii,jj,gg
        real,intent(in) :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:,:),SigmaA(:,:)
        integer,intent(in) :: BCL,BCR,BCB,BCT, material(:,:)
        integer, intent(in) :: out,I,J,K,M,G
        real,intent(in) :: toli
        integer,intent(in) :: maxinner
        
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
        
        if (toli <= 0.0d0) then
            write(out,'(a)') "Error: convergence criterion must be greater than 0"
            stop 7
        end if

        if (maxinner < 1 ) then
            write(out,'(a)') "Error: maximum # of iterations has to be grater than 1"
            stop 8
        end if

    end subroutine Input_check