program solver_n
    implicit none
    integer :: out, in
    integer :: I,J,K,M,G,BCL,BCR,BCB,BCT
    real, allocatable :: dx(:),dy(:),mu(:),eta(:),w(:),SigmaT(:,:),SigmaA(:,:),SigmaS(:,:),nuSigmaF(:,:),Sigma1_2(:,:),phi(:,:,:)
    integer,allocatable :: material(:,:)
    integer :: count0, count1, rate
    real :: keff
    real :: diff
    real :: toli, tolo
    integer :: h
    integer :: maxinner,maxouter
    integer :: ii,jj,gg
    character(len=100) :: fname


    call system_clock(count=count0)
    call system_clock(count_rate = rate)

    open(newunit=out, file="Final", status="replace", action = "write")
    call Version_data(out)
    fname = 'dummy_input.txt'
    call read_input(fname,out,in,I,J,K,G,M,h,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w, & 
                            SigmaT,SigmaA,nuSigmaF,Sigma1_2,material,toli,tolo,maxinner,maxouter)                        
    call refine_grid(I,J,h,out,material,dx,dy)
    call Input_check(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,BCL,BCR,BCB,BCT,material,toli,maxinner)
    call Input_echo(out,I,J,K,G,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,SigmaA,nuSigmaF,Sigma1_2,material)

    allocate(phi(I,J,G))

    call Power_iteration(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,Sigma1_2,nuSigmaF,material,toli,maxinner,tolo,maxouter,phi,keff)

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
    write(out,'(/,a,f16.10)') "k-effective = ", keff

    call write_vtk('flux.vtk',out,I,J,G,dx,dy,phi,material)

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

    subroutine read_input(fname,out,in,I,J,K,G,M,h,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w, &
                    SigmaT,SigmaA,nuSigmaF,Sigma1_2,material,toli,tolo,maxinner,maxouter) ! fname holds the name of the file is it supposed to open, e.g. 'dummy_input.txt'
    implicit none
    character(len=*), intent(in) :: fname 
    character(len=300) :: line
    character(len=20) :: card
    logical :: ifxst
    integer :: c,ii,jj,mm,o
    integer :: in,Nsn
    real(8) :: u8(36), e8(36), w8(36)

    integer, intent(in) :: out
    integer, intent(out) :: I,J,K,G,M
    integer,intent(out) :: h
    integer,intent(out) :: BCL,BCR,BCB,BCT
    real, allocatable, intent(out) :: dx(:),dy(:),mu(:),eta(:),w(:)
    real, allocatable, intent(out) :: SigmaT(:,:),SigmaA(:,:),nuSigmaF(:,:),Sigma1_2(:,:)
    integer,allocatable,intent(out) :: material(:,:)
    real, intent(out) :: toli,tolo
    integer, intent(out) :: maxinner,maxouter

   !--- check if input file exists
    inquire(file=fname, exist=ifxst)
    if (.not.ifxst) stop 'input does not exist'

   !--- open input file for reading
    open(newunit=in,file=fname,status='old', action='read')

    do 
        read(in,'(a)',end=300) line 

        c=index(line,'#')   !reads at what index in the line # lies
        if (c.gt.0) line(c:)=''  !if line does not start with #, it means that there are comments after, so string is replaced by empty after that #
        if (line.eq.' ') cycle

        read(line,*) card   !reads the first word before the blank

        if (card.eq.'ncell') then
            read (line,*) card, I, J
            allocate(dx(I),dy(J))
        elseif (card.eq.'xdim') then
            read(line,*) card, (dx(ii), ii=1,I)
        elseif (card.eq.'ydim') then
            read(line,*) card, (dy(jj), jj=1,J)
        elseif (card.eq.'Nsn') then
            read(line,*) card, Nsn
            K = Nsn*(Nsn+2)/8            ! directions per quadrant (S8 -> 10)

            u8 = 0.0d0
            e8 = 0.0d0
            w8 = 0.0d0
            call tsncon(u8,e8,w8,Nsn,K)  ! fills first K entries, weights sum to 0.25
            allocate(mu(K),eta(K),w(K))
            mu  = real(u8(1:K))
            eta = real(e8(1:K))
            w   = real(w8(1:K))
        elseif (card.eq.'h') then
            read(line,*) card, h
        elseif (card.eq.'bc') then
            read(line,*) card, BCL,BCR,BCB,BCT
        elseif (card.eq.'ngroup') then
            read(line,*) card, G
        elseif (card.eq.'nmat') then
            read(line,*) card, M
            allocate(SigmaT(M,G),SigmaA(M,G),nuSigmaF(M,G),Sigma1_2(M,G))
        elseif (card.eq.'siga') then
            read(line,*) card, o, (SigmaA(mm,o), mm=1,M)
        elseif (card.eq.'nusigf') then
            read(line,*) card, o, (nuSigmaF(mm,o), mm=1,M)
        elseif (card.eq.'sigt') then
            read(line,*) card, o, (SigmaT(mm,o), mm=1,M)
        elseif (card.eq.'sig1_2') then
            read(line,*) card, o, (Sigma1_2(mm,o), mm=1,M)
        elseif (card.eq.'matmap') then
            allocate(material(I,J))
            do jj=1,J
                read(in,*) (material(ii,jj), ii=1,I)
            enddo
        elseif(card.eq.'tol') then
            read(line,*) card, toli, tolo
        elseif(card.eq.'maxiter') then
            read(line,*) card, maxinner, maxouter
        endif
    enddo
    300 continue 
    close(in)

    write(out,'(a)') "Input has been read, variables assigned, matrices created"
    write(out,'(a,i0,a,i0,a)') "Quadrature: S", Nsn, " (", K, " directions per quadrant)"
    end subroutine read_input




      subroutine tsncon (u,e,w,n,mm)
      implicit none 

      integer, intent(in)  :: n
      integer, intent(in)  :: mm
      real(8), intent(out) :: u(36), e(36), w(36)

!--- local

      integer :: i, j, k, m
      integer :: l, la
      real(8) :: t
!
!---  do the usual Twotran Quadrature set
!
      m=(n*(n+2))/8
      i=n/2

      if (n.eq.2) then
        w(1)=1.0d0
        u(1)=1.0d0/sqrt(3.0d0)    ! .57735027d+0
        e(1)=u(1)
        go to 240
      elseif (n.eq.4) then
        u(2)=0.90444905d+0
        u(1)=0.30163878d+0
        w(1)=0.33333333d+0
        w(2)=w(1)
        w(3)=w(1)
      elseif (n.eq.6) then
        u(1)=0.23009194d+0
        u(2)=0.68813432d+0
        u(3)=0.94557676d+0
        w(1)=0.16944656d+0
        w(2)=0.16388677d+0
        w(3)=w(1)
        w(4)=w(2)
        w(5)=w(2)
        w(6)=w(1)
      elseif (n.eq.8) then
        u(1)=0.19232747d+0
        u(2)=0.57735027d+0
        u(3)=0.79352178d+0
        u(4)=0.96229948d+0
        w(1)=0.11678847d+0
        w(2)=0.09325523d+0
        w(3)=w(2)
        w(4)=w(1)
        w(5)=w(2)
        w(6)=0.09010320d+0
        w(7)=w(2)
        w(8)=w(2)
        w(9)=w(2)
        w(10)=w(1)
      elseif (n.eq.10) then
        u(1)=0.16962228d+0
        u(2)=0.50714192d+0
        u(3)=0.69686020d+0
        u(4)=0.84500612d+0
        u(5)=0.97080202d+0
        w(1)=0.089842043d+0
        w(2)=0.067288705d+0
        w(3)=0.055780071d+0
        w(4)=w(2)
        w(5)=w(1)
        w(6)=w(2)
        w(7)=0.053133809d+0
        w(8)=w(7)
        w(9)=w(2)
        w(10)=w(3)
        w(11)=w(7)
        w(12)=w(3)
        w(13)=w(2)
        w(14)=w(2)
        w(15)=w(1)
      elseif (n.eq.12) then
        u(1)=0.15395746d+0
        u(2)=0.45769112d+0
        u(3)=0.62869660d+0
        u(4)=0.76225828d+0
        u(5)=0.87568027d+0
        u(6)=0.97600932d+0
        w(1)=0.07332178d+0
        w(2)=0.05266740d+0
        w(3)=0.04161495d+0
        w(4)=w(3)
        w(5)=w(2)
        w(6)=w(1)
        w(7)=w(2)
        w(8)=0.03895667d+0
        w(9)=0.03249018d+0
        w(10)=w(8)
        w(11)=w(2)
        w(12)=w(3)
        w(13)=w(9)
        w(14)=w(9)
        w(15)=w(3)
        w(16)=w(3)
        w(17)=w(8)
        w(18)=w(3)
        w(19)=w(2)
        w(20)=w(2)
        w(21)=w(1)
      elseif (n.eq.14) then
        u(1)=0.14238965d+0
        u(2)=0.42048076d+0
        u(3)=0.57735027d+0
        u(4)=0.69990185d+0
        u(5)=0.80398498d+0
        u(6)=0.89605866d+0
        u(7)=0.97951538d+0
        w(1)=0.062171628d+0
        w(2)=0.043325697d+0
        w(3)=0.033217605d+0
        w(4)=0.031837060d+0
        w(5)=w(3)
        w(6)=w(2)
        w(7)=w(1)
        w(8)=w(2)
        w(9)=0.030486324d+0
        w(10)=0.024545116d+0
        w(11)=w(10)
        w(12)=w(9)
        w(13)=w(2)
        w(14)=w(3)
        w(15)=w(10)
        w(16)=0.019984453d+0
        w(17)=w(10)
        w(18)=w(3)
        w(19)=w(4)
        w(20)=w(10)
        w(21)=w(10)
        w(22)=w(4)
        w(23)=w(3)
        w(24)=w(9)
        w(25)=w(3)
        w(26)=w(2)
        w(27)=w(2)
        w(28)=w(1)
      elseif (n.eq.16) then
        u(1)=0.13344572d+0
        u(2)=0.39119433d+0
        u(3)=0.53689687d+0
        u(4)=0.65075610d+0
        u(5)=0.74746822d+0
        u(6)=0.83302700d+0
        u(7)=0.91058181d+0
        u(8)=0.98203079d+0
        w(1)=0.05415425d+0
        w(2)=0.03679653d+0
        w(3)=0.02777273d+0
        w(4)=0.02580284d+0
        w(5)=w(4)
        w(6)=w(3)
        w(7)=w(2)
        w(8)=w(1)
        w(9)=w(2)
        w(10)=0.02494275d+0
        w(11)=0.01962325d+0
        w(12)=0.01879762d+0
        w(13)=w(11)
        w(14)=w(10)
        w(15)=w(2)
        w(16)=w(3)
        w(17)=w(11)
        w(18)=0.01544801d+0
        w(19)=w(18)
        w(20)=w(11)
        w(21)=w(3)
        w(22)=w(4)
        w(23)=w(12)
        w(24)=w(18)
        w(25)=w(12)
        w(26)=w(4)
        w(27)=w(4)
        w(28)=w(11)
        w(29)=w(11)
        w(30)=w(4)
        w(31)=w(3)
        w(32)=w(10)
        w(33)=w(3)
        w(34)=w(2)
        w(35)=w(2)
        w(36)=w(1)
      else
        write(*,'(a,i0)') 'FATAL: No such Sn constant set available, n = ', n
        stop 2
      endif
!
      k=i+1
      do j=2,i
        la=i+1-j
        do l=1,la
          u(k)=u(l)
          k=k+1
        enddo
      enddo
      k=1
      do j=1,i
        la=i+1-j
        do l=1,la
          e(k)=u(j)
          k=k+1
        enddo
      enddo
  240 continue
      do j=1,m
        w(j)=.25d0*w(j)
      enddo
!
!     reorder sn library
!
      k=1
      l=mm
  260 t=u(k)
      u(k)=u(l)
      u(l)=t
      t=e(k)
      e(k)=e(l)
      e(l)=t
      t=w(k)
      w(k)=w(l)
      w(l)=t
      k=k+1
      l=l-1
      if (k.lt.l) go to 260
!
      return
      end subroutine tsncon

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

        subroutine write_vtk(vtkname,out,I,J,G,dx,dy,phi,material)
        !Writes the cell-averaged scalar flux (every group) and the material
        !map to a legacy VTK rectilinear grid file for VisIt / ParaView.
        !Based on editvtkpn (Palmtag).
        implicit none
        character(len=*), intent(in) :: vtkname
        integer, intent(in) :: out, I, J, G
        real, intent(in) :: dx(:), dy(:), phi(:,:,:)
        integer, intent(in) :: material(:,:)
        integer :: vtk, ii, jj, gg
        real :: xx

        write(out,'(/,2a)') "Writing flux to VTK file: ", trim(vtkname)
        open(newunit=vtk, file=vtkname, status='replace', action='write')

        !--- VTK header
        write(vtk,'(a)') '# vtk DataFile Version 2.0'
        write(vtk,'(a)') 'Sn transport solver - cell-averaged scalar flux'
        write(vtk,'(a)') 'ASCII'
        write(vtk,'(a)') 'DATASET RECTILINEAR_GRID'
        write(vtk,'(a,i0,1x,i0,a)') 'DIMENSIONS ', I+1, J+1, ' 1'

        !--- cell edge coordinates in x (running sum of dx)
        write(vtk,'(a,i0,a)') 'X_COORDINATES ', I+1, ' float'
        xx = 0.0
        write(vtk,'(es15.7)') xx
        do ii = 1,I
            xx = xx + dx(ii)
            write(vtk,'(es15.7)') xx
        end do

        !--- cell edge coordinates in y (running sum of dy)
        write(vtk,'(a,i0,a)') 'Y_COORDINATES ', J+1, ' float'
        xx = 0.0
        write(vtk,'(es15.7)') xx
        do jj = 1,J
            xx = xx + dy(jj)
            write(vtk,'(es15.7)') xx
        end do

        !--- single z plane -> 2D
        write(vtk,'(a)') 'Z_COORDINATES 1 float'
        write(vtk,'(a)') '0.0'

        !--- cell data: one scalar field per energy group
        write(vtk,'(a,i0)') 'CELL_DATA ', I*J
        do gg = 1,G
            write(vtk,'(a,i0,a)') 'SCALARS flux_g', gg, ' float'
            write(vtk,'(a)') 'LOOKUP_TABLE default'
            do jj = 1,J
                do ii = 1,I
                    write(vtk,'(es15.7)') phi(ii,jj,gg)
                end do
            end do
        end do

        !--- material map (useful for checking geometry/orientation)
        write(vtk,'(a)') 'SCALARS material int'
        write(vtk,'(a)') 'LOOKUP_TABLE default'
        do jj = 1,J
            do ii = 1,I
                write(vtk,'(i0)') material(ii,jj)
            end do
        end do

        close(vtk)
    end subroutine write_vtk



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


    subroutine Power_iteration(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,Sigma1_2,nuSigmaF, &
                            material,toli,maxinner,tolo,maxouter,phi,keff)
    implicit none
    integer, intent(in) :: out,I,J,K,M,G,maxinner,maxouter
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