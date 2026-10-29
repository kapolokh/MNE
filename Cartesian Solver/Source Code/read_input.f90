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