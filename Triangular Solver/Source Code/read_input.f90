subroutine read_input(fname,out,I,J,K,G,M,h,BCL,BCR,BCH,mu,eta,w,SigmaT,SigmaA,nuSigmaF,Sigma1_2,source, &
        material,toli,tolo,div,maxinner,maxouter,K_tot) ! fname holds the name of the file is it supposed to open, e.g. 'dummy_input.txt'
    implicit none
    character(len=*), intent(in) :: fname 
    character(len=300) :: line
    character(len=20) :: card
    logical :: ifxst
    integer :: c,ii,jj,mm,o
    integer :: in,Nsn
    
    integer, intent(in) :: out
    integer, intent(out) :: I,J,G,M,K, K_tot
    integer,intent(out) :: h
    integer,intent(out) :: BCL,BCR,BCH
    real(8), allocatable, intent(out) :: mu(:),eta(:),w(:)
    real(8), allocatable, intent(out) :: SigmaT(:,:),SigmaA(:,:),nuSigmaF(:,:),Sigma1_2(:,:),source(:,:)
    integer,allocatable,intent(out) :: material(:,:)
    real(8), intent(out) :: toli,tolo,div
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
        elseif (card.eq.'h') then
            read(line,*) card, h
        elseif (card.eq.'Nsn') then
            read(line,*) card, Nsn
            call hexsncon(Nsn,mu,eta,w,K_tot)  ! fills first K entries, weights sum to 0.25
            K=K_tot/6
        elseif (card.eq.'div') then
            read(line,*) card, div
        elseif (card.eq.'bc') then
            read(line,*) card, BCL,BCR,BCH !BC Left, BC Right, BC Horizontal
        elseif (card.eq.'ngroup') then
            read(line,*) card, G
        elseif (card.eq.'nmat') then
            read(line,*) card, M
            allocate(SigmaT(M,G),SigmaA(M,G),nuSigmaF(M,G),Sigma1_2(M,G),source(M,G))
        elseif (card.eq.'siga') then
            read(line,*) card, o, (SigmaA(mm,o), mm=1,M)
        elseif (card.eq.'nusigf') then
            read(line,*) card, o, (nuSigmaF(mm,o), mm=1,M)
        elseif (card.eq.'sigt') then
            read(line,*) card, o, (SigmaT(mm,o), mm=1,M)
        elseif (card.eq.'sig1_2') then
            read(line,*) card, o, (Sigma1_2(mm,o), mm=1,M)
        elseif (card.eq.'source') then
            read(line,*) card, o, (source(mm,o), mm=1,M)
        elseif (card.eq.'matmap') then
            allocate(material(I,J))
            do ii=1,I
                read(in,*) (material(ii,jj), jj=1,J)
            enddo
        elseif(card.eq.'tol') then
            read(line,*) card, tolo, toli
        elseif(card.eq.'maxiter') then
            read(line,*) card, maxouter, maxinner
        endif
    enddo
    300 continue 
    close(in)

    
    write(out,'(a)') "Input has been read, variables assigned, matrices created"

    write(out,'(a,i0)') "Hexagonal Quadrature: S",Nsn
    write(out,'(a,i0)') "Directions per sextant: ",K
    write(out,'(a,i0)') "Total angular directions: ",K_tot
    write(out,'(a,f12.8)') "Sum of quadrature weights: ",sum(w)

    end subroutine read_input