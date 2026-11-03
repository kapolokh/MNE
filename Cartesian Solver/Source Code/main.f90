program main
    implicit none
    integer :: mkdir_status
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

    call execute_command_line('mkdir -p ../Output')

    open(newunit=out, file="../Output/Final.txt", status="replace", action="write")
    call Version_data(out)
    fname = '../Input/base_input.txt'
    call read_input(fname,out,in,I,J,K,G,M,h,BCL,BCR,BCB,BCT,dx,dy,mu,eta,w, & 
                            SigmaT,SigmaA,nuSigmaF,Sigma1_2,material,toli,tolo,maxinner,maxouter)                        
    call refine_grid(I,J,h,out,material,dx,dy)
    call Input_check(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,BCL,BCR,BCB,BCT,material,toli,maxinner)
    call Input_echo(out,I,J,K,G,mu,eta,w,BCL,BCR,BCB,BCT,dx,dy,SigmaT,SigmaA,nuSigmaF,Sigma1_2,material)

    allocate(phi(I,J,G))

    call Power_iteration(out,I,J,K,M,G,dx,dy,mu,eta,w,SigmaT,SigmaA,Sigma1_2,nuSigmaF,material,toli,maxinner, &
        tolo,maxouter,phi,keff,BCB,BCR,BCL,BCT)

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

    call write_vtk('../Output/flux.vtk',out,I,J,G,dx,dy,phi,material)

    call system_clock(count=count1)
    diff = real(count1-count0)/real(rate)
    write(out,'(/,a,f10.4)') "Run Time (s):", diff

    close(out)

    contains

    include 'Version_data.f90'
    include 'read_input.f90'
    include 'refine_grid.f90'
    include 'Input_check.f90'
    include 'Input_echo.f90'
    include 'Power_Iteration.f90'
    include 'inner.f90'
    include 'sweep.f90'
    include 'ddsolve.f90'
    include 'tsncon.f90'
    include 'write_vtk.f90'

end program main