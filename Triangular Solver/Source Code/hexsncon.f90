subroutine hexsncon(Nsn,mu,eta,wt,K_tot)
    implicit none
    !Inputs
    integer, intent(in) :: Nsn
    !Outputs
    real(8), allocatable, intent(out) :: mu(:), eta(:), wt(:)
    integer, intent(out) :: K_tot
    !Local
    real(8), allocatable :: xi(:), ang(:), w(:)
    integer :: i,n,s,K
    real(8) :: pi,gamma

    pi=acos(-1.0d0)
    K = Nsn*(Nsn+2)/8
    K_tot = 6*K

    allocate (xi(K),ang(K),w(K))
    allocate (mu(K_tot),eta(K_tot), wt(K_tot))

    select case (Nsn)
    
    case(2) !K=1
        xi(1)=0.5773503d0
        ang(1)=pi/6.0d0
        w(1)=1.0d0
    case(4) !K=3
        xi(1:K)= (/0.3399810d0,0.3399810d0,0.8611363d0/)
        ang(1:K) = (/pi/12.0d0,pi/4.0d0,pi/6.0d0/)
        w(1:K) = (/0.3260726d0,0.3260726d0,0.3478548d0/)
    case(6) !K=6
        xi(1:K) = (/0.2386192d0,0.2386192d0,0.2386192d0,0.6612094d0,0.6612094d0,0.9324695d0/)
        ang(1:K) = (/pi/18.0d0,pi/6.0d0,5.0d0*pi/18.0d0,pi/12.0d0,pi/4.0d0,pi/6.0d0/)
        w(1:K) = (/0.1559713d0,0.1559713d0,0.1559713d0,0.1803808d0,0.1803808d0,0.1713245d0/)
    case(8)
        xi(1:K) = (/0.1834346d0,0.1834346d0,0.1834346d0,0.1834346d0,0.5255324d0,0.5255324d0, &
                    0.5255324d0,0.7966665d0,0.7966665d0,0.9602899d0/)
        ang(1:K) = (/pi/24.0d0,pi/8.0d0,5.0d0*pi/24.0d0,7.0d0*pi/24.0d0,pi/18.0d0,pi/6.0d0, &
                    5.0d0*pi/18.0d0,pi/12.0d0,pi/4.0d0,pi/6.0d0/)
        w(1:K) = (/0.0906709d0,0.0906709d0,0.0906709d0,0.0906709d0,0.1045689d0,0.1045689d0, &
                    0.1045689d0,0.1111905d0,0.1111905d0,0.1012285d0/)
    end select

    do s=0,5
        do i=1,K
            n=s*K+i
            gamma = ang(i)+s*pi/3.0d0
            mu(n) = sqrt(1.0d0-xi(i)**2)*cos(gamma)
            eta(n) = sqrt(1.0d0-xi(i)**2)*sin(gamma)
            wt(n) = w(i)/6.0d0
        enddo
    enddo

end subroutine hexsncon