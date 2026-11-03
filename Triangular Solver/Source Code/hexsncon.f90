subroutine hexsncon()
    implicit none
    !Inputs
    integer, intent(in) :: Nsn
    !Outputs
    !Local

    select case (Nsn)
    
    case(2)
        xi(1)=0.5773503
        ang(1)=pi/6.0d0
        wt(1)=1.0d0
    case(4)
        xi(1:3)= (/0.3399810d0,0.3399810d0,0.8611363d0,/)
        ang(1:3) = (/pi/12.0d0,pi/4.0d0,pi/6.0d0/)
        wt(1:3) = (/0.3260726d0,0.3260726d0,0.3478548d0/)






end subroutine hexsncon