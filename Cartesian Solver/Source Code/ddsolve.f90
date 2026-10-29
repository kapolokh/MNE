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