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