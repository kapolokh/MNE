Compiling the code: 1) Switch to directory containing the Code 3 files.
                    2) Compile with: gfortran solver_n.f90 -o solver.exe

Executing the code: 
Execute original input: 
                        1) In solver_n.f90 change the file name that is opened to Code4_Orig_Output.
                        2) Run ./solver.exe < input.txt

Status of the code: Operational

Expected Input:
The code establishes the framework for solving the transport equation using the discrete ordinates approximation. 
The expected input is the number of computational cells in x and y direction (I,J) as well as the size of the cell (dx,dy).
The program reads the number of angles per octant (K) and assigns values for mu, eta and w variables. 
By reading the number of materials (M) it is possible to create a material map and assign appropriate total and scattering cross sections per each material. 
The boundary conditions are defined as such, where 0 means vacuum and 1 means reflective.
Using desired number of computational cells, the material and source maps are created.

Resulting output:
1) Version information
2) Input has been successfully read and validated
3) Discrete Ordinates/Octants values printed
4) Boundary Conditions printed
5) Information with properties of each cell such as i,j,material,dx,dy,SigmaT,SigmaS and source are printed
6) Discrete Ordinate Method Solution (featuring scalar flux values)
7) Run Time of the code printed

Limitations:
Input structure has to be maintained. The BC used are only vacuum and are not designed/tested for reflective boundary conditions.























