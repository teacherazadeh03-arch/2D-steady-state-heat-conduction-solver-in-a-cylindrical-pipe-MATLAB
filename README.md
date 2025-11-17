# 2D-steady-state-heat-conduction-solver-in-a-cylindrical-pipe-MATLAB-Gauss-Seidel iterative method
This MATLAB code is a 2D steady-state heat conduction solver in a cylindrical pipe, including 3D visualization of temperature, axial gradients, and radial heat flux on the pipe wall. Here’s the detailed explanation:
Assumptions
Steady-state heat conduction (no time-dependence).
Cylindrical pipe with uniform thermal conductivity (k).
Neglects internal heat generation inside the material.
Symmetry in angular direction (θ), so 2D r–z model is sufficient.
Linearized convection for outer wall (Robin BC).
Optional inner wall convection (controlled by use_inner_convection).

What the code does

1- Solves 2D steady-state heat conduction in cylindrical coordinates (r–z plane).
Uses finite difference discretization in radial (r) and axial (z) directions.
Iterates with Gauss-Seidel method until convergence.
2- Boundary conditions handled:
Inner wall (r = ri): Can be either insulated or convective (use_inner_convection flag).
Outer wall (r = ro): Convective (Robin boundary) with heat transfer coefficient h and ambient temperature Tinf.
Axial ends (z = 0 and z = L): Dirichlet BCs, hot temperature at inlet (Thot) and cold at outlet (Tcold).
3-After convergence, it calculates:
Radial heat flux at outer wall
Axial thermal gradient on outer wall: 
Maps these to a 3D cylinder for visualization using surf.
4-Plots include:
Outer wall temperature.
Axial temperature gradient (
Radial heat flux on the outer wall.
Full r–z temperature field.
