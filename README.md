# Heat2D_Per_Lib
##Introductions:##
We present a family of integral- equation-based solvers  for the heat equation, reaction-diffusion,systems, the unsteady Stokes equation and the incompressible Navier-Stokes equations in two space dimensions.
commom: the adaptive tree structure (with updated adapt refinement and coarsening routines)

bfgt&pfgt:new verison of adaptive FGT library (box version and point version)

hpots:2D inhomogeneous heat eqn solver, adaptive version

    heat_box_per.f: full solver of inhomogeneous heat eqn, with periodic bc on a unit box,including inhomogeneous heat eqn solver(heat2dperfm3),Semilinear Heat Equation(heat2d_am24),reaction-diffusion system(heat2d_sys_am24)
    
    hpots2dadap.f:  Subroutines for the evaluation of heat potentials
    
unsteady_stokes_box2d_per: unsteady stokes equation solver with periodic BC,adaptive version

    unsteady_stokes2d_box_per_order4.f: 4th order Adams-Bashforth methods for unsteady-stokes eq

ns_solver: Navier-Stokes eqs solver,adaptive version,predict-correction

    ns2d_box_per_order4.f: 2nd-4th order Adams-Bashforth methods for unsteady-stokes eq,predict-correction
    
utils: utility routines (Chebyshev routines, Legendre routines, fft, gamma function evaluation, etc.)

###How to run###
makefile

