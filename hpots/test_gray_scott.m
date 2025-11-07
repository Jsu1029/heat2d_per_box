% this test program calls chebfun to test 
% the solution of gray-scott equation
% to find suitable parameters and initial conditions
% to demonstrate the neccessity of adaptivity
%
% chebfun/spin2:
% Solve stiff PDEs in 2D periodic domains, Fourier spectral method and 
% exponential integrators.
%
% chebfun/spinop2:
% SPINOP2   Class for representing the spatial part of 2D differential operators 
% for time-dependent PDEs.
% SPINOP2 is a class for representing the spatial part S of a time-dependent 
% PDE of the form u_t = S(u) = Lu + N(u) in 2D, where L is a linear 
% operator and N is a nonlinear operator. 
%
% solves the G-S equation
% u_t = ep1*Laplace(u)+b(1-u)-uv^2
% v_t =ep2*Laplace(v)-dv+uv^2
% with periodic boundary condition
%
% Fourier spectral method in space, exponential integrator in time 
%
ep1 = 0.00002; ep2 = 0.00001;
b = 0.04; d = 0.1;
dom = [-1 1 -1 1]; 
x = chebfun('x',dom(1:2)); 
tspan = [0 3500];
%
% initialize the spinop2 operator struct
S = spinop2(dom,tspan);
% linear and nonlinear part of the spatial operator
S.lin = @(u,v) [ep1*lap(u); ep2*lap(v)];
S.nonlin = @(u,v) [b*(1-u)-u.*v.^2;-d*v+u.*v.^2];
% initial condition
% chebfun2v: 2D, vector
S.init = chebfun2v(@(x,y) 1-exp(-80*((x+.05).^2+(y+.02).^2)), ...
                   @(x,y) exp(-80*((x-.05).^2+(y-.02).^2)),dom);
% call spin2 to solve, make sure you understand the interface
% u= spin2(S, N, dt) where N specifies the number of spatial grid pts
% dt is the time step
%tic, u = spin2(S,200,2,'plot','off');
tic, u = spin2(S,200,2);
%
plot(u{2}), view(0,90), axis equal, axis off
time_in_seconds = toc
