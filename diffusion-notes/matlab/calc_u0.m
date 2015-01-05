function u=calc_u0
%----------------------------------------------------------------
% GAJ 05/01/2015
% Computes an initial displacement u(x,0) for the discretised
% domain initialised via init_domain().
% The global amplitude A needs to be pre-specified.
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
global A x j
u=A*sin(2*pi*x(j+1));
