function u=u0_sin
%----------------------------------------------------------------
% GAJ 05/01/2015
% Computes an initial sinusoidal displacement u(x,0) for the 
% discretised domain initialised via init_domain().
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
global W x j
u=sin(2*pi/W*x(j+1));
