function u=u0_zero
%----------------------------------------------------------------
% GAJ 07/01/2015
% Computes a zero displacement u(x,0) for the 
% discretised domain initialised via init_domain().
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
global L
u=zeros(L,1);
