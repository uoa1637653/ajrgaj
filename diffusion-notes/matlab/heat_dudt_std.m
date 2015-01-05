function u_t=heat_dudt_std(t,u)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Computes the temporal derivative (RHS) of
% the heat eq., u_t=u_xx.
% Uses the standard finite difference approximation for u_xx.
%
% Assumes pre-initialisation of the discretised domain and 
% spatio-differential operators via init_domain().
%----------------------------------------------------------------
% Temporal derivative of u(x,t):
%display(t)
%display(u)
global H D2
u_t=D2*u/H^2;
%display(u_t)
%----------------------------------------------------------------
