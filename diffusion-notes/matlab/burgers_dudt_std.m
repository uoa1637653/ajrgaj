function u_t=burgers_dudt_std(t,u)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Computes the temporal derivative (RHS) of
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Uses standard finite difference approximations for u_x and u_xx.
%
% Assumes pre-initialisation of the discretised domain and 
% spatio-differential operators via init_domain().
%----------------------------------------------------------------
% Temporal derivative of u(x,t):
%display(t)
%display(u)
global H D2 MD
u_t=D2*u/H^2-u.*(MD*u/H);
%display(u_t)
%----------------------------------------------------------------
