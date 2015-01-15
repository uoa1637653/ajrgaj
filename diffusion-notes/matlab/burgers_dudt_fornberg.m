function u_t=burgers_dudt_fornberg(t,u)
%----------------------------------------------------------------
% GAJ 07/01/2015
% Computes the temporal derivative (RHS) of
% Burgers' eq. in Fornberg form with unit viscosity, 
% u_t=u_xx-(1/3)uu_x-(2/3)(u^2)_x/2.
% Uses standard finite difference approximations for u_x and u_xx.
%
% Assumes pre-initialisation of the discretised domain and 
% spatio-differential operators via init_domain().
%----------------------------------------------------------------
% Temporal derivative of u(x,t):
%display(t)
%display(u)
global H D2 MD
u_t=D2*u/H^2-(u.*(MD*u)+MD*(u.^2))/(3*H);
%display(u_t)
%----------------------------------------------------------------
