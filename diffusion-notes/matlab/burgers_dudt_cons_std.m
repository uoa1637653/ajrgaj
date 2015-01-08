function u_t=burgers_dudt_cons_std(t,u)
%----------------------------------------------------------------
% GAJ 07/01/2015
% Computes the temporal derivative (RHS) of
% Burgers' eq. in conservative form with unit viscosity, 
% u_t=u_xx-(u^2)_x/2.
% Uses standard finite difference approximations for u_x and u_xx.
%
% Assumes pre-initialisation of the discretised domain and 
% spatio-differential operators via init_domain().
%----------------------------------------------------------------
% Temporal derivative of u(x,t):
%display(t)
%display(u)
global H D2 MD
u_t=D2*u/H^2-MD*(u.^2)/(2*H);
%display(u_t)
%----------------------------------------------------------------
