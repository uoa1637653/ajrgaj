function u_t=burgers_dudt_holistic(t,u)
%----------------------------------------------------------------
% GAJ 07/01/2015
% Computes the temporal derivative (RHS) of
% Burgers' eq. with unit viscosity, 
% u_t=u_xx-uu_x, using only first-order terms of the holistic
% approximation.
%
% Assumes pre-initialisation of the discretised domain and 
% spatio-differential operators via init_domain().
%----------------------------------------------------------------
% Temporal derivative of u(x,t):
global H D2 MD S
u_t=S*D2*u/H^2-S*(u.*(MD*u)+MD*(u.^2))/(3*H);
%----------------------------------------------------------------
