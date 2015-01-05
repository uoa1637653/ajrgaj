function [v,lam]=burgers_eigs_std
%----------------------------------------------------------------
% GAJ 05/01/2015
% Checks stability of
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Uses standard finite difference approximations for u_x and u_xx.
% Solves on a [0,1] periodic, discretised domain.
%----------------------------------------------------------------
% Set up domain and operators:
init_domain;
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing equilibrium at t=0:
global L;
small=1e-6;
fun=@(v) burgers_dudt_std(0,small*v)/small;
% Note: Ignore error generated (randomly) by complex values.
[vec,d]=eigs(fun,L,L);
lam=diag(d);
%----------------------------------------------------------------
% Plot eigenvectors versus domain, grouped by wavenumbers.
plot_eigs(vec,lam);
