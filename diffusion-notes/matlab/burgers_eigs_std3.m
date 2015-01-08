function [vec,lam]=burgers_eigs_std3
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
% Find eigenvalues/vectors by perturbing u0 at t=0:
[vec,lam]=calc_eigs(u0_saw(), @burgers_dudt_std);
%----------------------------------------------------------------
% Plot eigenvectors versus domain, grouped by wavenumbers.
plot_eigs(vec,lam);
