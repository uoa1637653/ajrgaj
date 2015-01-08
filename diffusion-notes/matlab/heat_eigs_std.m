function [vec,lam]=heat_eigs_std
%----------------------------------------------------------------
% GAJ 05/01/2015
% Checks the stability of the heat eq., u_t=u_xx,
% using the standard finite difference approximation for u_xx.
% Solves on a [0,1] periodic, discretised domain.
%----------------------------------------------------------------
% Set up domain and operators:
init_domain;
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing equilibrium at t=0:
[vec,lam]=calc_eigs(u0_zero(), @heat_dudt_std);
%----------------------------------------------------------------
% Plot eigenvectors versus domain, grouped by wavenumbers.
plot_eigs(vec,lam);
