function [vec,lam]=burgers_eigs_std_zero(period)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Checks stability of
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Uses standard finite difference approximations for u_x and u_xx.
% Solves on a periodic, discretised domain.
%----------------------------------------------------------------
% Set up domain and operators:
init_domain(period);
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing equilibrium at t=0:
[vec,lam]=calc_eigs(u0_zero(), @burgers_dudt_std);
%----------------------------------------------------------------
% Plot eigenvectors versus domain, grouped by wavenumbers.
plot_eigs(vec,lam);
