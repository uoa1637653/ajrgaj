function [vec,lam]=heat_eigs_std_saw(period,amp)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Checks the stability of the heat eq., u_t=u_xx,
% using the standard finite difference approximation for u_xx.
% Solves on a periodic, discretised domain.
%----------------------------------------------------------------
% Set up domain and operators:
init_domain(period);
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing u0 at t=0:
[vec,lam]=calc_eigs(amp*u0_saw(), @heat_dudt_std);
%----------------------------------------------------------------
% Plot eigenvectors versus domain, grouped by wavenumbers.
plot_eigs(vec,lam);
