function [vec,lam]=heat_eigs_std2
%----------------------------------------------------------------
% GAJ 05/01/2015
% Checks the stability of the heat eq., u_t=u_xx,
% using the standard finite difference approximation for u_xx.
% Solves on a [0,1] periodic, discretised domain.
%----------------------------------------------------------------
% Set up domain and operators:
init_domain;
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing u0 at t=0:
u0=calc_u0();
small=1e-6;
fun=@(v) (heat_dudt_std(0,u0+small*v)-heat_dudt_std(0,u0))/small;
global L;
% Note: Ignore error generated (randomly) by complex values.
[vec,d]=eigs(fun,L,L);
lam=diag(d);
%----------------------------------------------------------------
% Plot eigenvectors versus domain, grouped by wavenumbers.
plot_eigs(vec,lam);
