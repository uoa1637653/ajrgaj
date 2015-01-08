function [vec,lam]=calc_eigs(u0,dudt)
%----------------------------------------------------------------
% GAJ 07/01/2015
% Checks the stability of an initial solution to a dynamical eq.
% on a [0,1] periodic, discretised domain (pre-initialised
% by init_domain()).
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing u0 at t=0:
small=1e-6;
fun=@(v) (dudt(0,u0+small*v)-dudt(0,u0))/small;
global L;
% Note: Ignore error generated (randomly) by complex values.
[vec,d]=eigs(fun,L,L);
lam=diag(d);
