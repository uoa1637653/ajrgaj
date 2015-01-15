function [nums,amps,lams]=burgers_bifurc_fornberg_sin(period,amp_start,amp_end,amp_step)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Checks bifurcation stability of
% Burgers' eq. in Fornberg form with unit viscosity, u_t=u_xx-(uu_x+(u^2)x)/3.
% Uses standard finite difference approximations for u_x and u_xx.
% Solves on a periodic, discretised domain.
%----------------------------------------------------------------
% Set up domain and operators:
init_domain(period);
%----------------------------------------------------------------
% Find eigenvalues/vectors by perturbing u0 at t=0:
n=0; % Number of eigenvalues with imaginary part.
nums=[];
amps=[];
lams=[];
for A=amp_start:amp_step:amp_end
  [~,lam]=calc_eigs(A*u0_sin(), @burgers_dudt_fornberg);
  m=sum(abs(imag(lam))>0);
  if (m > n)
      %display(A);
      %display(lam);
      n = m;
      nums=[nums,n];
      amps=[amps,A];
      lams=[lams,lam];
  end
end
%----------------------------------------------------------------
