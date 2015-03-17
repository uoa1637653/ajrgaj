function u=u0_hat
%----------------------------------------------------------------
% GAJ 17/03/2015
% Computes an initial sinusoidal displacement u(x,0) for the 
% discretised domain initialised via init_domain().
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
global x L
u=zeros(L,1);
idx=1:floor(L/2);
u(idx)=x(idx+1)/pi;
idx=ceil(L/2):L;
u(idx)=2-x(idx+1)/pi;
