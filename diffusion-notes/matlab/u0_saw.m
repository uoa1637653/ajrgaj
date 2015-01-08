function u=u0_saw
%----------------------------------------------------------------
% GAJ 05/01/2015
% Computes an initial sawtooth displacement u(x,0) for the 
% discretised domain initialised via init_domain().
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
global L
u=zeros(L,1);
u(1:2:L)=1;
u(2:2:L)=-1;
if mod(L,2)==1, u(L)=0; end
