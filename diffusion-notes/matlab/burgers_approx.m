function burgers_approx
% GAJ 24/12/2014
% Burgers' eq. with unit viscosity, u_t=u_xx-uu_x.
% Solve on a [0,1] periodic domain.
% Use finite difference approximations for u_x and u_xx.
%----------------------------------------------------------------
% Set up domain:
global L x H j te ye ie
L=19;                 % Periodicity
x=linspace(0,2*pi,L+1)'; % Grid-points
H=x(2)-x(1);          % Grid-point spacing
j=1:L;                % Grid-point indices, j=0 equiv. j=L
% Set up operator mu*delta (centred difference) as a matrix:
global md
md=zeros(L,L);
for k=1:L
  md(k,mod(k,L)+1)=1/2;
  md(k,mod(L+k-2,L)+1)=-1/2;
end
% Set up operator delta^2 (centred 2nd difference) as a matrix:
global d2
d2=-2*eye(L,L);
for k=1:L
  d2(k,mod(k,L)+1)=1;
  d2(k,mod(L+k-2,L)+1)=1;
end
% Set up operator S=(1+delta^2/6)^(-1) as a matrix:
global S
S=inv(eye(L,L)+d2/6);
%----------------------------------------------------------------
% Solve PDE from t=0 to t=T:
global A
A = 50;
T=2;
if 0, [t,u]=ode15s(@dudt,[0 T],u0());
else opts=odeset('Events',@largey);
[t,u,te,ye,ie]=ode15s(@dudt,[0 T],u0(),opts);
te=te
ye=ye
ie=ie
end
u=[u(:,L) u]; % Add j=0 point from j=L
surf(x,t,asinh(u));
v=sum(abs(diff(sign(diff(u')))))/4;
display(v);
%----------------------------------------------------------------
% Initialisation of u(x,t) at t=0:
function u=u0()
global A L x
u=A*sin(x(2:L+1));
%----------------------------------------------------------------
% Temporal derivative of u(x,t):
function u_t=dudt(t,u)
%display(t)
%display(u)
global H d2 md S
u_t=S*(d2*u/H^2-1/3*(u.*(md*u)+md*(u.*u))/H);
display(u_t)
%----------------------------------------------------------------
function [val,isfin,dirn]=largey(t,y)
isfin=1;
val=-max(abs(y))+1e3;
dirn=0;
