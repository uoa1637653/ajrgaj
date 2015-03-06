function init_domain(period)
%----------------------------------------------------------------
% GAJ 05/01/2015
% Initialises differential operators on a [0,2\pi] periodic domain,
% as LxL matrices.
% Discretises the domain into L+1 points, labelled j=0..L,
% where x_0=0, x_L=1 and u(x_0,t)=u(x_L,t) by definition. 
% Hence, only indices j=1..L are required.
%----------------------------------------------------------------
% Set up domain:
global W L x H j
W=2*pi;                  % Width
L=period;                % Periodicity
x=linspace(0,W,L+1)'; % Grid-points
H=x(2)-x(1);             % Grid-point spacing
j=1:L;                   % Grid-point indices, j=0 equiv. j=L
% Set up operator mu*delta (centred difference) as a matrix:
global MD
MD=zeros(L,L);
for k=1:L
  MD(k,mod(k,L)+1)=1/2;
  MD(k,mod(L+k-2,L)+1)=-1/2;
end
% Set up operator delta^2 (centred 2nd difference) as a matrix:
global D2
D2=-2*eye(L,L);
for k=1:L
  D2(k,mod(k,L)+1)=1;
  D2(k,mod(L+k-2,L)+1)=1;
end
% Set up holistic operator S=(1+delta^2/6)^(-1) as a matrix:
global Sinv II
II=eye(L,L);
Sinv=(II+D2/6);
global SD2
SD2=6*(II-Sinv);
