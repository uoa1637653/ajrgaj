on div;
let cos(~x)*cos(~~y)=>1/2*(cos(x+y)+cos(x-y));
let sin(~x)*sin(~~y)=>1/2*(cos(x-y)-cos(x+y));
let cos(~x)*sin(~~y)=>1/2*(sin(x+y)-sin(x-y));

% Solve u_t=u_xx-u.u_x using variational perturbations...
% => Lu=N(u) where Lu=u_t-u_xx, N(u)=-u_vec{a}.vec{g}-u.u_x.
% Let u~u^(1)+u^(2)+..., vec{a}_t=vec{g}~vec{g}^(1)+vec{g}^(2)+...

% Follow mode 1 satisfying Lu^(1)=0:
u1:=a*exp(-t)*sin(x)+b*exp(-t)*cos(x);

n2:=-df(u1,a)*g1-df(u1,b)*h1-u1*df(u1,x);
n2s:=int(n2*exp(t)*sin(x),x);
n2c:=int(n2*exp(t)*cos(x),x);
sub(x=2*pi,n2s)-sub(x=0,n2s);
sub(x=2*pi,n2c)-sub(x=0,n2c);

% Solve manually:
g1:=0; h1:=0;
u2:=exp(-2*t)*(-1/4*(a^2-b^2)*sin(2*x)-1/2*a*b*cos(2*x));
res2:=df(u2,t)-df(u2,x,x)-n2;

n3:=-df(u2,a)*g1-df(u2,b)*h1-df(u1,a)*g2-df(u1,b)*h2-df(uu1*df(u1,x);
n3s:=int(n3*exp(t)*sin(x),x);
n3c:=int(n3*exp(t)*cos(x),x);
sub(x=2*pi,n3s)-sub(x=0,n3s);
sub(x=2*pi,n3c)-sub(x=0,n3c);

% Solve manually:
g2:=-1/8*a*(a^2+b^2)*exp(-2*t);
h2:=-1/8*b*(a^2+b^2)*exp(-2*t);
