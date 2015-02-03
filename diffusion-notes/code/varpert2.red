on div;
let cos(~x)^2=>1/2*(cos(2*x)+1);
let sin(~x)^2=>1/2*(1-cos(2*x));
let cos(~x)*cos(~~y)=>1/2*(cos(x+y)+cos(x-y));
let sin(~x)*sin(~~y)=>1/2*(cos(x-y)-cos(x+y));
let cos(~x)*sin(~~y)=>1/2*(sin(x+y)-sin(x-y));

% Solve u_t=u_xx-u.u_x using variational perturbations about u=0,
% where u~u^(1)+u^(2)+..., vec{a}_t=vec{g}~vec{g}^(1)+vec{g}^(2)+...
% Step 1. Let u^(1)(x,t)=\vec{a}(t).[sin,cos](x).
u1:=a*sin(x)+b*cos(x);
% Step 2. Let Lu=u-u_xx so that Lu^(1)=0.
% Step 3. Derive Lu=N(u) where N(u)=epsilon.u-u_vec{a}.vec{g}-u.u_x.
n2:=epsilon*u1-df(u1,a)*g1-df(u1,b)*h1-u1*df(u1,x);
n2s:=int(n2*sin(x),x);
n2c:=int(n2*cos(x),x);
g1slv:=sub(x=2*pi,n2s)-sub(x=0,n2s);
h1slv:=sub(x=2*pi,n2c)-sub(x=0,n2c);
g1:=-sub(g1=0,g1slv)/(sub(g1=1,g1slv)-sub(g1=0,g1slv));
h1:=-sub(h1=0,h1slv)/(sub(h1=1,h1slv)-sub(h1=0,h1slv));
% solve by hand:
u2:=-a*b/5*cos(2*x)+1/10*(b^2-a^2)*sin(2*x);
res2:=u2-df(u2,x,x)-n2;

% Repeat:
n3:=epsilon*u2-df(u1,a)*g2-df(u1,b)*h2-df(u2,a)*g1-df(u2,b)*h1-u2*df(u1,x)-u1*df(u2,x);
n3s:=int(n3*sin(x),x);
n3c:=int(n3*cos(x),x);
g2slv:=sub(x=2*pi,n3s)-sub(x=0,n3s);
h2slv:=sub(x=2*pi,n3c)-sub(x=0,n3c);
g2:=-sub(g2=0,g2slv)/(sub(g2=1,g2slv)-sub(g2=0,g2slv));
h2:=-sub(h2=0,h2slv)/(sub(h2=1,h2slv)-sub(h2=0,h2slv));
