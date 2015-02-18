on div;
let cos(~x)^2=>1/2*(cos(2*x)+1);
let sin(~x)^2=>1/2*(1-cos(2*x));
let cos(~x)*cos(~~y)=>1/2*(cos(x+y)+cos(x-y));
let sin(~x)*sin(~~y)=>1/2*(cos(x-y)-cos(x+y));
let cos(~x)*sin(~~y)=>1/2*(sin(x+y)-sin(x-y));

% Solve u_t=u_xx-u.u_x using variational perturbations about u=0,
% where u~u^(1)+u^(2)+..., vec{a}_t=vec{g}~vec{g}^(2)+vec{g}^(3)+...
% Step 1. Let u^(1)(x,t)=\vec{a}(t).[sin,cos](x).
u1:=a*sin(x)+b*cos(x);
% Step 2. Let Lu=u+u_xx so that Lu^(1)=0.
% Step 3. Derive Lu=N(u) where N(u)=epsilon.u+u_vec{a}.vec{g}+u.u_x.
n2:=epsilon*u1+df(u1,a)*g2+df(u1,b)*h2+u1*df(u1,x);
n2s:=int(n2*sin(x),x);
n2c:=int(n2*cos(x),x);
g2slv:=sub(x=2*pi,n2s)-sub(x=0,n2s);
h2slv:=sub(x=2*pi,n2c)-sub(x=0,n2c);
g2:=-sub(g2=0,g2slv)/(sub(g2=1,g2slv)-sub(g2=0,g2slv));
h2:=-sub(h2=0,h2slv)/(sub(h2=1,h2slv)-sub(h2=0,h2slv));

maxK:=10;
% Define solver for Lu=N,
% noting that Lsin(kx)=(1-k^2)sin(kx), etc.:
procedure Lsolve nn;
begin scalar un,ck,sk;
  un:=0;
  for k:=2:maxK do begin
    ck:=coeffn(nn,cos(k*x),1);
    sk:=coeffn(nn,sin(k*x),1);
	un:=un+(ck*cos(k*x)+sk*sin(k*x))/(1-k^2);
  end;
  return un;
end;

u2:=Lsolve(n2);
res2:=u2+df(u2,x,x)-n2;

% Repeat:
n3:=epsilon*u2+df(u1,a)*g3+df(u1,b)*h3+df(u2,a)*g2+df(u2,b)*h2+u2*df(u1,x)+u1*df(u2,x);
n3s:=int(n3*sin(x),x);
n3c:=int(n3*cos(x),x);
g3slv:=sub(x=2*pi,n3s)-sub(x=0,n3s);
h3slv:=sub(x=2*pi,n3c)-sub(x=0,n3c);
g3:=-sub(g3=0,g3slv)/(sub(g3=1,g3slv)-sub(g3=0,g3slv));
h3:=-sub(h3=0,h3slv)/(sub(h3=1,h3slv)-sub(h3=0,h3slv));
u3:=Lsolve(n3);
res3:=u3+df(u3,x,x)-n3;

% Repeat:
n4:=epsilon*u3+df(u1,a)*g4+df(u1,b)*h4+df(u2,a)*g3+df(u2,b)*h3+df(u3,a)*g2+df(u3,b)*h2+u3*df(u1,x)+u1*df(u3,x)+u2*df(u2,x);
n4s:=int(n4*sin(x),x);
n4c:=int(n4*cos(x),x);
g4slv:=sub(x=2*pi,n4s)-sub(x=0,n4s);
h4slv:=sub(x=2*pi,n4c)-sub(x=0,n4c);
g4:=-sub(g4=0,g4slv)/(sub(g4=1,g4slv)-sub(g4=0,g4slv));
h4:=-sub(h4=0,h4slv)/(sub(h4=1,h4slv)-sub(h4=0,h4slv));
u4:=Lsolve(n4);
res4:=u4+df(u4,x,x)-n4;
