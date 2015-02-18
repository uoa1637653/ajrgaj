on div; off allfac;
let cos(~x)^2=>1/2*(cos(2*x)+1);
let sin(~x)^2=>1/2*(1-cos(2*x));
let cos(~x)*cos(~~y)=>1/2*(cos(x+y)+cos(x-y));
let sin(~x)*sin(~~y)=>1/2*(cos(x-y)-cos(x+y));
let cos(~x)*sin(~~y)=>1/2*(sin(x+y)-sin(x-y));

% Solve u_t=u_xx-u.u_x using variational perturbations about u=0,
% where u~u^(1)+u^(2)+..., vec{a}_t=vec{g}~vec{g}^(2)+vec{g}^(3)+...
% Define M(u,g)=u_xx-u_vec{a}.vec{g}-u.u_x to solve M(u,g)=0:
%operator MM;
%depend MM,x,a,b;
%let MM(~u)=>df(u,x,x)-df(u,a)*g-df(u,b)*h-u*df(u,x);
procedure MM(u,g,h); df(u,x,x)-df(u,a)*g-df(u,b)*h-u*df(u,x);

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

procedure Msolve(uhat,ghat,hhat);
begin scalar nn,nn_s,nn_c,g_slv,h_slv,gn,hn,un;
  % Compute N^(n)=N(u^(n-1),vec{g}^(n-1)), where M(u,vec{g})=Lu-N(u,vec{g}), 
  % and enforce solvability conditions on vec{g}^(n):
%  nn:=-sub(g=ghat+g,h=hhat+h,MM(uhat));
  nn:=-MM(uhat,ghat+g,hhat+h);
  nn_s:=int(nn*sin(x),x);
  nn_c:=int(nn*cos(x),x);
  g_slv:=sub(x=2*pi,nn_s)-sub(x=0,nn_s);
  h_slv:=sub(x=2*pi,nn_c)-sub(x=0,nn_c);
  gn:=-sub(g=0,g_slv)/(sub(g=1,g_slv)-sub(g=0,g_slv));
  hn:=-sub(h=0,h_slv)/(sub(h=1,h_slv)-sub(h=0,h_slv));
  % Solve Lu^(n)=N^(n) for u^(n):
  un:=Lsolve(sub(g=gn,h=hn,nn));
  return {un,gn,hn};
end;

% Choose u^(1)(x,t)=\vec{a}(t).[sin,cos](x),
% such that Lu^(1)=0 where Lu=u+u_xx.
u1:=a*sin(x)+b*cos(x);
uhat:=u1;
ghat:=0;
hhat := 0;

let { a^5=>0, b^5=>0, aa^9=>0 };
maxN:=6;
maxK:=2^maxN;
for i:=1:maxN do
begin
  alist:=Msolve(uhat,ghat,hhat);
  un:=first(alist);
  gn:=second(alist);
  hn:=third(alist);
  uhat:=uhat+un;
  ghat:=ghat+gn;
  hhat:=hhat+hn;
  write aadot:=((a*ghat+b*hhat)/aa where a^2+b^2=>aa^2);
end;

end;
