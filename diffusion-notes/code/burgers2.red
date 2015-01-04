%% Key: 
%% hh := H = X_j - X_{j-1}
%% xi := xi_j = (x - X_{j-1}) / H
%% uu = U_j
%% p := E^+ = sigma, right-shift
%% m := E^- = sigma^{-1}, left-shift
%% d2 := delta^2 = p + m - 2
%% ss := S = (1 + delta^2/6)^{-1}
%% md := mu*delta = (p - m)/2

on div; off allfac; on revpri;
factor hh, epsilon, gamma;

depend xi, j;
depend uu, j, t;
depend gg, j;

operator p, m, d2, md, ss;
linear p, m, d2, md, ss;

%% Expansions
let p(~z,j) => z + md(z,j) + d2(z,j)/2,
    m(~z,j) => z - md(z,j) + d2(z,j)/2;

%% Independence
let md(1,j) => 0,
    d2(1,j) => 0,
	ss(1,j) => 1;

%% Canonical orderings
let md(ss(~z,j),j) => ss(md(z,j),j);
let d2(ss(~z,j),j) => ss(d2(z,j),j);
let md(d2(~z,j),j) => d2(md(z,j),j);

%% Invariants:
% Next follows from direct expansion, and is required for p(m(z)) = m(p(z)) = z:
let md(md(~z,j),j) => d2(z,j) + d2(d2(z,j),j)/4;
% Next follows from definition of S:
let ss(d2(~z,j),j) => 6*(z-ss(z,j));
% Next follows from expanding mu*delta z^2:
let d2(~y,j)*md(~z,j) => md(z^2,j) - 2*z*md(z,j) when y=z;
% Next two follow from p(y*m(z)) = p(y)*z and m(y*p(z)) = m(y)*z:
let md(~y*md(~z,j),j) => 1/2*(1/2*d2(y*d2(z,j),j) + y*d2(z,j) - z*d2(y,j) + d2(y*z,j));
let md(~y*d2(~z,j),j) => d2(y*md(z,j),j) + 2*(y*md(z,j) + z*md(y,j) - md(y*z,j));
% Next follows from either m(z)^2 = m(z^2) or p(z)^2 = p(z^2):
let md(~z,j)^2 => 1/2*d2(z^2,j) - 1/4*d2(z,j)^2 - z*d2(z,j);

% Temporo-spatial composition:
let df(uu,t) => gg;
operator !~f;
let df(~f(~~z,j),t) => f(df(z,t),j);
let df(~z,x) => df(z,xi)/hh;
let df(~z,x,2) => df(z,xi,2)/hh^2;

% Initiate approximations:
u0 := xi*uu + (1-xi)*m(uu,j);
u := u0;
g := 0;

let gamma^2 => 0, epsilon^2 => 0;
for iter := 1:2 do begin

% Check internal boundary conditions:
amp := sub(xi=1,u) - uu;                  % u|X_j = U_j
cty := sub(xi=0,p(u,j)) - sub(xi=1,u);    % [u]_j = 0
ux := df(u,x)$
jmp := sub(xi=0,p(ux,j)) - sub(xi=1,ux)
       - (1-gamma)*sub(xi=1,d2(u,j))/hh;  % [u']_j = (1-gamma)/H*delta^2 U_j
pde := -sub(gg=g,df(u,t)) + df(ux,x) - epsilon*u*ux;

% Satisfy solvability condition, <v0,pde> = 0, where v0 := xi + p(1-xi,j), to obtain g_n;
% ensure internal boundary conditions are met.
% (Note: Use temporary variables to avoid weird error in integration):
pde_xi := pde*xi$         
pde_1mxi := (1-xi)*pde$
slv := (int(pde_xi,xi,0,1) + p(int(pde_1mxi,xi,0,1),j))*hh + jmp;
% Update g from error in solvability:
gn := ss(slv,j)/hh;
% Update u by solving pde = 0 for u := u + u_n:
tn := xi*gn + (1-xi)*m(gn,j) - pde;
un := hh^2*int(int(tn,xi),xi);
% Impose integration constants to satsify u_n|X_{j-1} = 0, u_n|X_j = 0:
un := un - sub(xi=1,un)*xi;
% Update iteration:
u := u + un;
g := g + gn;

end;

% Check internal boundary conditions:
amp := sub(xi=1,u) - uu;                  % u|X_j = U_j
cty := sub(xi=0,p(u,j)) - sub(xi=1,u);    % [u]_j = 0
ux := df(u,x)$
jmp := sub(xi=0,p(ux,j)) - sub(xi=1,ux)
       - (1-gamma)*sub(xi=1,d2(u,j))/hh;  % [u']_j = (1-gamma)/H*delta^2 U_j
pde := -sub(gg=g,df(u,t)) + df(ux,x) - epsilon*u*ux;

% Debugging...
% Start with previous terms (which have been checked via computer and by hand):
u00 := coeffn(coeffn(u,gamma,0),epsilon,0)$
u10 := coeffn(coeffn(u,gamma,1),epsilon,0)$
g10 := coeffn(coeffn(g,gamma,1),epsilon,0)$
u01 := coeffn(coeffn(u,gamma,0),epsilon,1)$
g01 := coeffn(coeffn(g,gamma,0),epsilon,1)$
% Computer help with 'hand' solution:
depend g11, uu;
u11dd_g := sub(gg=g11,df(u00,t))+sub(gg=g10,df(u01,t))+sub(gg=g01,df(u10,t))+u00*df(u10,x)+u10*df(u00,x)$
depend c11, uu;
u11d_cg := hh*int(u11dd_g,xi)+c11$
u11_cg := hh*int(u11d_cg,xi)$
% Compute c11 that makes u11=0 at xi=1:
c11_g := c11-sub(xi=1,u11_cg)/hh;
u11d_g := sub(c11=c11_g,u11d_cg);
u11_g := sub(c11=c11_g,u11_cg);
% Check identity:
df(u11_g,x)-u11d_g;
% Compute g11 that makes [u']=0:
jmp11_g := p(sub(xi=0,u11d_g),j)-sub(xi=1,u11d_g)$
g11h := ss(sub(g11=0,jmp11_g),j)/hh$
% Check identity:
sub(g11=g11h,jmp11_g);
u11h := sub(g11=g11h,u11_g)$
u11dh := sub(g11=g11h,u11d_g)$
df(u11h,x)-u11dh;
% Check hand solution matches computer:
g11c := coeffn(coeffn(g,gamma,1),epsilon,1)$
g11h-g11c;
u11c := coeffn(coeffn(u,gamma,1),epsilon,1)$
u11h-u11c;

% Determine why jump condition is seemingly not satisfied...
% Recheck jump:
jmp11h := p(sub(xi=0,u11dh),j)-sub(xi=1,u11dh);
sub(xi=1,u11dh)-sub(g11=g11h,xi=1,u11d_g);
sub(xi=0,u11dh)-sub(g11=g11h,xi=0,u11d_g);
p(sub(xi=0,u11dh),j)-p(sub(g11=g11h,xi=0,u11d_g),j);
p(sub(g11=g11h,xi=0,u11d_g),j)-sub(g11=g11h,p(sub(xi=0,u11d_g),j));
u11d_g_only := u11d_g - sub(g11=0,u11d_g)$
p(sub(g11=g11h,xi=0,u11d_g_only),j)-sub(g11=g11h,p(sub(xi=0,u11d_g_only),j));
% Noted that sub(xi=0,u11d_g_only) has only 3 terms: g11, md(g11) and d2(g11):
p(g11h,j)-sub(g11=g11h,p(g11,j));
p(d2(g11h,j),j)-sub(g11=g11h,p(d2(g11,j),j));
p(md(g11h,j),j)-sub(g11=g11h,p(md(g11,j),j)); %% Non-zero!!!
% Noted that only p(md()) difference is non-zero:
md(g11h,j)-sub(g11=g11h,md(g11,j));
% Test g11h term by term:
depend z,uu;
n:= 1$ zh:=uu^2$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:= 2$ zh:=uu*ss(uu,j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:= 3$ zh:=ss(uu*ss(uu,j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:= 4$ zh:=ss(uu*ss(md(uu,j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:= 5$ zh:=ss(uu*md(uu,j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:= 6$ zh:=ss(md(uu,j)*ss(uu,j),j)$ diff6:=p(md(zh,j),j)-sub(z=zh,p(md(z,j),j)); %% Non-zero!!!
n:= 7$ zh:=ss(md(uu,j)*ss(md(uu,j),j),j)$ diff7:=p(md(zh,j),j)-sub(z=zh,p(md(z,j),j)); %% Non-zero!!
n:= 8$ zh:=ss(uu*d2(uu,j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:= 9$ zh:=ss(d2(uu,j)*ss(md(uu,j),j),j)$ diff9:=p(md(zh,j),j)-sub(z=zh,p(md(z,j),j)); %% Non-zero!!
n:=10$ zh:=ss(ss(uu*ss(uu,j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=11$ zh:=ss(ss(uu*ss(md(uu,j),j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=12$ zh:=ss(ss(uu*md(uu,j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=13$ zh:=ss(ss(ss(uu,j)*md(uu,j),j),j)$ diff13:=p(md(zh,j),j)-sub(z=zh,p(md(z,j),j)); %% Non-zero!!
n:=14$ zh:=ss(ss(ss(uu,j)*d2(uu,j),j),j)$ diff14:=p(md(zh,j),j)-sub(z=zh,p(md(z,j),j)); %% Non-zero!!
n:=15$ zh:=ss(ss(ss(uu*md(uu,j),j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=16$ zh:=ss(ss(ss(md(uu^2,j),j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=17$ zh:=ss(ss(md(uu*ss(uu,j),j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=18$ zh:=ss(ss(md(uu*ss(md(uu,j),j),j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=19$ zh:=ss(ss(md(uu^2,j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=20$ zh:=ss(ss(uu^2,j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=21$ zh:=uu*ss(md(uu,j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=22$ zh:=ss(md(uu*ss(uu,j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=23$ zh:=ss(md(uu*ss(md(uu,j),j),j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=24$ zh:=ss(md(uu^2,j),j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));
n:=25$ zh:=ss(uu^2,j)$ p(md(zh,j),j)-sub(z=zh,p(md(z,j),j));

% Accumulate terms across all potential invariants:
terms := diff6+diff7+diff9+diff13+diff14;

% Apply invariants from ssrelations.red:
let ss(md(ss(uu,j)*uu,j),j) =>
  ss(ss(md(uu,j),j)*uu,j)
  - 1/2*ss(md(uu,j)*ss(uu,j),j)
  - 3/2*ss(md(uu,j),j)*uu
  + 3/2*ss(md(uu^2,j),j);

let ss(md(uu,j)*ss(md(uu,j),j),j) =>
  18*uu**2 
  - 9*ss(uu,j)*uu 
  + 6*ss(ss(uu,j)*uu,j) 
  - 3/2*ss(d2(uu,j)*uu,j) 
  - 2*ss(md(ss(md(uu,j),j)*uu,j),j) 
  - 15*ss(uu**2,j);

let ss(d2(uu,j)*ss(md(uu,j),j),j) =>
  - 6*ss(md(uu,j)*uu,j) 
  + 3*ss(md(uu,j)*ss(uu,j),j) 
  - 3*ss(md(uu,j),j)*uu 
  + 3*ss(md(uu**2,j),j);

let ss(md(ss(md(uu,j),j)*uu,j),j) =>
  - 6*ss(uu,j)*uu
  + 3*ss(ss(uu,j)*uu,j)
  - 1/2*ss(d2(uu,j)*ss(uu,j),j) 
  - 6*ss(uu**2,j) 
  + 9*uu**2;

end;