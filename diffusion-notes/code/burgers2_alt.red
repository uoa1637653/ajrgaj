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

% Spatial consistency (functions shifted in j take same form):
%%let md(xi,j) => xi*md(1,j),
%%    md(~z*xi,j) => xi*md(z,j),
%%    md(xi^~n,j) => xi^n*md(1,j),
%%	md(~~z*xi^~n,j) => xi^n*md(z,j);
%%let d2(xi,j) => xi*d2(1,j),
%%    d2(~z*xi,j) => xi*d2(z,j),
%%    d2(xi^~n,j) => xi^n*d2(1,j),
%%	d2(~~z*xi^~n,j) => xi^n*d2(z,j);
	
%on list;

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
pde := -sub(gg=g,df(u,t)) + df(ux,x) - epsilon*df(u^2,x)/2;

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
pde := -sub(gg=g,df(u,t)) + df(ux,x) - epsilon*df(u^2,x)/2;

end;