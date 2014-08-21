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
depend g, j;

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

%% Simplifications
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

% Temporal composition
let df(uu,t) => g;
operator !~f;
let df(~f(~~z,j),t) => f(df(z,t),j);

%on list;

% Initiate approximations:
u0 := xi*uu + (1-xi)*m(uu,j);
u := u0;
g := 0;

% Check internal boundary conditions:
amp := sub(xi=1,u) - uu;                  % u|X_j = U_j
cty := sub(xi=0,p(u,j)) - sub(xi=1,u);    % [u]_j = 0
ux := df(u,xi)/hh$
jmp := sub(xi=0,p(ux,j)) - sub(xi=1,ux)
       - (1-gamma)*sub(xi=1,d2(u,j))/hh;  % [u']_j = (1-gamma)/H*delta^2 U_j
pde := -df(u,t) + df(ux,xi)/hh - epsilon*u*ux;

% Satisfy solvability condition, <v0,pde> = 0, for g := g + g_n, 
% where v0 := xi + p(1-xi,j):
tmp_j := pde*xi$         % Use temporary variables to avoid
tmp_jp1 := pde - tmp_j$  % weird error in integration.
slv := (int(tmp_j,xi,0,1) + p(int(tmp_jp1,xi,0,1),j))*hh;
% Recover term <v0,d^2u/dx^2>, lost due to u0 being piecewise linear.
slv := slv + jmp;
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

% Check internal boundary conditions:
let gamma^2 => 0, epsilon^2 => 0;
amp := sub(xi=1,u) - uu;                  % u|X_j = U_j
cty := sub(xi=0,p(u,j)) - sub(xi=1,u);    % [u]_j = 0
ux := df(u,xi)/hh$
jmp := sub(xi=0,p(ux,j)) - sub(xi=1,ux)
       - (1-gamma)*sub(xi=1,d2(u,j))/hh;  % [u']_j = (1-gamma)/H*delta^2 U_j
pde := -df(u,t) + df(ux,xi)/hh - epsilon*u*ux;

% Satisfy solvability condition, <v0,pde> = 0, for g := g + g_n, 
% where v0 := xi + p(1-xi,j):
tmp_j := pde*xi$         % Use temporary variables to avoid
tmp_jp1 := pde - tmp_j$  % weird error in integration.
slv := (int(tmp_j,xi,0,1) + p(int(tmp_jp1,xi,0,1),j))*hh;
% Recover term <v0,d^2u/dx^2>, lost due to u0 being piecewise linear.
slv := slv + jmp;
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

% Check internal boundary conditions:
let gamma^2 => 0, epsilon^2 => 0;
amp := sub(xi=1,u) - uu;                  % u|X_j = U_j
cty := sub(xi=0,p(u,j)) - sub(xi=1,u);    % [u]_j = 0
ux := df(u,xi)/hh$
jmp := sub(xi=0,p(ux,j)) - sub(xi=1,ux)
       - (1-gamma)*sub(xi=1,d2(u,j))/hh;  % [u']_j = (1-gamma)/H*delta^2 U_j
pde := -df(u,t) + df(ux,xi)/hh - epsilon*u*ux;

end;