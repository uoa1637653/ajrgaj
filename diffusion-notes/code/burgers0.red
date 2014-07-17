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
factor hh, xi, epsilon, gamma, ss;

depend xi, j;
depend uu, j, t;
depend g, j;

operator p, m, d2, md, ss;
linear p, m, d2, md, ss;

% Independent
let p(1,j) => 1;

% Inverses
let p(m(~z,j),j) => z;
let m(p(~z,j),j) => z;
let p(~z,j)^~n => p(z^n,j);
let m(~z,j)^~n => m(z^n,j);
let p(m(~x,j)*~y,j) => x*p(y,j);

% Canonical ordering
let md(ss(~z,j),j) => ss(md(z,j),j);
let d2(ss(~z,j),j) => ss(d2(z,j),j);
let md(d2(~z,j),j) => d2(md(z,j),j);

% Simplification, using S delta^2 = 6(1-S):
let ss(d2(d2(~z,j),j),j) => 6*(d2(z,j)-ss(d2(z,j),j));
let ss(d2(md(~z,j),j),j) => 6*(md(z,j)-ss(md(z,j),j));
let ss(md(md(~z,j),j),j) => (3*d2(z,j)-ss(d2(z,j),j))/2;

% Temporal composition
let df(uu,t) => g;
operator !~f;
let df(~f(~~z,j),t) => f(df(z,t),j);

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
slv := slv + gamma*d2(uu,j)/hh;
% Tidy up:
slv := (slv where {
   p(~z,j) => 1 + md(z,j) + d2(z,j)/2,
   m(~z,j) => 1 - md(z,j) + d2(z,j)/2
});
% Update g from error in solvability:
gn := ss(slv,j)/hh;
% Update u by solving pde = 0 for u := u + u_n:
tn := xi*gn + (1-xi)*m(gn,j) - pde;
un := hh^2*int(int(tn,xi),xi);
% Impose integration constants to satsify u_n|X_{j-1} = 0, u_n|X_j = 0:
un := un - sub(xi=1,un)*xi;
% Tidy up:
un := (un where {
   p(~z,j) => 1 + md(z,j) + d2(z,j)/2,
   m(~z,j) => 1 - md(z,j) + d2(z,j)/2
});

% Tidy up.
%let ss(d2(~z,j),j) => 6*(1-ss(z,j));

end;