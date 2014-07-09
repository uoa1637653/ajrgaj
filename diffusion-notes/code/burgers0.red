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

operator p, m, d2, md;
linear p, m, d2, md;
let p(m(uu,j),j) => uu;
let p(uu,j)^~n => p(uu^n,j);
let m(uu,j)^~n => m(uu^n,j);
let p(m(uu^~n,j),j) => uu^n;
let p(m(uu,j)*uu,j) => uu*p(uu,j);

let df(uu,t) => g;
let df(p(uu,j),t) => p(g,j);
let df(m(uu,j),t) => m(g,j);

u0 := xi*uu + (1-xi)*m(uu,j);
%v0 := xi + p(1-xi,j);


u := u0;
g := 0;

ux := df(u,xi)/hh$
amp := sub(xi=1,u) - uu;
cty := sub(xi=0,p(u,j)) - sub(xi=1,u);
jmp := sub(xi=0,p(ux,j)) - sub(xi=1,ux)
       - (1-gamma)*sub(xi=1,d2(u,j))/hh;
pde := -df(u,t) + df(ux,xi)/hh - epsilon*u*ux;
% Compute solvability function, int(v0*pde,x,-infty,infty):
% Use temporary variables to avoid weird error in integration.
tmp_j := pde*xi$ 
tmp_jp1 := pde - tmp_j$
slv := gamma*d2(uu,j)/hh + (int(tmp_j,xi,0,1) + p(int(tmp_jp1,xi,0,1),j))*hh;
gn := ss*slv/hh;

% Tidy up.
let p(~z,j) => 1 + md(z,j) + d2(z,j)/2;
let m(~z,j) => 1 - md(z,j) + d2(z,j)/2;
%let p(uu^~n,j) => 1 + md(uu^n,j) + d2(uu^n,j)/2;
%let m(uu^~n,j) => 1 - md(uu^n,j) + d2(uu^n,j)/2;

end;