on div; off allfac; on revpri;
factor hh,uu,cc,gamma;

% Interpret x in terms of xi_j...
depend xi, x;
let { df(~u,x) => df(u,xi)/hh,
	  int(~u,x) => int(u,xi)*hh };
% and 1-xi_j (to match hand-integraton):
depend omxi, x;
let { df(omxi,xi) => -1,
	  int(omxi,xi) => -omxi^2/2,
	  int(omxi^~~p,xi) => -omxi^(p+1)/(p+1),
	  int(xi*omxi,xi) => -xi*omxi^2/2 - omxi^3/6,
	  int(xi*omxi^~~p,xi) => -xi*omxi^(p+1)/(p+1) - omxi^(p+2)/((p+1)*(p+2)),
	  int(xi^~q*omxi,xi) => -xi^q*imxi^2/2 + q*int(xi^(q-1)*omxi^2,xi)/2,
	  int(xi^~q*omxi^~~p,xi) => -xi^q*omxi^(p+1)/(p+1) + q*int(xi^(q-1)*omxi^(p+1),xi)/(p+1)
	};

% Specify initial approximation (em := sigma^-1)...
u0 := (xi + omxi*em)*uu;
% and its dual (ep := sigma):
v0 := xi + ep*omxi;
let { ep*em => 1 };

% Define end-points of jth interval...
operator xxj;
let xxj(~u) => sub({xi=1,omxi=0}, u);
operator xxjm1;
let xxjm1(~u) => sub({xi=0,omxi=1}, u);
% and the integral across the interval:
operator iij;
let iij(~u) => xxj(int(u,x)) - xxjm1(int(u,x));

% Define manifold operator...
operator omm;
let omm(~u,~~g) => df(u,x,x) - df(u,uu)*g - cc*df(u,x);
% and corresponding solvability operator:
operator oss;
let oss(~u,~~g) => gamma*xxj(delta^2*u)/hh - iij(df(u,uu)*v0)*g - cc*(ep-1)*iij(u)/hh;

let {gamma^3 => 0, cc^3 => 0 };
% Start with initial approximation:
u := u0;
g := 0;

% Derive next approximation:
sres := oss(u,g+g1);
g1 := -sub(g1=0, sres) / df(sres, g1);
g := g + g1;
mres := -omm(u,g);
u1 := int(int(mres,x),x);
u1 := u1 - xxjm1(u1);
u1 := u1 - xi*xxj(u1);
u := u + u1;

% Derive next approximation:
sres := oss(u,g+g2);
let {gamma*g2 => 0, cc*g2 => 0 };
g2 := -sub(g2=0, sres) / df(sres, g2);
g := g + g2;
mres := -omm(u,g);
%u2 := int(int(mres,x),x);
%u2 := u2 - xxjm1(u2);
%u2 := u2 - xi*xxj(u2);
%u := u + u2;

let { ep => 1+mu*delta+delta^2/2,
      em => 1-mu*delta+delta^2/2,
	  mu^2 => 1+delta^2/4,
	   1+delta^2/6 => ss^-1 };
%let { (sigma-sigma^-1)/2 => mu*delta,
%       sigma+sigma^-1-2 => delta^2,
%	   -1+sigma^-1+mu*delta => delta^2/2,
%	   1+delta^2/6 => ss^-1 };
