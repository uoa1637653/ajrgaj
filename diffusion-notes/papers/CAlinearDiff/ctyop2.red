on div; off allfac; on revpri;
factor hh,uu,cc,gamma;
depend xi,x;
depend uu,t;
let df(~u,x) => df(u,xi)/hh;
let df(~u,t) => df(u,uu)*g;
let int(~u,x) => int(u,xi)*hh;
operator xxj;
let xxj(~u) => sub(xi=1,u);
operator xxjm1;
let xxjm1(~u) => sub(xi=0,u);
operator iij;
let iij(~u) => xxj(int(u,x))-xxjm1(int(u,x));

u0 := (xi+(1-xi)*em)*uu;
v0 := xi+ep*(1-xi);

let {gamma^3 => 0, cc^3 => 0 };
u := u0;
inner_rhs := gamma*xxj(delta^2*u)/hh-cc*(ep-1)*iij(u)/hh;
g := inner_rhs / iij(df(u,uu)*v0);
pde_rhs := df(u,t) + cc*df(u,x);
u := int(int(pde_rhs,x),x);
u := u - xxjm1(u-u0);
u := u - xi*xxj(u-u0);

let { ep => 1+mu*delta+delta^2/2,
      em => 1-mu*delta+delta^2/2,
	  mu^2 => 1+delta^2/4 };

%let { (sigma-sigma^-1)/2 => mu*delta,
%       sigma+sigma^-1-2 => delta^2,
%	   -1+sigma^-1+mu*delta => delta^2/2,
%	   1+delta^2/6 => ss^-1 };
