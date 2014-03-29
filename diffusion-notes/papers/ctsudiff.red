Comment derive slow manifold model for the dynamics of
linear diffusion on interval [-1,1].  Here ensure
continuity at mid-point but form basis with
discontinuous derivative at the mid-point. Find good
convergence in gamma, but not spectacular. Here the
update rule is not optimal as takes twice as many
iterations as should be necessary. For simplicity use
the int() function although gets very slow for large
problems.  Tony Roberts, June 2010;
on div; off allfac; on revpri;
factor gamma,a;
let { df(sign(x),x)=>0, sign(~x)^2=>1 };
xpos:={sign(x)=>+1}$
xneg:={sign(x)=>-1}$

depend a,t;
let df(a,t)=>g;
u:=a*(1-x*sign(x));
g:=0;

let gamma^9=>0;
for iter:=1:99 do begin
write
	resd:=df(u,t)-df(u,x,2);
	ux:=df(u,x)$
write
	resb:=sub(x=0,(ux where xpos)-(ux where xneg))*(1+gamma/4)
		-(1-gamma)*(sub(x=+1,ux)-sub(x=-1,ux));
	
	g:=g+(gd:=-int(resd,x,-1,1)/3+resb);
	u:=u+int(int(resd+gd*(1-x*sign(x)),x),x)$
	u:=u-x*sign(x)*sub(x=+1,u);
	
	if {resd,resb}={0,0} then write iter:=1000000+iter;
end;

resa:=a-sub(x=0,u);
resc:=sub(x=0,(u where xpos)-(u where xneg));
resl:=sub(x=-1,u);
resr:=sub(x=+1,u);


end;