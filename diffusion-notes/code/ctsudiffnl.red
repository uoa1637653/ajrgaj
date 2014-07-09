%Comment derive slow manifold model for the dynamics of
%nonlinear diffusion on interval [-1,1].  Here ensure
%continuity at mid-point but form basis with discontinuous
%derivative at the mid-point. For simplicity use the int()
%function.  Parameter dd does not make any difference in this
%problem, so far.  But the derivative coupling is not what we
%are using now?  Tony Roberts, Mar 2014;
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
for iter:=1:9 do begin
write
	resd:=df(u,t)-u*df(u,x,2);
	ux:=df(u,x)$
write
    resb := sub(x=0,(ux where xpos)-(ux where xneg))
        -(1-gamma)*sub(x=0,-2*u);	
	g:=g+(gd:=-int(resd,x,-1,1)*dd+resb*a);
	write
	u:=u+int(int(resd*sign(x)/x+gd,x),x)/a$
	u:=u-x*sign(x)*sub(x=+1,u);
	
	if {resd,resb}={0,0} then write iter:=1000000+iter;
end;

resa:=a-sub(x=0,u);
resc:=sub(x=0,(u where xpos)-(u where xneg));
resl:=sub(x=-1,u);
resr:=sub(x=+1,u);


end;