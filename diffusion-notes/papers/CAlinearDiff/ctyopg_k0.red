on div; off allfac; on revpri;
factor hh,uu,cc,d;
ep:=1+mu*delta+delta^2/2;
em:=1-mu*delta+delta^2/2;
let { mu^2=>1+delta^2/4
    , ss*delta^2=>6-6*ss };
depend xi,x; 
let df(~a,x)=>df(a,xi)/hh;
operator linv; linear linv;
let { linv(xi^~~p,xi)=>(xi^(p+2)-xi)/(p+1)/(p+2)
    , linv(1,xi)=>(xi^2-xi)/2 };
depend uu,t; 
let df(uu,t)=>g;
g:=0;
u:=xi*uu+(1-xi)*em*uu;
let { gamma^7=>0, cc^3=>0 };
for it:=1:99 do begin 
    pde:= -df(u,t)+df(u,x,x)-cc*df(u,x);
    amp:=sub(xi=1,u)-uu;
    cty:=sub(xi=0,ep*u)-sub(xi=1,u);
    hux:=hh*df(u,x)$
    jmp:=-sub(xi=0,ep*hux)+sub(xi=1,hux)
        +(1-gamma)*sub(xi=1,ep*u-2*u+em*u);
    write lengthres:=map(length(~a),{pde,amp,cty,jmp});
    g:=g+(gd:=-ss*jmp/hh^2);
    u:=u-linv(pde-(xi+(1-xi)*em)*gd,xi)*hh^2;
    showtime;
    if {pde,amp,cty,jmp}={0,0,0,0} then write it:=it+10000;
end;

let kappa^16 => 0;
%% Expansion of S\tilde{u}=\tilde{u} about kappa=0.
%% see taylor(3/(2+cos(kappa)),kappa,0,10); etc.
%% for higher terms, use taylor(3/(2+cos(kappa))-lower_terms, ...);
ssk := 1+kappa^2/6+kappa^4/72+kappa^6/2160-17*kappa^8/362880-11*kappa^10/1209600
       -1079*kappa^12/1437004800-18197*kappa^14/784604620800;
%% Expansion of sin(kappa) about kappa=0 where \mu\delta\tilde{u}=is\tilde{u} 
sk := kappa-kappa^3/6+kappa^5/120-kappa^7/5040+kappa^9/362880-kappa^11/39916800+kappa^13/6227020800-kappa^15/1307674368000;
%% Partial sums
pg := (1+gamma+gamma^2+gamma^3+gamma^4+gamma^5+gamma^6)*g;
factor i;
%% Expand about kappa=0 with dimensionalisation adjustment.
pgk0 := sub(ss=ssk,delta=i,mu=sk,hh=1,uu=1,pg);

for n:=0:6 do begin;
	write n;
	write coeffn(pgk0,gamma,n);
end;
on rounded;
print_precision 4;
for n:=0:6 do begin;
	write n;
	write coeffn(pgk0,gamma,n);
end;

%% End script.
end;