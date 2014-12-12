% try to find nontrivial relations between operators
% AJR Dec 2014

on div; off allfac; on revpri;
operator ss;
operator md;
operator d2;
operator uu; 
%jmp:=a0*uu+a1*ss(uu,j)+a2*d2(uu,j)+a3*md(uu,j)+a4*ss(d2(uu,j),j);
nas:=13;
jmp:=a0*ss(ss(md(uu,j),j)*uu,j)
    +a1*ss(md(uu,j)*ss(uu,j),j)
    +a2*ss(md(uu,j),j)*uu
    +a3*ss(md(ss(uu,j)*uu,j),j)
    +a4*ss(md(uu^2,j),j)
    +a5*uu^2
    +a6*ss(uu,j)*uu
    +a7*ss(ss(uu,j)*uu,j)
    +a8*ss(md(uu,j)*ss(md(uu,j),j),j)
    +a9*ss(d2(uu,j)*uu,j)
    +a10*ss(md(ss(md(uu,j),j)*uu,j),j)
    +a11*ss(uu^2,j)
    +a12*ss(md(uu,j)*uu,j)
    +a13*ss(d2(uu,j)*ss(md(uu,j),j),j);
jmpy:=(jmp where { uu=>uu(k) });

%%%% assume periodicity ll (should be even for now)
ll:=6;
let { uu(k+~d)=>uu(k+d-ll) when d>ll/2
    , uu(k-~d)=>uu(k-d+ll) when d>=ll/2 };
matrix ssm(ll,ll);
for i:=1:ll do begin
  ssm(i,i):=4/6;
  if i<ll then ssm(i,i+1):=ssm(i+1,i):=1/6;
  if i=ll then ssm(1,i):=ssm(i,1):=1/6;
end;
ssm:=1/ssm; % get the matrix of the ss operator
checkEquals6uuk:=(ss(uu(k-1)+4*uu(k)+uu(k+1),j) where
    ss(~a,j)=>ssm(1,1)*a+(for i:=1:ll/2 sum (ssm(i+1,1)*sub(k=k+i,a)))
                      +(for i:=1:(ll-1)/2 sum (ssm(i+1,1)*sub(k=k-i,a)))
    );

%% Can we find some nontrivial relations, first expand
jmpy:=(jmpy where {md(~a,j)=>(sub(k=k+1,a)-sub(k=k-1,a))/2
    , d2(~a,j)=>sub(k=k+1,a)-2*a+sub(k=k-1,a) 
    , ss(~a,j)=>ssm(1,1)*a+(for i:=1:ll/2 sum (ssm(i+1,1)*sub(k=k+i,a)))
                      +(for i:=1:(ll-1)/2 sum (ssm(i+1,1)*sub(k=k-i,a)))
    });
%% form equations
as:=for i:=0:nas collect mkid(a,i);
eqns:=for i:=-ll/2:ll/2 join coeffn(jmpy,uu(k+i),2).(
      for j:=-ll/2:ll/2 collect 
      coeffn(coeffn(jmpy,uu(k+i),1),uu(k+j),1));
soln:=solve(eqns,as);
factor arbcomplex; on list;
jmp:=sub(soln,jmp);

showtime;
end;
