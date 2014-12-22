% try to find nontrivial relations between operators
% AJR Dec 2014

on div; off allfac; on revpri;
operator ss;
operator md;
operator d2;
operator uu; 
%jmp:=a0*uu+a1*ss(uu,j)+a2*d2(uu,j)+a3*md(uu,j)+a4*ss(d2(uu,j),j);
nas:=21;
jmp:=a0*ss(ss(md(uu(k),j),j)*uu(k),j)
    +a1*ss(md(uu(k),j)*ss(uu(k),j),j)
    +a2*ss(md(uu(k),j),j)*uu(k)
    +a3*ss(md(ss(uu(k),j)*uu(k),j),j)
    +a4*ss(md(uu(k)^2,j),j)
    +a5*uu(k)^2
    +a6*ss(uu(k),j)*uu(k)
    +a7*ss(ss(uu(k),j)*uu(k),j)
    +a8*ss(md(uu(k),j)*ss(md(uu(k),j),j),j)
    +a9*ss(d2(uu(k),j)*uu(k),j)
    +a10*ss(md(ss(md(uu(k),j),j)*uu(k),j),j)
    +a11*ss(uu(k)^2,j)
    +a12*ss(md(uu(k),j)*uu(k),j)
    +a13*ss(d2(uu(k),j)*ss(md(uu(k),j),j),j)
    +a14*ss(ss(ss(uu(k),j)*uu(k),j),j)
    +a15*ss(ss(ss(md(uu(k),j),j)*uu(k),j),j)
    +a16*ss(ss(md(uu(k),j)*ss(uu(k),j),j),j)
    +a17*ss(ss(d2(uu(k),j)*ss(uu(k),j),j),j)
    +a18*ss(ss(md(ss(uu(k),j)*uu(k),j),j),j)
    +a19*ss(ss(md(ss(md(uu(k),j),j)*uu(k),j),j),j)
    +a20*ss(ss(md(uu(k)^2,j),j),j)
    +a21*ss(ss(uu(k)^2,j),j)
    ;
%jmpy:=(jmp where { uu=>uu(k) });
jmpy:=jmp;

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

off list;
off nat;
out "invariants.out";
i1c := coeffn(jmp,arbcomplex(1),1);
i2c := coeffn(jmp,arbcomplex(2),1);
i3c := coeffn(jmp,arbcomplex(3),1);
i4c := coeffn(jmp,arbcomplex(4),1);
i5c := coeffn(jmp,arbcomplex(5),1);
shut "invariants.out";
on nat;

% Check invariants reduce jmp to zero:
i1h :=  -ss(ss(md(uu(k),j),j)*uu(k),j) 
  + 1/2*ss(md(uu(k),j)*ss(uu(k),j),j) 
  + 3/2*ss(md(uu(k),j),j)*uu(k) 
  + ss(md(ss(uu(k),j)*uu(k),j),j) 
  - 3/2*ss(md(uu(k)**2,j),j);

i2h :=  -18*uu(k)**2 
  + 9*ss(uu(k),j)*uu(k) 
  - 6*ss(ss(uu(k),j)*uu(k),j) 
  + ss(md(uu(k),j)*ss(md(uu(k),j),j),j) 
  + 3/2*ss(d2(uu(k),j)*uu(k),j) 
  + 2*ss(md(ss(md(uu(k),j),j)*uu(k),j),j) 
  + 15*ss(uu(k)**2,j);


i3h := 6*ss(md(uu(k),j)*uu(k),j) 
  - 3*ss(md(uu(k),j)*ss(uu(k),j),j) 
  + ss(d2(uu(k),j)*ss(md(uu(k),j),j),j) 
  + 3*ss(md(uu(k),j),j)*uu(k) 
  - 3*ss(md(uu(k)**2,j),j);

i4h := 12*ss(ss(uu(k),j)*uu(k),j) 
  - 6*ss(ss(ss(uu(k),j)*uu(k),j),j) 
  + ss(ss(d2(uu(k),j)*ss(uu(k),j),j),j) 
  + 2*ss(ss(md(ss(md(uu(k),j),j)*uu(k),j),j),j) 
  + 12*ss(ss(uu(k)**2,j),j) 
  - 18*ss(uu(k)**2,j);

i5h := 3/2*ss(ss(md(uu(k),j),j)*uu(k),j) 
  - ss(ss(ss(md(uu(k),j),j)*uu(k),j),j) 
  + 1/2*ss(ss(md(uu(k),j)*ss(uu(k),j),j),j) 
  + ss(ss(md(ss(uu(k),j)*uu(k),j),j),j) 
  - 3/2*ss(ss(md(uu(k)**2,j),j),j);

% Note: i5h == ss(i1h,j).
i1c-i1h;
i2c-i2h;
i3c-i3h;
i4c-i4h;
i5c-i5h;

jmps := jmp$
jmps := (jmps where 
ss(md(ss(uu(k),j)*uu(k),j),j) =>
  ss(ss(md(uu(k),j),j)*uu(k),j)
  - 1/2*ss(md(uu(k),j)*ss(uu(k),j),j)
  - 3/2*ss(md(uu(k),j),j)*uu(k)
  + 3/2*ss(md(uu(k)^2,j),j)
);
jmps := (jmps where 
ss(md(uu(k),j)*ss(md(uu(k),j),j),j) =>
  18*uu(k)**2 
  - 9*ss(uu(k),j)*uu(k) 
  + 6*ss(ss(uu(k),j)*uu(k),j) 
  - 3/2*ss(d2(uu(k),j)*uu(k),j) 
  - 2*ss(md(ss(md(uu(k),j),j)*uu(k),j),j) 
  - 15*ss(uu(k)**2,j)
);
jmps := (jmps where 
ss(d2(uu(k),j)*ss(md(uu(k),j),j),j) =>
  - 6*ss(md(uu(k),j)*uu(k),j) 
  + 3*ss(md(uu(k),j)*ss(uu(k),j),j) 
  - 3*ss(md(uu(k),j),j)*uu(k) 
  + 3*ss(md(uu(k)**2,j),j)
);
jmps := (jmps where 
ss(ss(md(ss(md(uu(k),j),j)*uu(k),j),j),j) =>
  - 6*ss(ss(uu(k),j)*uu(k),j) 
  + 3*ss(ss(ss(uu(k),j)*uu(k),j),j) 
  - 1/2*ss(ss(d2(uu(k),j)*ss(uu(k),j),j),j) 
  - 6*ss(ss(uu(k)**2,j),j) 
  + 9*ss(uu(k)**2,j)
);
%% 5: erk!!
%jmps := (jmps where 
%ss(ss(md(ss(uu(k),j)*uu(k),j),j),j) =>
%  - 3/2*ss(ss(md(uu(k),j),j)*uu(k),j) 
%  + ss(ss(ss(md(uu(k),j),j)*uu(k),j),j) 
%  - 1/2*ss(ss(md(uu(k),j)*ss(uu(k),j),j),j) 
%  + 3/2*ss(ss(md(uu(k)**2,j),j),j)
%);

end;
