% numerically find linear solutions of potential well problem
A=-1
H=pi
V=@(x) A*cos(2*x);
Npts=9

x=linspace(0,H,Npts);
j=2:Npts
xj=x(j-1)
h=diff(x(1:2))
A0=sparse(j,j-1,1/h^2,Npts+1,Npts+1) ...
  +sparse(j,j+1,1/h^2,Npts+1,Npts+1) ...
  +sparse(j,j,-2/h^2-V(xj),Npts+1,Npts+1);
A0=full(A0)
B=full(sparse(j,j,1,Npts+1,Npts+1))

ws=[];
phis=linspace(0,pi,21);
for phi=phis
A=A0;
A([1,Npts+1],[1,Npts+1])=eye(2);
A(1,Npts)=-exp(-1i*phi);
A(Npts+1,2)=-exp(+1i*phi);
w=eig(A,B);
w(isinf(w))=[]
ws=[ws;w'];
end
normimagfreq=norm(imag(ws))
plot(phis,real(ws),'*')
