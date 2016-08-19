% numerically find linear solutions of potential well problem
A=-1
H=pi
V=@(x) A*cos(2*x);
Npts=19

x=linspace(0,H,Npts);
j=2:Npts;
xj=x(j-1);
h=diff(x(1:2));
A0=sparse(j,j-1,1/h^2,Npts+1,Npts+1) ...
  +sparse(j,j+1,1/h^2,Npts+1,Npts+1) ...
  +sparse(j,j,-2/h^2-V(xj),Npts+1,Npts+1);
A0=full(A0);
B=full(sparse(j,j,1,Npts+1,Npts+1));

ws=[];
%phis=linspace(0,pi,21);
h=0.1,phis=(0:3)*h
for phi=phis
A=A0;
A([1,Npts+1],[1,Npts+1])=eye(2);
A(1,Npts)=-exp(-1i*phi);
A(Npts+1,2)=-exp(+1i*phi);
w=eig(A,B);
w(isinf(w))=[];
ws=[ws;-sort(-real(w))'];
end
normimagfreq=norm(imag(ws))
plot(phis,asinh(real(ws)),'*')
tickx=[0.5 1 2];
tickx=sort([tickx 0 -tickx -10*tickx -100*tickx]');
set(gca,'Ytick',asinh(tickx),'YtickLabel',num2str(tickx))

% find vector of derivatives
w=ws(:,1);
w=[1 0 0 0
-2/h^2 2/h^2 0 0
6/h^4 -8/h^4 2/h^4 0
-20/h^6 30/h^6 -12/h^6 2/h^6]*w;
derivCoeffs=[1 0 0 0 
 0 1 -h^2/12 h^4/90
 0 0 1 -h^2/6
 0 0 0 1]*w
