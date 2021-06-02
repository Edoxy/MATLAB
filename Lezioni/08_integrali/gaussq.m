function [d,z]=gaussq(ic,n,alpha,beta,kpts,endpts)
%  calcolo di zeri e pesi di polinomi ortogonali
%   ic = 1, Polinomi di Legendre w(x) =1
%   ic = 2, Polinomi di Chebychev di I tipo w(x)=1/sqrt(1-x2)
%   ic = 3, Polinomi di Chebychev di II tipo w(x) = sqrt(1-x^2)
%   ic = 4, Polinomi di Hermite w(x) = exp(-x^2)
%   ic = 5, Polinomi di Jacobi w(x) =(1-x)^alpha*(1+x)^beta
%   ic = 6, Polinomi di Laguerre w(x) = exp(-x)
%
%
[a1,b1,mu1]=class1(ic,n,alpha,beta);
if kpts == 1,
   [gam]=solve(endpts(1),n,a1,b1);
   a1(n)=gam*b1(n-1)^2+endpts(1);
elseif kpts==2,
   [gam]=solve(endpts(1),n,a1,b1);
   [gam1]=solve(endpts(2),n,a1,b1);
   t1=((endpts(1)-endpts(2))/(gam1-gam));
   b1(n-1)=sqrt(t1);
   a1(n)=gam*t1+endpts(1);
end


[d,e,z,ier]=gausq2(n,a1,b1);
z=z.^2*mu1;


