function [x,ier] = punto_fisso(g,x0,nmax,tol)
%
% [x,ier] = punto_fisso(g,x0,nmax,tol)
%
%   g   :   funzione non lineare di cui si vuole 
%           calcolare un punto fisso
%   x0  :   approssimazione iniziale
%   nmax:   numero massimo di iterazioni consentite
%   x   :   vettore contenente le approssimazioni calcolate
%   ier :   indicatore del criterio d'arresto utilizzato, vale 0 se raggiunte nmax iterazioni, 
%           vale 1 se soddisfatta la tolleranza relativa tol
%

ier = 0;
x(1) = x0;
for n = 1:nmax
    x(n+1) = g(x(n));
    if abs(x(n+1)-x(n)) <= tol*abs(x(n+1))
        ier = 1;
        break
    end
end