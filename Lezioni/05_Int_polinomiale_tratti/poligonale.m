function [p,ier] = poligonale(x,y,t)
%
% [p,ier] = poligonale(x,y,t)
%
% x : vettore contenente le ascisse (in ordine crescente) dei dati di interpolazione
% y : vettore contenente le ordinate dei dati di interpolazione
% t : punto in cui si vuole valutare la poligonale
%
% p : vettore contenente il valore che la poligonale assume in t
% ier : vale 0 se t non appartiene a [x(1),x(n)],
%       vale 1 se t appartiene a [x(1),x(n)] e, conseguentemente, è stato 
%       calcolato il valore della poligonale in t
%
n = length(x);
if t < x(1) || t > x(n)
    p = NaN;
    ier = 0;
    return
end
i = 1;
while t > x(i+1)
    i = i+1;
end
i
% t appartiene all'intervallo (x(i),x(i+1)]
p = (t-x(i+1))/(x(i)-x(i+1))*y(i)+(t-x(i))/(x(i+1)-x(i))*y(i+1);
ier = 1;