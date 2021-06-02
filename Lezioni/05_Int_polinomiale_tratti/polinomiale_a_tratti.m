function [p,ier] = polinomiale_a_tratti(x,y,d,t)
%
% [p,ier] = polinomiale_a_tratti(x,y,d,t)
%
% x : vettore contenente le ascisse (in ordine crescente) dei dati di interpolazione
%     contenente n = N*d+1 elementi
% y : vettore contenente le ordinate dei dati di interpolazione
% d : grado locale della funzione polinomiale a tratti
% t : punto in cui si vuole valutare la funzione polinomiale a tratti
%
% p : vettore contenente il valore che la funzione polinomiale a tratti assume in t
% ier : vale 0 se t non appartiene a [x(1),x(n)],
%       vale 1 se t appartiene a [x(1),x(n)] e, conseguentemente, è stato 
%       calcolato il valore della poligonale in t
%
n = length(x);
if rem(n-1,d) ~= 0
   p = NaN;
   ier = -1;
   return
end
if t < x(1) || t > x(n)
    p = NaN;
    ier = 0;
    return
end
i = 1;
while t > x(i*d+1)
    i = i+1;
end
% t appartiene all'intervallo (x((i-1)*d+1),x(i*d+1)]
%c = polyfit(x((i-1)*d+1:i*d+1),y((i-1)*d+1:i*d+1),d);
%p = polyval(c,t);
c = difdiv(x((i-1)*d+1:i*d+1),y((i-1)*d+1:i*d+1));
p = interp(x((i-1)*d+1:i*d+1),c,t);
ier = 1;
