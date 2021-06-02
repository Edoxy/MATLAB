function s = spline_vincolata_partizione_uniforme(x,y,di,df,t)
%
% s = spline_vincolata_partizione_uniforme(x,y,di,df,t)
%
% x : vettore contenente le ascisse (equidistanti) dei dati di interpolazione
% y : vettore contenente le ordinate dei dati di interpolazione
% di : valore della derivata f'(x(1))
% df : valore della derivata f'(x(end))
% t : vettore contenente i punti dell'intervallo di interpolazione 
%     in cui si vuole valutare la spline
%
% s : vettore contenente i valori che la spline vincolata assume nei punti memorizzati in t
%
n = length(x)-1;
h = x(2)-x(1);
d = [2*h 4*h*ones(1,n-1) 2*h];
c = h*ones(1,n);
b(1) = 6*(y(2)-y(1))/h-6*di;
b(2:n) = 6*(y(3:n+1)-y(2:n))./h-6*(y(2:n)-y(1:n-1))./h;
b(n+1) = 6*df-6*(y(n+1)-y(n))/h;
M = gauss_tridiag_nopiv(d,c,b);
s = zeros(1,length(t));
for k = 1:length(t)
    i = 1;
    while t(k) > x(i+1)
          i = i+1;
    end
    % t(k) appartiene all'intervallo [ x(i),x(i+1) )
    s(k) = ((x(i+1)-t(k))^3*M(i)+(t(k)-x(i))^3*M(i+1))/(6*h)+...
           ((y(i+1)-y(i))/h+h/6*(M(i)-M(i+1)))*(t(k)-x(i))+y(i)-h^2/6*M(i);
end


