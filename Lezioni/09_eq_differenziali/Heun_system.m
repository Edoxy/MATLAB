function [x,y] = Heun_system(f,x_0,y_0,x_N,N)
%
% Uso:     [x,y] = Heun_system(f,x_0,y_0,x_N,N)
% Scopo:   risolve un sistema di m equazioni differenziali 
%          y'(x) = f(x,y(x)) con condizioni iniziali y(x_0)=y_0
%          mediante il metodo di Heun
% Input:   f = funzione nota f(x,y) del problema 
%        x_0 = punto iniziale
%        y_0 = vettore (colonna) contenente gli m valori iniziali
%        x_N = punto finale
%          N = controlla l'ampiezza del passo h
% Output:  x = vettore (colonna) contenente N+1 nodi dell'intervallo 
%              [x_0,x_N], equidistanti con passo h=(x_N-x_0)/N 
%          y = matrice N+1 x m contenente le approssimazioni di y_1(x),
%              y_2(x),...,y_m(x) nei punti del vettore x
%
h = (x_N-x_0)/N;
x = (x_0:h:x_N)';
m = length(y_0);
y = zeros(m,N+1);
y(:,1) = y_0;
for n = 1:N
    f1 = f(x(n),y(:,n));
    f2 = f(x(n)+h,y(:,n)+h*f1); 
    y(:,n+1) = y(:,n)+h/2*(f1+f2);
end
y = y';