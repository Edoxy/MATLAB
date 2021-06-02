function [x,y] = Eulero_esp(f,x_0,y_0,x_N,N)
%
% Uso:     [x,y] = Eulero_esp(f,x_0,y_0,x_N,N)
% Scopo:   risolve l'equazione differenziale y'(x)=f(x,y(x)) 
%          con condizione iniziale y(x_0)=y_0,
%          mediante il metodo di Eulero esplicito
% Input:   f = funzione nota f(x,y) del problema 
%        x_0 = punto iniziale
%        y_0 = valore iniziale
%        x_N = punto finale
%          N = controlla l'ampiezza del passo h
% Output:  x = vettore (colonna) contenente N+1 nodi dell'intervallo
%              [x_0,x_N], equidistanti con passo h=(x_N-x_0)/N 
%          y = vettore colonna contenente le N+1 approssimazioni di y(x),
%              nei punti del vettore x
%
h = (x_N-x_0)/N;
x = (x_0:h:x_N)';
y = zeros(N+1,1);
y(1) = y_0;
for n = 1:N
   fn = f(x(n),y(n)); 
   y(n+1) = y(n)+h*fn;
end