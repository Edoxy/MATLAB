function  p = interp(x,a,z)
%
% p = interp(x,a,z)
%
% x : vettore contenente le ascisse dei dati di interpolazione
% a : vettore contenente i coefficienti della rappresentazione 
%     di Newton del polinomio interpolante, ordinati dal termine 
%     di grado 0 a quello di grado massimo (output di difdiv)
% z : vettore di punti in cui si vuole valutare il polinomio
%
% p : vettore contenente le valutazioni del polinomio nei
%     punti memorizzati in z
%
n = length(x)-1;
p = a(n+1)*ones(size(z));
for i = n:-1:1
    p = p.*(z-x(i))+a(i); 
end