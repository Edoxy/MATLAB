function  a = difdiv(x,y)
%
% a = difdiv(x,y)
%
% x : vettore contenente le ascisse dei dati di interpolazione
% y : vettore contenente le ordinate dei dati di interpolazione
%
% a : vettore contenente i coefficienti della rappresentazione 
%     di Newton del polinomio interpolante, ordinati dal termine 
%     di grado 0 a quello di grado massimo
%
n = length(x)-1;
for i = 1:n
    for j = n+1:-1:i+1
        y(j) = (y(j)-y(j-1))/(x(j)-x(j-i)); 
    end
end
a = y;
