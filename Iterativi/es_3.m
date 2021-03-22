clear all
clc
%devo redere A simmetrica definita positiva
%basta prendere B non singolare e A = B'*B è simmetrica definita positiva
B = rand(5);
A = B'*B;
b = sum(A, 2);
n = length(b);
x = zeros(n, 1);
kmax= 10000;
toll = 1.0e-07;
[x, k, ier] = gradiente(A, b, x, toll, kmax)
%molto mal condizionata infatti:
cond(A)