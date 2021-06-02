% esercizio_11_1
clear
clc
% dati equazione e metodo
C = input('fornisci C = ');
f = @(t,P) C*P;
sol = @(x) exp(C*x);
t0 = 0;
P0 = 1;
T = 20;
z = linspace(t0,T);
N = 10000;
[x,y] = Eulero_esp(f,t0,P0,T,N);
semilogy(x,sol(x),'g',x,y,'b','linewidth',2)
legend('soluzione esatta','soluzione approssimata') 