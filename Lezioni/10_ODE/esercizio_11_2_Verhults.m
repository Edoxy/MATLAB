% esercizio_11_2
clear
clc
% dati equazione e metodo
C = 2;
B = 100; 
f = @(t,P) C*P*(1-P/B);
t0 = 0;
P0 = 200;  %P0 = 2;
sol = @(x) B./(1+(B-P0)/P0*exp(-C*x));
T = 20;
z = linspace(t0,T);
for N = 20:20:400
    h = (T-t0)/N
    [x,y] = Eulero_esp(f,t0,P0,T,N);
    plot(x,sol(x),'g',x,y,'b','linewidth',2)
    legend('soluzione esatta','soluzione approssimata') 
    pause
end