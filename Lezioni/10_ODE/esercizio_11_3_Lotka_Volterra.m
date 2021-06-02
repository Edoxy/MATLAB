% esercizio_11_3
clear
close all
clc
% dati equazione e metodo
C1 = 1;
C2 = -1;
d1 = 1;
d2 = -1;
f = @(t,P) [C1*P(1)+d2*P(1)*P(2); C2*P(2)+d1*P(1)*P(2)];
t0 = 0;
P0 = [2;2];  % P0 = [1.2;1.2];
T = 10;
z = linspace(t0,T);
% punto di equilibrio stabile
PES = [-C2/d1,-C1/d2];
for N = 100:200:600
    [x,y] = Eulero_esp_system(f,t0,P0,T,N);
    figure(1)
    plot(x,y(:,1),'r',x,y(:,2),'b','linewidth',2)
    legend('P_1- preda','P_2- predatore') 
    title('soluzioni approssimate')
    figure(2)
    for i = 1:length(x)
        plot(PES(1),PES(2),'.r',y(i,1),y(i,2),'g.',y(1:i,1),y(1:i,2),'g','linewidth',2,'markersize',40)
        m1 = min(y(:,1));
        M1 = max(y(:,1));
        m2 = min(y(:,2));
        M2 = max(y(:,2));
        axis([m1 M1 m2 M2])
        pause(0.01)
    end
    hold on
    % si rappresenta il campo vettoriale
    [X,Y] = meshgrid(m1:0.1:M1,m2:0.1:M2);
    u = C1*X+d2*X.*Y;
    v = C2*Y+d1*X.*Y;
    norma = 1;%sqrt(u.^2+v.^2);
    quiver(X,Y,u./norma,v./norma,'linewidth',2)    
    title('piano delle fasi') 
    pause
    hold off
end
options = odeset('AbsTol',1.0e-08,'RelTol',1.0e-08);
[x,y] = ode45(f,[t0 T],P0,options);
figure(3)
plot(x,y(:,1),'r',x,y(:,2),'b','linewidth',2)
legend('P_1- preda','P_2- predatore') 
title('soluzioni approssimate')
figure(4)
for i = 1:length(x)
    plot(PES(1),PES(2),'.r',y(i,1),y(i,2),'g.',y(1:i,1),y(1:i,2),'g','linewidth',2,'markersize',40)
    m1 = min(y(:,1));
    M1 = max(y(:,1));
    m2 = min(y(:,2));
    M2 = max(y(:,2));
    axis([m1 M1 m2 M2])
    pause(0.01)
end
hold on 
[X,Y] = meshgrid(m1:0.1:M1,m2:0.1:M2);
u = C1*X+d2*X.*Y;
v = C2*Y+d1*X.*Y;
norma = 1;%sqrt(u.^2+v.^2);
quiver(X,Y,u./norma,v./norma,'linewidth',2)    
title('piano delle fasi') 
hold off
