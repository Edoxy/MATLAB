%esercizio_10_1_pendolo_semplice
clear all
close all
clc
g = 9.81;
l = 1;
f = @(t,y) [y(2); -g/l*sin(y(1))];
y0 = [pi/4;0];
T = 5;
tspan = [0,T];
options = odeset('AbsTol',1.0e-8,'RelTol',1.0e-8,'stats','on');
[t,y] = ode45(f,tspan,y0,options);
%N = 5000;
%[t,y] = Eulero_esp_system(f,tspan(1),y0,tspan(2),N);
%N = 500;
%[t,y] = Runge_Kutta_4_system(f,tspan(1),y0,tspan(2),N);
figure(1)
plot(t,y(:,1),'linewidth',2)
title('spostamento')
figure(2)
plot(t,y(:,2),'linewidth',2)
title('velocità')
figure(3)
plot(y(:,1),y(:,2),'linewidth',2)
title('piano delle fasi')
pause
figure(4)
for i = 1:length(t)
    xp1 = l*sin(y(i,1));  %l*cos(-pi/2+y(i))
    yp1 = -l*cos(y(i,1)); %l*sin(-pi/2+y(i))   
    plot(0,0,'r.','markersize',40)
    hold on
    plot([0 xp1],[0 yp1],'k',xp1,yp1,'.g','linewidth',3,'markersize',40)     
    axis([-2*l 2*l -2*l 2*l])
    title([num2str(t(i)),' s'])
    pause(0.01)
    hold off    
end
