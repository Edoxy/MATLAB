%esercizio_10_3_pendolo_doppio
clear all
close all
clc
g = 9.81;
l = 1;
f = @(t,y) [6*((2*y(3)-3*cos(y(1)-y(2))*y(4))/(16-9*cos(y(1)-y(2))^2));...
            6*((8*y(4)-3*cos(y(1)-y(2))*y(3))/(16-9*cos(y(1)-y(2))^2));...
         -1/2*(6*((2*y(3)-3*cos(y(1)-y(2))*y(4))/(16-9*cos(y(1)-y(2))^2))*...
         6*((8*y(4)-3*cos(y(1)-y(2))*y(3))/(16-9*cos(y(1)-y(2))^2))*sin(y(1)-y(2))+3*g/l*sin(y(1)));... 
         -1/2*(-6*((2*y(3)-3*cos(y(1)-y(2))*y(4))/(16-9*cos(y(1)-y(2))^2))*...
         6*((8*y(4)-3*cos(y(1)-y(2))*y(3))/(16-9*cos(y(1)-y(2))^2))*sin(y(1)-y(2))+g/l*sin(y(2)))];
T = 10;
tspan = [0,T];
y0 = [pi/2;pi/2;0;0];
options = odeset('AbsTol',1.0e-8,'RelTol',1.0e-8,'stats','on');
[t,y] = ode45(f,tspan,y0,options);
N = length(t)-1;
%N = 250;
%[t,y] = Eulero_esp_system(f,tspan(1),y0,tspan(2),N);
%[t,y] = Runge_Kutta_4_system(f,tspan(1),y0,tspan(2),N);
figure(1)
plot(t,l*sin(y(:,1)),'r','linewidth',2)
title('spostamento ascissa m1')
figure(2)
plot(t,-l*cos(y(:,1)),'b','linewidth',2)
title('spostamento ordinata m1')
pause
figure(3)
for i = 1:N+1
    xp1 = l*sin(y(i,1));
    yp1 = -l*cos(y(i,1));
    xp2 = xp1 + l*sin(y(i,2));
    yp2 = yp1 - l*cos(y(i,2));
    plot(0,0,'r.',[0 xp1 xp2],[0 yp1 yp2],'k',xp1,yp1,'.g',xp2,yp2,'.g','linewidth',3,'markersize',40)
    axis([-2*l 2*l -2*l 2*l])
    title([num2str(t(i)),' s'])
    pause(0.001)
end