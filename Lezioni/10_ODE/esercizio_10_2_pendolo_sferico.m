%esercizio_10_2_pendolo_sferico
clear all
close all
clc
g = 9.81;
l = 1;
f = @(t,y) [y(2); y(4)^2*sin(y(1))*cos(y(1))-g/l*sin(y(1)); 
            y(4); -2*y(2)*y(4)*cot(y(1))];
%y0 = [sqrt(2)/2;sqrt(2)/2;1;0];
%T = 10;
%y0 = [0.2;0;0;1];
%T = 221;
y0 = [1;0;0;0.17];
T = 63.2;
tspan = [0 T];
options = odeset('AbsTol',1.0e-8,'RelTol',1.0e-8,'stats','on');
[x,y] = ode45(f,tspan,y0,options);
%N = 10^5;
%[x,y] = Runge_Kutta_4_system(f,tspan(1),y0,tspan(end),N);
x1 = l*sin(y(:,1)).*cos(y(:,3));
x2 = l*sin(y(:,1)).*sin(y(:,3));
x3 = -l*cos(y(:,1));
figure(1)
plot(y(:,1),y(:,2),'linewidth',2)
title('piano delle fasi (theta)')
figure(2)
plot(y(:,3),y(:,4),'linewidth',2)
title('piano delle fasi (phi)')
figure(3)
plot3(x1,x2,x3,'linewidth',2)
pause
for i = 1:length(x)
     plot3(x1(1:i),x2(1:i),x3(1:i),'b',x1(i),x2(i),x3(i),'r.','markersize',24,'linewidth',2)
     axis([min(x1) max(x1) min(x2) max(x2) min(x3) max(x3)])
     pause(.001)
end