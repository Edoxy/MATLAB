clear all
close all
clc
format short e
f = @(x,y) [y(2); -y(1)/(sqrt(y(1)^2+y(3)^2))^3; 
            y(4); -y(3)/(sqrt(y(1)^2+y(3)^2))^3];
x_0 = 0;
y_0 = [0.5;0.1;0;1];
%x_N = 3;
x_N = 1.23;
N = 1000;%70000;%1000;%70000;
[x,y] = eulero_esp_system(f,x_0,y_0,x_N,N);
%options = odeset('AbsTol',1.0e-8,'RelTol',1.0e-8,'stats','on');
%[x,y] = ode45(f,[x_0 x_N],y_0,options);
plot(y(:,1),y(:,3),'linewidth',2)
pause
for i = 1:length(x)
    plot(y(1:i,1),y(1:i,3),'b',y(i,1),y(i,3),'r.','linewidth',2,'markersize',30)
    axis([min(y(:,1)) max(y(:,1)) min(y(:,3)) max(y(:,3))])
    pause(0.001)
end
   