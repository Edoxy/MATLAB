% Laib07: sistemi di equazioni non lineari
clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_7_1*****************')
disp('***********************************************')

% dati sistema e metodo
F1 = @(x) [x(1)^2+2*x(1)*x(2)+x(3); x(2)^3+x(3)^2; x(1)*x(3)-1];
J1 = @(x) [2*(x(1)+x(2)) 2*x(1) 1; 0 3*x(2)^2 2*x(3); x(3) 0 x(1)];
n_max = 20;
tol = 1.0e-10;
exact = [1;-1;1];
x0 = [0.5;-0.5;0.1];
[x,ier] = newton_system(F1,J1,x0,tol,n_max);
err = norm(exact-x);
disp(['err = ',num2str(err)])

% il metodo di Newton è sensibile alla scelta del punto iniziale 
x0 = [-0.5;-0.5;-0.1]; %converge a un'altra soluzione 
[x,ier] = newton_system(F1,J1,x0,tol,1000)

x0 = [-5;0;-5]; %converge a un'altra soluzione diversa dalle precedenti
[x,ier] = newton_system(F1,J1,x0,tol,1000)

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_7_2*****************')
disp('***********************************************')

d1 = 4.5;
d2 = 6;
r1 = 5;
r2 = 9;
syms a b
f1 = d1*cos(b)+d2*cos(a+b)-r1;
f2 = d1*sin(b)+d2*sin(a+b)-r2;
figure(1)
fimplicit(f1,[-3*pi,3*pi],'b','linewidth',2)
hold on
fimplicit(f2,[-3*pi,3*pi],'r','linewidth',2)
grid on
pause
figure(2)
fimplicit(f1,[-pi,pi],'b','linewidth',2)
hold on
fimplicit(f2,[-pi,pi],'r','linewidth',2)
grid on
pause
F2 = @(x) [d1*cos(x(2))+d2*cos(x(1)+x(2))-r1;
           d1*sin(x(2))+d2*sin(x(1)+x(2))-r2];
J2 = @(x) [-d2*sin(x(1)+x(2))  -d1*sin(x(2))-d2*sin(x(1)+x(2));
            d2*cos(x(1)+x(2))   d1*cos(x(2))+d2*cos(x(1)+x(2))];
   
x10 = [-0.4; 1.3];
x20 = [0.25; 0.75];

%options = optimoptions('fsolve','OptimalityTolerance',1.0e-10);
[x1_m,FVAL1,EXITFLAG1,OUTPUT1] = fsolve(F2,x10)%,options);
[x2_m,FVAL2,EXITFLAG2,OUTPUT2] = fsolve(F2,x20)%,options);

nmax = 20;
tol = 1.0e-4;
[x1,ier1] = newton_system(F2,J2,x10,tol,nmax);
[x2,ier2] = newton_system(F2,J2,x20,tol,nmax);

figure(3)
plot([0 d1*cos(x1(2)) d1*cos(x1(2))+d2*cos(x1(1)+x1(2))],[0 d1*sin(x1(2)) d1*sin(x1(2))+d2*sin(x1(1)+x1(2))],'b','linewidth',3)
hold on
plot([0 d1*cos(x2(2)) d1*cos(x2(2))+d2*cos(x2(1)+x2(2))],[0 d1*sin(x2(2)) d1*sin(x2(2))+d2*sin(x2(1)+x2(2))],'g','linewidth',3)
grid on
legend('soluzione 1','soluzione 2')

err1 = norm(x1_m-x1)/norm(x1_m);
err2 = norm(x2_m-x2)/norm(x2_m);
disp('errori relativi associati a due soluzioni in [-pi,pi] calcolate con il metodo di Newton')
disp(['err1 = ',num2str(err1)])
disp(['err2 = ',num2str(err2)])
disp(' ')
h = 2^(-30);
%J2_app = @(x) [d2*(cos(x(1)+h+x(2))-cos(x(1)+x(2)))/h   d1*(cos(x(2)+h)-cos(x(2)))/h + d2*(cos(x(1)+h+x(2))-cos(x(1)+x(2)))/h;  
%               d2*(sin(x(1)+h+x(2))-sin(x(1)+x(2)))/h    d1*(sin(x(2)+h)-sin(x(2)))/h + d2*(sin(x(1)+h+x(2))-sin(x(1)+x(2)))/h];
F21 = @(x1,x2) d1*cos(x2)+d2*cos(x1+x2)-r1;
F22 = @(x1,x2) d1*sin(x2)+d2*sin(x1+x2)-r1;
J2_app = @(x) [(F21(x(1)+h,x(2))-F21(x(1),x(2)))/h (F21(x(1),x(2)+h)-F21(x(1),x(2)))/h;
    (F22(x(1)+h,x(2))-F22(x(1),x(2)))/h (F22(x(1),x(2)+h)-F22(x(1),x(2)))/h];
[x1app,ier1app] = newton_system(F2,J2_app,x10,tol,nmax);
[x2app,ier2app] = newton_system(F2,J2_app,x20,tol,nmax);

err1app30 = norm(x1_m-x1app)/norm(x1_m);
err2app30 = norm(x2_m-x2app)/norm(x2_m);

disp('errori relativi associati a due soluzioni in [-pi,pi] calcolate con il metodo di Newton')
disp('ove si usa un''approssimazione della matrice jacobiana') 
disp(['h = 2^-30, err1app = ',num2str(err1app30)])
disp(['h = 2^-30, err2app = ',num2str(err2app30)])
disp('si ottengono errori aventi lo stesso ordine di grandezza di quelli precedenti')
disp(' ')

h = 2^(-50);
%J2_app = @(x) [d2*(cos(x(1)+h+x(2))-cos(x(1)+x(2)))/h   d1*(cos(x(2)+h)-cos(x(2)))/h + d2*(cos(x(1)+h+x(2))-cos(x(1)+x(2)))/h;  
%               d2*(sin(x(1)+h+x(2))-sin(x(1)+x(2)))/h    d1*(sin(x(2)+h)-sin(x(2)))/h + d2*(sin(x(1)+h+x(2))-sin(x(1)+x(2)))/h];
F21 = @(x1,x2) d1*cos(x2)+d2*cos(x1+x2)-r1;
F22 = @(x1,x2) d1*sin(x2)+d2*sin(x1+x2)-r1;
J2_app = @(x) [(F21(x(1)+h,x(2))-F21(x(1),x(2)))/h (F21(x(1),x(2)+h)-F21(x(1),x(2)))/h;
    (F22(x(1)+h,x(2))-F22(x(1),x(2)))/h (F22(x(1),x(2)+h)-F22(x(1),x(2)))/h];

[x1app,ier1] = newton_system(F2,J2_app,x10,tol,nmax);
[x2app,ier2] = newton_system(F2,J2_app,x20,tol,nmax);

err1app50 = norm(x1_m-x1app)/norm(x1_m);
err2app50 = norm(x2_m-x2app)/norm(x2_m);

disp('errori relativi associati a due soluzioni  in [-pi,pi] calcolate con il metodo di Newton')
disp('ove si usa un''approssimazione della matrice jacobiana teoricamente più accurata!!!') 
disp(['h = 2^-50, err1app = ',num2str(err1app50)])
disp(['h = 2^-50, err2app = ',num2str(err2app50)])
disp('si ottengono errori aventi ordine di grandezza maggiore a causa della cancellazione numerica')
disp(' ')

% formule di prostaferesi per eliminare la cancellazione numerica
% che si verifica al numeratore dei rapporti incrementali.

cmc1 = @(x) -2/h*sin(h/2)*sin((2*x(1)+h+2*x(2))/2);
cmc2 = @(x) -2/h*sin(h/2)*sin((h+2*x(2))/2);
sms1 = @(x) 2/h*cos((2*x(1)+h+2*x(2))/2)*sin(h/2);
sms2 = @(x) 2/h*cos((h+2*x(2))/2)*sin(h/2);

J2_appc = @(x) [d2*cmc1(x)   d1*cmc2(x) + d2*cmc1(x);  
                d2*sms1(x)   d1*sms2(x) + d2*sms1(x)];

[x1app,ier1] = newton_system(F2,J2_appc,x10,tol,nmax);
[x2app,ier2] = newton_system(F2,J2_appc,x20,tol,nmax);

err1app50 = norm(x1_m-x1app)/norm(x1_m);
err2app50 = norm(x2_m-x2app)/norm(x2_m);

disp('errore relativo associato a due soluzioni  in [-pi,pi] calcolate con il metodo di Newton')
disp('ove si usa un''approssimazione della matrice jacobiana teoricamente e numericamente più accurata') 
disp(['h = 2^-50, err1appc = ',num2str(err1app50)])
disp(['h = 2^-50, err2appc = ',num2str(err2app50)])
disp('la cancellazione numerica è stata eliminata e si ottengono errori aventi lo stesso ordine di grandezza di quelli che si ottengono con la matrice jacobiana')

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_7_3*****************')
disp('***********************************************')

x = -3:.01:3;
y = x;
[X,Y] = meshgrid(x,y);
Z = peaks(X,Y);
meshc(X,Y,Z)
axis([-3 3 -3 3 -10 10]);
syms x y
z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) - 1/3*exp(-(x+1).^2 - y.^2);
zx = diff(z,x);
zy = diff(z,y);
% tentativo (che fallisce!) di ottenere una soluzione simbolica del sistema
M = solve(zx,zy);

% si procede quindi numericamente
F3 = @(x) [(exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2))/3 + 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2) + exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2) - 6*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 - 2*x(1)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5);
      (2*x(2)*exp(- (x(1) + 1)^2 - x(2)^2))/3 + 50*x(2)^4*exp(- x(1)^2 - x(2)^2) - 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)*(x(1) - 1)^2 - 2*x(2)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5)];

x0 = [0.4;1.5]; 
options = optimoptions('fsolve','OptimalityTolerance',1.0e-10);
[sol,FVAL,EXITFLAG,OUTPUT] = fsolve(F3,x0,options);

z = @(x,y) 3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) - 1/3*exp(-(x+1).^2 - y.^2);
hold on
% rappresentazione grafica del massimo
plot3(sol(1),sol(2),z(sol(1),sol(2)),'b*','linewidth',3)
MAX = [sol(1),sol(2),z(sol(1),sol(2))]

% la matrice Jacobiana è stata calcolata con il comando simbolico jacobian!!
J3 = @(x) [ (2*exp(- (x(1) + 1)^2 - x(2)^2))/3 + 6*exp(- (x(2) + 1)^2 - x(1)^2) + 60*x(1)*exp(- x(1)^2 - x(2)^2) - 6*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 - 2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) - (exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2)^2)/3 + 12*x(1)^2*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 + 4*x(1)^2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) - 12*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2) - 4*x(1)*exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2),  4*x(1)*x(2)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) - (2*x(2)*exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2))/3 - 2*x(2)*exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2) - 100*x(1)*x(2)^4*exp(- x(1)^2 - x(2)^2) - 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2)*(2*x(2) + 2) + 6*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)*(x(1) - 1)^2; 4*x(1)*x(2)*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) - (2*x(2)*exp(- (x(1) + 1)^2 - x(2)^2)*(2*x(1) + 2))/3 - 2*x(2)*exp(- x(1)^2 - x(2)^2)*(30*x(1)^2 - 2) - 100*x(1)*x(2)^4*exp(- x(1)^2 - x(2)^2) - 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(1) - 2)*(2*x(2) + 2) + 6*x(1)*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)*(x(1) - 1)^2, (2*exp(- (x(1) + 1)^2 - x(2)^2))/3 - 6*exp(- (x(2) + 1)^2 - x(1)^2)*(x(1) - 1)^2 - 2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) + 200*x(2)^3*exp(- x(1)^2 - x(2)^2) - 200*x(2)^5*exp(- x(1)^2 - x(2)^2) - (4*x(2)^2*exp(- (x(1) + 1)^2 - x(2)^2))/3 + 4*x(2)^2*exp(- x(1)^2 - x(2)^2)*(10*x(1)^3 - 2*x(1) + 10*x(2)^5) + 3*exp(- (x(2) + 1)^2 - x(1)^2)*(2*x(2) + 2)^2*(x(1) - 1)^2];
tol = 1.0e-10;
nmax = 700;
% il metodo di Newton è sensibile alla scelta del valore iniziale x0
x0 = [0.4;1.5];  %prova con x0 = [0.5;1.5];
[x3,ier] = newton_system(F3,J3,x0,tol,nmax)
x3