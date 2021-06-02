% Laib04: approssimazione, interpolazione
clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_4_1*****************')
disp('***********************************************')

a = 0;
b = 1;
ord = 2:20;
for n = 1:length(ord)
    x = linspace(a,b,ord(n));
    A = vander(x);
    K1(n) = cond(A,1);
    K2(n) = cond(A,2);
    Kinf(n) = cond(A,inf);
end
figure(1)
semilogy(ord,K1,'*',ord,K2,'o',ord,Kinf,'+','linewidth',3)
legend('K_1','K_2','K_\infty')
xlabel('ordine sistema')
ylabel('condizionamento')
for n = 1:length(ord)
    x = linspace(a,b,ord(n));
    A = vander(x);
    t = sum(A,2);
    z = A\t;
    err(n) = norm(z-ones(ord(n),1))/norm(ones(ord(n),1));
end
figure(2)
semilogy(ord,err,'o','linewidth',3)
xlabel('ordine sistema')
ylabel('errore')

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_4_2*****************')
disp('***********************************************')

a = 0;
b = 1;
f = @(x) sin(pi*x);
z = linspace(a,b);
ord = 2:20;
for n = 1:length(ord)
    x = linspace(a,b,ord(n));
    y = f(x);
    c = polyfit(x,y,ord(n)-1);
    p = polyval(c,z);
    plot(z,f(z),'r',x,y,'g*',z,p,'b','linewidth',3)
    axis([0 1 0 1])
    pause
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause
clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_4_3*****************')
disp('***********************************************')

a = -5;
b = 5;
f = @(x) 1./(1+x.^2);
z = linspace(a,b);
for n = 5:5:20
    x = linspace(a,b,n+1);
    y = f(x);
    c = polyfit(x,y,n);
    p = polyval(c,z);
    err_p = max(abs(f(z)-p));
    
    fprintf('n = %d, err_p = %e\n',n,err_p);
    plot(x,y,'g*',z,f(z),'r',z,p,'linewidth',3)
    hold on
    legend('dati','funzione di Runge','polinomio interpolante')
    pause
end

disp('**************************************************************************************')
pause
close all

a = 1;
b = 2;
f = @(x) 1./(1+x.^2);
z = linspace(a,b);
for n = 5:5:20
    x = linspace(a,b,n+1);
    y = f(x);
    c = polyfit(x,y,n);
    p = polyval(c,z);
    err_p = max(abs(f(z)-p));
    
    fprintf('n = %d, err_p = %e\n',n,err_p);
    plot(x,y,'g*',z,f(z),'r',z,p,'linewidth',3)
    hold on
    legend('dati','funzione di Runge','polinomio interpolante')
    pause 
end

disp('**************************************************************************************')
pause
close all
a = -5;
b = 5;
f = @(x) 1./(1+x.^2);
z = linspace(a,b);
for n = 5:5:20
    u = -cos((2*[1:n+1]-1)*pi/(2*(n+1)));
    x = (b-a)/2*u+(b+a)/2;
    y = f(x);
    c = polyfit(x,y,n);
    p = polyval(c,z);
    err_p = max(abs(f(z)-p));
    fprintf('n = %d, err_p = %e\n',n,err_p);
    plot(x,y,'g*',z,f(z),'r',z,p,'linewidth',3)
    hold on
    legend('dati','funzione di Runge','polinomio interpolante')
    pause
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_4_4*****************')
disp('***********************************************')

format short e
f = @(x) sqrt(x);
x = 1.8:0.1:2.2;
y = [1.341641,1.378405,1.414214,1.449138,1.483240];
a = difdiv(x,y);
z = 1.85:0.1:2.15;
p = interp(x,a,z);
err = abs(f(z)-p)
zz = linspace(1.8,2.2);
pzz = interp(x,a,zz);
figure(1)
plot(x,y,'g*',zz,pzz,'b',zz,f(zz),'r','linewidth',2)
figure(2)
semilogy(z,err,'k*','linewidth',2)

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_4_6*****************')
disp('***********************************************')

x = 0:3;
y = [1 4 8 16];
c = polyfit(exp(x),y,3);
z = linspace(0,3);
p = polyval(c,exp(z));
plot(x,y,'g*',z,p,'b','linewidth',2)

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
