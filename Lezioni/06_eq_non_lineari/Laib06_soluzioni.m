% Laib06: equazioni e sistemi non lineari
clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_6_1*****************')
disp('***********************************************')

syms x % variabile x simbolica
% costruzione simbolica del polinomio 
% p(x) = (x-1)(x-2)...(x-20)
e = prod(x-[1:20]); 
e = expand(e); 
c = coeffs(e);

% in c i coefficienti sono ordinati inversamente 
% rispetto alla rappresentazione polinomiale di Matlab 
% (c(1)x^20 + c(2)x^19 + ...)
c = flip(c);

% calcolo degli zeri del polinomio p(x) con 
% i coefficienti ''esatti''
x_exact = roots(c);
x_exact = flip(x_exact);

% ritorno da simbolico ad aritmetica reale
c = double(c); 

% calcolo degli zeri del polinomio p(x) con
% i coefficienti ''perturbati''
x = roots(c);

% errore relativo commesso su ciascuna radice
x_exact_d = double(x_exact);
err_rel = max((abs(x-x_exact_d))./((abs(x_exact_d))));
disp(['Il max errore relativo associato agli zeri del polinomio perturbato è pari a ',num2str(err_rel)])

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_6_2*****************')
disp('***********************************************')

nmax = 100;
tol = 1.0e-10;
%
disp('****************funzione 1****************') 
f = @(x) x.^2-2;
fd = @(x) 2*x;
figure
x = linspace(-2,2);
plot(x,f(x),'b',x,0*x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ');
[x,ier] = newton(f,fd,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)

% calcolo dell'ordine sperimentale di convergenza
osc = [0 log(e(2:end))./log(e(1:end-1))];
% stampa delle iterate calcolate, dei corrispondenti
% errori e dell'ordine sperimentale di convergenza
fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
pause
%
disp('****************funzione 2****************')
f = @(x) x.^3-x-1;
fd = @(x) 3*x.^2-1;
figure
x = linspace(-2,2);
plot(x,f(x),'b',x,0*x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ');
[x,ier] = newton(f,fd,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
pause
%
disp('****************funzione 3****************')
f = @(x) (x-2.^(-x)).^3;
fd = @(x) 3*(x-2.^(-x)).^2*(1+2.^(-x).*log(2));
figure
x = linspace(-2,2);
plot(x,f(x),'b',x,0*x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ');
[x,ier] = newton(f,fd,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
pause
%
disp('****************funzione 4****************')
f = @(x) exp(x)-2*x.^2;
fd = @(x) exp(x)-4*x;
figure
x = linspace(-2,2);
plot(x,f(x),'b',x,0*x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ')
[x,ier] = newton(f,fd,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_6_3*****************')
disp('***********************************************')

nmax = 100;
tol = 1.0e-10;
disp('****************funzione 1****************')
g = @(x) -sqrt(exp(x)/2);
figure
x = linspace(-2,2);
plot(x,g(x),'b',x,x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ')
[x,ier] = punto_fisso(g,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
pause
%
disp('****************funzione 2****************')
g = @(x) (2*x.^3+4*x.^2+10)./(3*x.^2+8*x);
figure
x = linspace(1,2);
plot(x,g(x),'b',x,x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ')
[x,ier] = punto_fisso(g,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_6_4*****************')
disp('***********************************************')

nmax = 100;
tol = 1.0e-10;
%
disp('****************funzione 1****************')
g = @(x) -log(x);
figure
x = linspace(0.2,0.8);
plot(x,g(x),'b',x,x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ')
[x,ier] = punto_fisso(g,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
% non converge |g'(csi)|>1 (le iterate sono complesse)
pause
%
disp('****************funzione 2****************')
g = @(x) exp(-x);
figure
x = linspace(0.2,0.8);
plot(x,g(x),'b',x,x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ');
[x,ier] = punto_fisso(g,x0,nmax,tol)
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
% converge |g'(csi)| è circa 0.607
pause
%
disp('****************funzione 3 (k=2) ****************')
g = @(x) (x+exp(-x))/2;
figure
x = linspace(0.2,0.8);
plot(x,g(x),'b',x,x,'r','linewidth',3)
grid on
x0 = input('fornire x0 = ');
[x,ier] = punto_fisso(g,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
% converge |g'(csi)| è circa 0.197
pause
%
disp('****************funzione 4 (k=3/2) ****************')
g = @(x) (x+2*exp(-x))/3;
figure
x = linspace(0.2,0.8);
plot(x,g(x),'b',x,x,'r','linewidth',3)
grid on
x0 = input('fornire x0 = ');
[x,ier] = punto_fisso(g,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end
% converge |g'(csi)| è circa 0.071
pause
%
disp('****************funzione 5****************')
% newton
f = @(x) x+log(x);
fd = @(x) 1+(1/x);
figure
x = linspace(0.2,0.8);
plot(x,f(x),'b',x,0*x,'r','linewidth',3)
grid on
x0 = input('fornire x0 = ');
[x,ier] = newton(f,fd,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_6_5*****************')
disp('***********************************************')

nmax = 100;
tol = 1.0e-10;
f = @(x) sqrt(pi)/2*erf(x)-1/2;
fd = @(x) exp(-x.^2);
figure
x = linspace(-2,2);
plot(x,f(x),'b',x,0*x,'r','linewidth',3)
grid on
x0 = input('fornisci x0 = ')
[x,ier] = newton(f,fd,x0,nmax,tol);
figure
n = length(x);
e = abs(x(2:n)-x(1:n-1));
semilogy(1:n-1,e,'linewidth',3)
pause
osc = [0 log(e(2:end))./log(e(1:end-1))];

fprintf('%10s\t%10s\t\t%10s\t%10s\r\n','n','approssimazione','errore','p');
for i = 1:n-1
    fprintf('%10.0d\t%6.16e\t%6.4d\t%6.4e\n',i,x(i),e(i),osc(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
