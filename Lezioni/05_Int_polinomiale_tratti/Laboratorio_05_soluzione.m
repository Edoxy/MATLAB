% Laib05: interpolazione, funzioni polinomiali a tratti
clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_5_1*****************')
disp('***********************************************')

a = 0;
b = 1;
n = 20;
x = linspace(a,b,n);
f = @(x) sin(x);
y = f(x);
z = linspace(a,b);
p = zeros(1,length(z));
for k = 1:length(z)
    [p(k),ier] = poligonale(x,y,z(k));
end
plot(x,y,'g*',z,f(z),'r',z,p,'b','linewidth',3)
hold on
p_matlab = interp1(x,y,z);
plot(z,p_matlab,'c','linewidth',3)
err_rel = norm(p-p_matlab,inf)/norm(p_matlab,inf);
disp(['L''errore relativo in norma infinito in 100 punti è pari a ',num2str(err_rel)])

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_5_2*****************')
disp('***********************************************')

format short e
a = 0;
b = 1;
n = 21;
x = linspace(a,b,n);
f = @(x) sin(x);
y = f(x);
z = linspace(a,b);

p1 = zeros(1,length(z));
p2 = zeros(1,length(z));
p4 = zeros(1,length(z));
for k = 1:length(z)
    [p1(k),ier] = polinomiale_a_tratti(x,y,1,z(k));
    [p2(k),ier] = polinomiale_a_tratti(x,y,2,z(k));
    [p4(k),ier] = polinomiale_a_tratti(x,y,4,z(k));
end
plot(x,y,'g*',z,f(z),'r',z,p1,'c',z,p2,'b',z,p4,'k','linewidth',3)
legend('dati','funzione','d=1','d=2','d=4')
err1 = norm(f(z)-p1,inf);
err2 = norm(f(z)-p2,inf);
err3 = norm(f(z)-p4,inf);
disp(['L''errore in norma infinito in 100 punti per d=1 è pari a ',num2str(err1)])
disp(['L''errore in norma infinito in 100 punti per d=2 è pari a ',num2str(err2)])
disp(['L''errore in norma infinito in 100 punti per d=4 è pari a ',num2str(err3)])

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_5_3*****************')
disp('***********************************************')

a = 0;
b = 1;
n = 20;
x = linspace(a,b,n);
f = @(x) sin(x);
y = f(x);
fd = @(x) cos(x);
di = fd(a);
df = fd(b);
z = linspace(a,b,100);
s = spline_vincolata_partizione_uniforme(x,y,di,df,z);
plot(x,y,'g*',z,f(z),'r',z,s,'b','linewidth',2)
hold on
s_matlab = spline(x,[di y df],z);
plot(z,s_matlab,'c','linewidth',2)
err_rel = norm(s-s_matlab,inf)/norm(s_matlab,inf);
disp(['L''errore relativo in norma infinito in 100 punti è pari a ',num2str(err_rel)])


disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_5_4*****************')
disp('***********************************************')

a = -5;
b = 5;
f = @(x) 1./(1+x.^2);
z = linspace(a,b);
for n = 5:5:15
    x = linspace(a,b,n);
    y = f(x);
    s = spline(x,y,z);
    err_s = max(abs(f(z)-s));    
    fprintf('n = %d, err_s = %e\n',n,err_s);
    plot(x,y,'g*',z,f(z),'r',z,s,'linewidth',3)
    legend('dati','funzione di Runge','spline interpolante')
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
disp('*****************esercizio_5_5*****************')
disp('***********************************************')

f = @(x) (1-x.^2).^(5/2);
fd = @(x) (5/2)*(1-x.^2).^(5/2-1).*(-2*x);
di = fd(-1);
df = fd(1);
z = linspace(-1,1);
fz = f(z);
for k = 2:5
    n = 2^k;
    x = -1+2*(0:n)/n;
    y = f(x);
    s = spline(x,y,z);
    s1 = spline(x,[di y df],z);
    figure(1)
    plot(x,y,'ko',z,fz,'r',z,s,'b',z,s1,'g','linewidth',3)
    legend('dati','f(x)','spline not-a-knot','spline vincolata')
    pause
    figure(2)
    semilogy(z,abs(s-fz),'b',z,abs(s1-fz),'g','linewidth',3)
    legend('errore spline not-a-knot','errore spline vincolata')
    err = norm(fz-s,inf);
    err1 = norm(fz-s1,inf);
    fprintf('n = %d, err = %e, err1 = %e\n',n,err,err1);
    pause
end    
% la spline vincolata fornisce un'approssimazione più
% accurata della spline not-a-knot in quanto, a differenza di
% quest'ultima, oltre alle condizioni di interpolazione, soddisfa
% due ulteriori condizioni che stabiliscono un legame
% con la funzione f

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

