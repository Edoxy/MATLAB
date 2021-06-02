% Laib08: integrali
clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_8_1*****************')
disp('***********************************************')

kmax = 128;
tol = 1.0e-10;
for j = 1:6
    disp(['integrale numero ',num2str(j)])
    switch j
      case 1 % esempio 1
           clearvars -except j kmax tol
           a = 0; 
           b = 1;
           f = @(x) exp(x);
           exact = exp(1)-1;
           % funzione di classe C^inf
      case 2 % esempio 2
           clearvars -except j kmax tol
           a = 0;
           b = 1;
           f = @(x) cos(x);
           exact = sin(1);
           % funzione di classe C^inf
      case 3 % esempio 3
           clearvars -except j kmax tol
           a = 0.01;
           b = 1.1;
           f = @(x) 1./x.^4;
           exact = 1/3*(10^6-1.1^(-3));
           % funzione di classe C^inf ma con singolarità vicina
      case 4 % esempio 4
           clearvars -except j kmax tol
           a = 0;
           b = 1;
           f = @(x) sqrt(x);
           exact = 2/3;
           % funzione di classe C^0
      case 5 % esempio 5
           clearvars -except j kmax tol
           a = 0;
           b = 1;
           f = @(x) sin(99*pi*x);
           x = linspace(a,b,1000);
           plot(x,f(x),'linewidth',3)
           exact = 2/(99*pi);
           % funzione di classe C^inf ma molto oscillante 
      case 6 % esempio 6
           clearvars -except j kmax tol
           a = 0;
           b = 1;
           f = @(x) sin(100*pi*x);
           x = linspace(a,b,1000);
           plot(x,f(x),'linewidth',3)
           exact = 0;
           % funzione di classe C^inf ma molto oscillante 
           % dispari rispetto a x =1/2
    end
    k = 1;
    for n = 2:2:kmax
        num(k) = n;
        [nodi,pesi] = gaussq(1,num(k),0,0,0,[0 0]');
        nodi_s = (b-a)/2*nodi+(b+a)/2;
        pesi_s = (b-a)/2*pesi;
        int(k) = pesi_s*f(nodi_s)';
        err(k) = abs(exact-int(k));%/abs(exact); 
        if err(k)<= tol
           break
        elseif n < kmax
           k= k+1;
        end
    end
    fprintf('%10s\t%10s\t\t%10s\r\n','n','approssimazione','errore');
    for i = 1:k
        fprintf('%10.0d\t%6.16e\t%6.4d\t\n',num(i),int(i),err(i)); 
    end
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
disp('*****************esercizio_8_2*****************')
disp('***********************************************')

kmax = 128;
tol = 1.0e-10;
f = @(x) x.^2.*sin(x.^2);
exact = 2^(1/4)/16*sqrt(pi*(sqrt(2)+2));

k = 1;
for n = 2:2:kmax
    num(k) = n;
    % formula di quadratura di Gauss-Hermite
    [nodi,pesi] = gaussq(4,n,0,0,0,[0 0]');
    int(k) = 1/2*pesi*f(nodi)';
    err(k) = abs(exact-int(k))/abs(exact); 
    if err(k)<= tol
       break
    elseif n < kmax
       k= k+1;
    end
end

fprintf('%10s\t%10s\t\t%10s\r\n','n','approssimazione','errore');
for i = 1:k
    fprintf('%10.0d\t%6.16e\t%6.4d\t\n',num(i),int(i),err(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_8_3*****************')
disp('***********************************************')

kmax = 128;
tol = 1.0e-10;
f = @(x) exp(-1)*(x+1).^(3/2);
exact = 5/(2*exp(1)) + 3/4*sqrt(pi)*(1-erf(1));

k = 1;
for n = 2:2:kmax
    num(k) = n;
    % formula di quadratura di Gauss-Laguerre
    [nodi,pesi] = gaussq(6,n,0,0,0,[0 0]');
    int(k) = pesi*f(nodi)';
    err(k) = abs(exact-int(k))/abs(exact); 
    if err(k)<= tol
       break
    elseif n < kmax
       k= k+1;
    end
end

fprintf('%10s\t%10s\t\t%10s\r\n','n','approssimazione','errore');
for i = 1:k
    fprintf('%10.0d\t%6.16e\t%6.4d\t\n',num(i),int(i),err(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_8_4*****************')
disp('***********************************************')

f = @(x) x.^9+x.^8+x.^7+x.^2+9;
exact = (1251*pi)/128;
% formula di quadratura di Gauss-Chebichev
n = 5;
int = pi/n*sum(f(cos((2*[1:n]-1)*pi/(2*n))));
err = abs(int-exact)/abs(exact);
disp(['L''errore relativo per n = ' num2str(n),' è ',num2str(err)])

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause

clear all
close all
clc
disp('***********************************************')
disp('*****************esercizio_8_5*****************')
disp('***********************************************')

clear all
close all
clc
f = @(x) (cos(x)).^2.*exp(sin(2*x));
for j = 1:2
    disp(['integrale numero ',num2str(j)])
    switch j
      case 1
           clearvars -except f
           a = 0;
           b = 2*pi;
           exact = 3.9774632605064206e+00;
      case 2
           clearvars -except f
           a = 0;
           b = 1;
           exact = 1.429777221309004;
    end
    k = 0;
    for r = 2:10
        k = k+1;
        n(k) = 2^r;
        % formula dei trapezi
        x = linspace(a,b,n(k));
        y = f(x);
        int_T(k) = (b-a)/(2*(n(k)-1))*(y(1)+2*sum(y(2:n(k)-1))+y(n(k)));
        err_T(k) = abs(exact-int_T(k))/abs(exact); 
        % formula di quadratura Gauss-Legendre
        [nodi,pesi] = gaussq(1,n(k),0,0,0,[0 0]');
        nodi_s = (b-a)/2*nodi+(b+a)/2;
        pesi_s = (b-a)/2*pesi;
        int_GL(k) = pesi_s*f(nodi_s)';
        err_GL(k) = abs(exact-int_GL(k))/abs(exact); 
    end
    fprintf('%5s\t%15s\t\t\t\t%15s\r\n','n','trapezi','Gauss-Legendre');
    for i = 1:k
        fprintf('%5.0d\t%6.16e\t%6.16e\t\n',n(i),int_T(i),int_GL(i)); 
    end
    semilogy(n,err_T,'r',n,err_GL,'b','linewidth',3)
    legend('trapezi','Gauss-Legendre')
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
disp('*****************esercizio_8_6*****************')
disp('***********************************************')

clear all
close all
clc
f = @(x,y) exp(x+y);
a = 0;
b = 2*pi;
k = 0;
for r = 2:10
    k = k+1;
    n(k) = 2^r;
    % formula dei trapezi
    t = linspace(a,b,n(k));
    % equazioni parametriche del dominio di integrazione
    x = 2*cos(t);
    y = 3*sin(t);
    jac = sqrt((-2*sin(t)).^2+(3*cos(t)).^2);
    z = f(x,y).*jac;
    int_T(k) = (b-a)/(2*(n(k)-1))*(z(1)+2*sum(z(2:n(k)-1))+z(n(k)));
end
fprintf('%5s\t%15s\r\n','n','trapezi');
for i = 1:k
    fprintf('%5.0d\t%6.16e\t\n',n(i),int_T(i)); 
end

disp('**********************************************FINE ESERCIZIO**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
