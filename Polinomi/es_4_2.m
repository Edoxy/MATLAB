clear all
close all
clc

%f = @(x) sin(pi*x);
f = @(x) 1./(1 +x.^2);
a = -5;
b = 5;
z = linspace(a, b);
for n = 5:5:80
    n
    %x = linspace(a, b, n);
    i = 1:n+1;
    zz = -cos((2 *i -1)*pi/(2*(n+1)));
    x = (b-a)/2*zz+(b+a)/2
    y = f(x);
    c = polyfit(x, y, n);
    p = polyval(c, z);
    plot(x, y, 'r*', z, f(z), 'r', z, p, 'b:')
    err = norm(f(x) - p)
    pause
end