clear all
close all
clc

f = @(x) sqrt(x);
x = [1.8 : 0.1 : 2.2];
y = [1.341641 1.378405 1.414214 1.449138 1.483240];

z = linspace(1.8, 2.2);
a = difdiv(x, y);
p = interp(x, a, z);

plot(x, y, 'r*', z, f(z), 'r', z, p, 'b:', 'linewidth', 2);
legend('dati', 'funz', 'pl. interpolato')
xm = [1.85, 1.95, 2.05, 2.15];
pxm = interp(x, a, xm);
err = abs(f(xm) -pxm)