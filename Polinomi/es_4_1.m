clear all
close all
clc

for n = 2 : 20
    x = linspace(0, 1, n);
    V = vander(x);
    norm1(n-1) = cond(V, 1);
    norm2(n-1) = cond(V, 2);
    norminf(n-1) = cond(V, inf);

    b = sum(V,2);
    x = V\b;
    err1(n-1) = norm(ones(n-1) - x, 1) /norm(ones(n-1), 1);
    err2(n-1) = norm(ones(n-1) - x, 2) /norm(ones(n-1), 2);
    errinf(n-1) = norm(ones(n-1) - x, inf) /norm(ones(n-1), inf);
end

figure(1);
semilogy([2:20], norm1, 'r*', [2:20], norm2, 'bo', [2:20], norminf, 'kd');
figure(2);
semilogy([2:20], err1, 'r*',[2:20], err1, 'bo',[2:20], err1, 'kd' );
