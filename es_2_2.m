clear all
close all
clc

n = 400
A = sprandsym(n, 0.02, rand(n, 1) * 100);
x1 = zeros(n, 1);
b = sum(A, 2);
toll = 1.0e-06;
kmax = 10000;

tic
[x, k, ier] = gauss_seidel_sparse(A, b, x1, toll, kmax);
t_GSS = toc
tic
[x, k, ier] = gauss_seidel(A, b, x, toll, kmax);
t_GS = toc