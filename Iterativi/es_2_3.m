clear all
close all
clc

% A = [1 -2 2; -1 1 -1; -2 -2 1];
% b = sum(A, 2);
% n = length(b);
% x0 = zeros(n, 1);
% kmax = 100;
% toll = 1.0e-07;

% [x, k, ier] = GS_1(A, b, x0, toll, kmax)
% [xj, kj, ier_j] = jacobi(A, b, x0, toll, kmax)

%%seconda matrice
A = [4 0 2/5; 0 5 2/5; 5/2 2 1];
b = sum(A, 2);
n = length(b);
x0 = zeros(n, 1);
kmax = 100;
toll = 1.0e-07;

[x, k, ier] = GS_1(A, b, x0, toll, kmax)
[xj, kj, ier_j] = jacobi(A, b, x0, toll, kmax)
