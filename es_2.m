A = hilb(3);
b = sum(A, 2);
x = zeros(3, 1);
toll = 1.0e-03;
kmax = 1000;
[x, k, ier] = jacobi(A, b, x, toll, kmax)
%min 1:18h lezione 19-3