function [x, k, ier] = jacobi(A, b, x0, toll, kmax)
    ier = 1;
    x = [];
    k = 0;
    n = length(b);
    D = diag(diag(A));
    rho = max(abs(eig(eye(n)-inv(D)*A)));
    if rho > 1
        ier = 2;
        return
    C = A-D;
    for K = 1:kmax
        x = D\(b-C*x0);
        if norm(x - x0, inf) <= toll*norm(x, inf)
            ier = 0;
            break
        end
        x0 = x;
    end
end