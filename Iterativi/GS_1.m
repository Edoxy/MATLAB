function [x, k, ier] = GS_1(A, b, x0, toll, kmax)
    ier = 1;
    x = [];
    k = 0;
    n = length(b);
    D = tril(A);
    %controlla della convergenza
    rho = max(abs(eig(eye(n)-inv(D)*A)))
    if rho > 1
        %se non converge la matrice allora interrompe il metodo
        ier = 2;
        return
    end
    C = A-D;
    for k = 1:kmax
        x = D\(b-C*x0);
        if norm(x - x0, inf) <= toll*norm(x, inf)
            ier = 0;
            break
        end
        x0 = x;
    end
end