function [x, k, ier] = gauss_seidel(A, b, x, tol, kmax)
    %myFun - Gauss_seidel for sparce matrixes
    %
    % Syntax: [x, k, ier] = Gauss_Seidel(A, b, x, tol, kmax)
    %
    % Long funzione che applica l'algoritmo di gauss saidel a una matrice sparsa
    n = length(b);
    d = diag(A);
    [c, r, a] = find(A' -diag(diag(A)));
    ier = 1

    for k = 1:kmax
        y = x(1);
        x(1) = (b(1) - A(1, 2:n) * x(2:n)) / A(1, 1);
        xmax = abs(x(1));
        emax = abs(x(1) - y);

        for i = 2:n
            x(i) = (b(i) - A(i, 1:i - 1) * x(1:i - 1) - ...
                A(i, i + 1:n) * x(i + 1:n)) / A(i, i);

            if abs(x(i)) > xmax
                xmax = abs(x(i));
            end

            if abs(x(i) - y) > emax
                emax = abs(x(i) - y);
            end

        end

        if emax < tol * xmax
            ier = 0;
            return
        end

    end

end
