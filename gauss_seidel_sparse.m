function [x, k, ier] = gauss_seidel(A, b, x, tol, kmax)
    %myFun - Gauss_seidel for sparce matrixes
    %
    % Syntax: [x, k, ier] = Gauss_Seidel(A, b, x, tol, kmax)
    %
    % Long funzione che applica l'algoritmo di gauss saidel a una matrice sparsa
    %1:00 h nella lezione del 12
    n = length(b);
    d = diag(A);
    [c, r, a] = find(A' -diag(d));
    ier = 1;

    for k = 1:kmax

        xmax = 0;
        emax = 0;
        count = 0;

        for i = 1:n
            y = x(i);

            count0 = count;
            
            while count +1 <= length(r) && r(count +1) == i
                count = count + 1;
            end

            s = a(count0 +1:count)'*x(c(count0 + 1:count));
            x(i) = (b(i) - s) / d(i);

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
