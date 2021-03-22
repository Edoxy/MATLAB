function [x, k, ier] = gauss_seidel_1(A, b, x, toll, kmax)
    n = length(b);
    ier = 1;
    d = diag(A);
    A_ = A - diag(d);

    for k = 1:kmax
        xmax = 0;
        emax = 0;

        for j = 1:n
            y = x(j);
            x(j) = (b(j) - sum(A_(j, :) * x, 2)) / d(j);

            if abs(x(j)) > xmax
                xmax = abs(x(j));
            end

            if abs(y - x(j)) > emax
                emax = abs(y - x(j));
            end

        end

        if emax < toll * xmax
            ier = 0;
            break;
        end

    end