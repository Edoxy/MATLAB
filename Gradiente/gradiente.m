function [x, k, ier] = gradiente(A, b, x, toll, kmax)
    ier = 1;
    r = b-A*x;
    for k =1:kmax
        %operazione più costosa
        s = A*r;
        alpha = r'*r/(r'*s);
        x = x + alpha * r; %x^(k+1) = x^(k) + lpha * r ===> alpha-
        %r = r - alpha * s; %ogni tanto meglio aggiornare r = b - A * x
        r = b - A * x;
        if norm(r) <= toll*norm(b)
            ier = 0;
            return
        end
    end
end