%come horner ma parte dalla rappresentazione di newton del polinomio e non quella monomiale
function p = interp(x, a, z)
    n= length(x) -1;
    p =a(n+1) * ones(size(z));%un vettore p le cui componenti sono tutte uguali a a(n+1)
    for i = n:-1:1
        p = p .* (z-x(i)) + a(i);
    end
end