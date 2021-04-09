%interpolazione del polinomio tramite la scrittura di newton
function a = difdiv(x, y)
    n = length(x) -1;
    for i = 1 : n
        for j = n+1:-1:i+1
            y(j) = (y(j) - y(j-1))/(x(j) - x(j-i));
        end
        a = y;
    end
end
%potremmo scrivere questo algoritmo con un solo ciclo for e non due innestati
