function d = distancias_rango(a,b,rango)
    N = length(rango);
    d = zeros(1,N);
    for i = 1:N
        d(i) = norm(a-(circshift(b,rango(i))));
    end;