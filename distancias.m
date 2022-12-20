function d = distancias(a,b)
    N = length(a);
    d = zeros(1,N);
    d(1) = norm(a-b);
    for i = 1:N-1
        d(i+1) = norm(a-(circshift(b,i)));
    end;