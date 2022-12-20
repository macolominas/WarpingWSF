function [frec,A] = frequency_and_amplitude(F,omega,c)
[~,N] = size(F);
frec = zeros(1,N);
A = zeros(1,N);
    for n = 1:N
        frec(n) = omega(c(n),n); %guardamos el valor del operador frecuencial de synchro
        A(n) = abs(F(c(n),n)); %estimamos la amplitud
    end;
A = 2*pi*A;%correcion normalizando por el máximo del filtro