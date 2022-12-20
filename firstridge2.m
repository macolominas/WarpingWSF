function s = firstridge2(c,F,fmax)
%INPUT: c (orginally extracted ridge)
%       F (short-time Fourier transform)
%       fmax (maximum frequency of the STFT, in number of bins)

F = abs(F).^2;
N = length(c);
s = zeros(size(c));

K = 8; %maximum number of divisions of the original ridge

for i = 1:N
    aux = zeros(1,K); %energies for each subdivision
    for k = 1:K
        M = floor(k*fmax/c(i));
        for m = 1:M
            aux(k) = aux(k) + F(round(m*c(i)/k),i);
        end;
    end;
    [~,a] = max(aux);
    s(i) = round(c(i)/a);
end;