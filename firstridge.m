function s = firstridge(F)
%INPUT: F (short-time Fourier transform)

N = length(F(1,:));
gamma = median(abs(real(F(:))));
F(abs(F)<gamma) = 0;
F = abs(F).^2;
s = zeros(1,N);

for i = 1:N
    [~,s(i)] = findpeaks(F(:,i),'npeaks',1);
end;