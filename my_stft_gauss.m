function [F sumaf] = my_stft_gauss(x,redun,sigma,fmax)
% FUNCTION to compute STFT with Gaussian window
% INPUT:    x: signal to be analyzed
%       redun: resolution of the STFT. redun = 1 computes the STFT with N/2
%       frequency bins for a signal of length N when fmax = 0.5.
%       sigma: parameter of the Gaussian window
%        fmax: maximum normalized frequency (usually 0.5)
% OUTPUT:   F: Short-Time Fourier Transform
%       sumaf: accumulation of the moduli of the filters in frequency
%       domain (useful for reconstruction).
N = length(x);
F = zeros(round(N*redun*fmax),N);
X = fft(x);
sumaf = zeros(size(X));

f = [(0:1/N:0.5) (-0.5+1/N:1/N:-1/N)];%nomalized frequencies

for h = 1/(N*redun) : 1/(N*redun) : fmax
    G = sqrt(pi/sigma)*exp(-pi^2*(f-h).^2/sigma);%filter
    F(round(h*N*redun),:) = (ifft(X.*(G)));
    sumaf = sumaf + abs(G);
end;

