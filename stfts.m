function [Fg,Ftg,Ftgp,Ft2g,Fgp] = stfts(s,sigma,ft,bt)
% stfts : computes 5 STFTs of a signal with different windows.
%
% INPUTS:
%   s: real or complex signal, must be of length 2^N
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
%
% OUTPUTS:
%   Fg: the short-time Fourier transform (STFT) with a Gaussian window g.
%   Ftg: STFT with tg.
%   Ftgp: STFT with tg'.
%   Ft2g: STFT with t^2g.
%   Fgp: STFT with g'.


n = length(s);
%s = s(:);

% Optional parameters
if nargin<5
    ft = 1:n/2;
    bt = 1:n;
end
nb = length(bt);
neta = length(ft);

S = fftshift(fft(s));
f = -n/2:n/2-1;

% Initialization
Fg = zeros(neta,nb);
Ftg = zeros(neta,nb);
Ftgp = zeros(neta,nb);
Ft2g = zeros(neta,nb);
Fgp = zeros(neta,nb);

for b = 1:length(ft) % For each row
    
    % Computes STFTs
    i = 0;
    g_hat = (sigma*sqrt(pi)/(1i*2*pi))^i*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
    tmp = ifft(S.*conj(g_hat));
    Fg(b,:) = tmp;
    
    i = 1;
    g_hat = (sigma*sqrt(pi)/(1i*2*pi))^i*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
    tmp = ifft(S.*conj(g_hat));
    Ftg(b,:) = tmp;
    
    i = 2;
    g_hat = (sigma*sqrt(pi)/(1i*2*pi))^i*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
    tmp = ifft(S.*conj(g_hat));
    Ft2g(b,:) = tmp;
   
    i = 0;
    tgp_hat = (-2*pi/sigma^2)*(sigma*sqrt(pi)/(1i*2*pi))^(i+1)*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i+1).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
    tmp = ifft(S.*conj(tgp_hat));
    Fgp(b,:) = tmp;
    
    i = 1;
    tgp_hat = (-2*pi/sigma^2)*(sigma*sqrt(pi)/(1i*2*pi))^(i+1)*...
            hermite_poly(sigma*sqrt(pi)*(f-ft(b)),i+1).*...
            exp(-sigma^2*pi*(f-ft(b)).^2);
    tmp = ifft(S.*conj(tgp_hat));
    Ftgp(b,:) = tmp;

end
end

function hp = hermite_poly(t,order)

h{1} = @(x)ones(1,length(x));
h{2} = @(x)2*x;
h{3} = @(x)4*x.^2   - 2;
h{4} = @(x)8*x.^3   - 12*x;
h{5} = @(x)16*x.^4  - 48*x.^2   + 12;
h{6} = @(x)32*x.^5  - 160*x.^3  + 120*x;
h{7} = @(x)64*x.^6  - 480*x.^4  + 720*x.^2      - 120;
h{8} = @(x)128*x.^7 - 1344*x.^5 + 3360*x.^3     - 1680*x;
h{9} = @(x)256*x.^8 - 3584*x.^6 + 13440*x.^4    - 13440*x.^2 + 1680;
h{10}= @(x)512*x.^9 - 9216*x.^7 + 48384*x.^5    - 80640*x.^3 + 30240*x;

hp = h{order+1}(t);

end