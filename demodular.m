function [y,M] = demodular(x,sigma,redun,fmax,jump)
N = length(x);
bt = 1:N; 
ft = 1/redun:1/redun:round(fmax*N)-1/redun; 
gamma = 0; 
[Fx,~,~,~,~,~,omega2] = sstn_test_modL(x,gamma,sigma,ft,bt);
c = exridge(Fx,0,0,jump);
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fx(c(i),i)); 
end;
phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;
t_m = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m,'spline');
y = interp1(1:N,x./A,m,'spline');
P = ((t_m(end)-t_m(1)+1)/(2*pi));
T = floor(N/P); 
P = floor(P);
M = zeros(P,T);
for i = 1:P
    M(i,:) = interp1(0.5*t_m/pi,y,i-1+[0:1/(T):1-1/(T)],'spline');
end;

