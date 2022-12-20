function [y,t_m,m,t_m_] = demodulation(x,t,sigma,redun,fmax,jump,flag)
N = length(x);
bt = 1:N;
ft = 1/redun:1/redun:round(fmax*N)-1/redun; 
gamma = 0;
[Fx,~,~,~,~,~,omega2] = sstn_test_modL(x,gamma,sigma,ft,bt);
c = exridge(Fx,0,0,jump);
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = 2*pi*abs(Fx(c(i),i)); 
end;
phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;
t_m_ = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m_,'spline');
a = interp1(t,1:N,0);
M = interp1(m,t_m_,a);
t_m = t_m_ - M;

if flag == 1
    y = interp1(1:N,x./A,m,'spline');
else
    y = interp1(1:N,x,m,'spline');
end