function y = trigon_reg(x,t,order)
%t = 0:1/length(x):1-1/length(x);
% t = 1:length(x);
modelfun = @(b,n)(sum(diag(b)*([cos(2*pi*((1:floor(length(b)/2))'*n));sin(2*pi*((1:floor(length(b)/2))'*n))])));
b_e = nlinfit(t,x,modelfun,ones(1,2*order)+0.1*randn(1,2*order));
y = sum( diag(b_e)*([cos(2*pi*(1:order)'*t);sin(2*pi*(1:order)'*t)]));  
