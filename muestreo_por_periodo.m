t_r = 100+rand(1,21);
N = 100;
t_r = 110+ [0:1/N:1-1/N];
y_r = interp1(0.5*t_m/pi,y,t_r);
figure;scatter(mod(t_r,1),y_r)
hold on
[a,b] = sort(mod(t_r,1));
plot(a,trigon_reg(y_r(b),a,11))