N = 6000;
t = 0:1/N:1-1/N;
rng(1)

phi = 2*pi*20*t + 10*cos(2*pi*2*t);
phi = 2*pi*50*t + 2*pi*20*t.^2;
phi = 2*pi*40*t + cos(2*pi*4*t);

a = cos(phi);
b = (1./(1 + exp(-250*(t-0.33))) - 1./(1 + exp(-250*(t-0.66)))).*cos(2*phi);
c = (1./(1 + exp(-250*(t-0.33)))).*cos(3*phi);
x = a+b+c;
xn = x + 0.5*std(x)*randn(size(x));

%-----STFT
N = length(x);
redun = 1;
fmax = 0.025;
sigma = 0.05;
bt = 1:N; %vector temporal
ft = 1/redun:1/redun:round(fmax*N)-1/redun; %vector frecuencial
gamma = 0; %umbral para la STFT
[Fxn,~,~,~,~,~,omega2] = sstn_test_modL(xn,gamma,sigma,ft,bt);
%------------

c = exridge(Fxn,0,0,10*redun);%ridge extraction
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fxn(c(i),i));
end;

phi_est = 2*pi*cumsum(frec)/N;

%-----resampling
sobremuestreo = 1;
t_m = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m,'spline');
y = interp1(1:N,xn,m,'spline');
[Fy,~,~,~,~,~,omega2] = sstn_test_modL(y,gamma,sigma,ft,bt);
%-------------

%-----segmentation
P = ((t_m(end)-t_m(1)+1)/(2*pi));
T = floor(N/P); 
P = floor(P);
M = zeros(P,T); 

for i = 1:P
    M(i,:) = interp1(0.5*(t_m)/pi,y,i-1+[1/T:1/(T):1],'spline');
end;
%---------------------

%----clusterizado y estimacion de formas de onda------
eva = evalclusters(M,'kmeans','CalinskiHarabasz','KList',[1:6]);
% K = eva.OptimalK;
K = 3;
I = kmeans(M,K,'replicates',50);
I = sort(I);
WSF = zeros(K,T); 
regWSF = zeros(K,T);
orden = 3; 
    WSF(i,:) = median(M(I==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T:1-1/T,orden);
end;
%------------------------------------------------------

%-----figures----------------------------------

figure;
subplot(2,3,1:3)
plot(t,xn,'k'); set(gca,'yticklabel',[]); text(0.025,5,'Noisy artificial signal')
hold on; plot(t,x,'r')
xlabel('time')
subplot(2,3,4);plot(0:1/T:1-1/T,regWSF(1,:),'k'); set(gca,'yticklabel',[]); xlabel('time'); text(0.5,0.8,'waveform no. 1')
%hold on; plot(0:2*pi/T:2*pi-2*pi/T,cos(t_m(1:T)+phi(1)-4*phi_est(1)),'r--')
hold on; plot(0:1/T:1-1/T,cos(t_m(1:T)+1-t_m(1)),'r--')
% plot(0:1/T:1-1/T,trigon_reg(0:1/T:1-1/T,1),'r--')
subplot(2,3,5);plot(0:1/T:1-1/T,regWSF(2,:),'k'); set(gca,'yticklabel',[]); xlabel('time'); text(0.5,2.4,'waveform no. 2')
%hold on; plot(0:2*pi/T:2*pi-2*pi/T,cos(t_m(1:T)+phi(1)-6*phi_est(1))+cos(2*(t_m(1:T)+phi(1)-6*phi_est(1)))+cos(3*(t_m(1:T)+phi(1)-6*phi_est(1))),'r--')
hold on; plot(0:1/T:1-1/T,cos(t_m(1:T)+1-t_m(1))+cos(2*(t_m(1:T))+2-2*t_m(1))+cos(3*(t_m(1:T))+3-3*t_m(1)),'r--')
subplot(2,3,6);plot(0:1/T:1-1/T,regWSF(3,:),'k'); set(gca,'yticklabel',[]); xlabel('time'); text(0.5,1.6,'waveform no. 3');
%hold on; plot(0:2*pi/T:2*pi-2*pi/T,cos(t_m(1:T)+phi(1)-4*phi_est(1))+cos(3*(t_m(1:T)+phi(1)-4*phi_est(1))),'r--')
hold on; plot(0:1/T:1-1/T,cos(t_m(1:T)+1-t_m(1))+cos(3*(t_m(1:T))+3-3*t_m(1)),'r--')

