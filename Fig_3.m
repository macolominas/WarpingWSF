N = 6000;
t = -0.5:2/N:1.5-2/N;

phi = 2*pi*40*t + cos(2*pi*4*t);

a = cos(phi);
b = (1./(1 + exp(-250*(t-0.33))) - 1./(1 + exp(-250*(t-0.66)))).*cos(2*phi);
c = (1./(1 + exp(-250*(t-0.33)))).*cos(3*phi);
x = a+b+c;

t_aux = 0:0.001:1-0.001;
wsf1 = cos(2*pi*t_aux+1);
wsf2 = cos(2*pi*(t_aux)+1) + cos(2*pi*2*(t_aux)+2) + cos(2*pi*3*(t_aux)+3);
wsf3 = cos(2*pi*t_aux+1) + cos(2*pi*3*t_aux+3);


amp_ruido = [0.0316];% .1 .3162];
J = 100; %cantidad de realizaciones
rango = -50:50; %rango de corrimientos

redun = 1;
fmax = 0.025;
sigma = 0.05;
bt = 1:N;
ft = 1/redun:1/redun:round(fmax*N)-1/redun;
gamma = 0; %umbral para la STFT

jump = 10*redun;

K = 3;%number of clusters
    
error1 = zeros(3,3,J);
error2 = zeros(3,3,J);
error3 = zeros(3,3,J);

svd_entropy = zeros(3,3,J);

cant_ciclos = 39;

for i = 1:3
    for j = 1:J
        
        ff = cumsum(randn(1,N));
    IFr = ff./max(abs(ff));% + pi/2;
    IFr = IFr - min(IFr);
    IF = smooth(IFr,80);
    
    aux = (cumsum(IF)')/N;

    phin = phi + aux - aux(501);
%     phin = phi;
    a = cos(phin);
    b = (1./(1 + exp(-250*(t-0.33))) - 1./(1 + exp(-250*(t-0.66)))).*cos(2*phin);
    c = (1./(1 + exp(-250*(t-0.33)))).*cos(3*phin);
    x = a+b+c;

    xn = x + amp_ruido(i)*std(x)*randn(size(x));    
    
    [y1,t_m1,m1] = demodulation(xn,t,sigma,redun,fmax,jump,0);
    
    M1 = segmentation(y1,t_m1,cant_ciclos,1000);
    
%---------------------

    I = kmeans(M1,K,'replicates',50);
    I = sort(I);
    WSF = zeros(K,1000); %formas de onda encontradas
    for k = 1:K
        WSF(k,:) = median(M1(I==k,:),1);%podemos usar mean, median, etc.
    end;

    error1(1,i,j) = min(distancias_rango(wsf1,WSF(1,:),rango));
    error2(1,i,j) = min(distancias_rango(wsf2,WSF(2,:),rango));
    error3(1,i,j) = min(distancias_rango(wsf3,WSF(3,:),rango));
    
    [~,s,~] = svd(M1);
    s = diag(s);
    s = s/sum(s);
    svd_entropy(1,i,j) = -sum(s.*log(s));
% %--------------------------------

    [y2,t_m2,m2] = demodulation(y1,.5*t_m1/pi,sigma,redun,fmax,jump,0);
    
    M2 = segmentation(y2,t_m2,cant_ciclos,1000);

%-------------------
    I = kmeans(M2,K,'replicates',50);
    I = sort(I);
    WSF = zeros(K,1000); %formas de onda encontradas
    for k = 1:K
        WSF(k,:) = median(M2(I==k,:),1);%podemos usar mean, median, etc.
    end;

    error1(2,i,j) = min(distancias_rango(wsf1,WSF(1,:),rango));
    error2(2,i,j) = min(distancias_rango(wsf2,WSF(2,:),rango));
    error3(2,i,j) = min(distancias_rango(wsf3,WSF(3,:),rango));
          
    [~,s,~] = svd(M2);
    s = diag(s);
    s = s/sum(s);
    svd_entropy(2,i,j) = -sum(s.*log(s));
% %--------------------------------

    [y3,t_m3,m3] = demodulation(y2,.5*t_m2/pi,sigma,redun,fmax,jump,0);
    
    M3 = segmentation(y3,t_m3,cant_ciclos,1000);


%---------------------
    I = kmeans(M3,K,'replicates',50);
    I = sort(I);
    WSF = zeros(K,1000); %formas de onda encontradas
    for k = 1:K
        WSF(k,:) = median(M3(I==k,:),1);%podemos usar mean, median, etc.
    end;

    error1(3,i,j) = min(distancias_rango(wsf1,WSF(1,:),rango));
    error2(3,i,j) = min(distancias_rango(wsf2,WSF(2,:),rango));
    error3(3,i,j) = min(distancias_rango(wsf3,WSF(3,:),rango));
    
    [~,s,~] = svd(M3);
    s = diag(s);
    s = s/sum(s);
    svd_entropy(3,i,j) = -sum(s.*log(s));
%--------------------------------
i
j
    end;    
end;

figure;
subplot(4,3,[1 4])
boxplot(squeeze(error1(:,1,:))')
title('SNR = 30 dB')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,1.25,'waveform no. 1','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,[2 5])
boxplot(squeeze(error1(:,2,:))')
title('SNR = 20 dB')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,1.9,'waveform no. 1','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,[3 6])
boxplot(squeeze(error1(:,3,:))')
title('SNR = 10 dB')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,3.7,'waveform no. 1','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,7)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,1,1,:)));
    hold on;
end;
plot(t_aux,wsf(1,:),'k--')
xlabel('time')
legend('1 iter','2 iter','3 iter','GT')
text(0.1,0.8,'waveform no. 1','fontsize',20)
set(gca,'fontsize',20)

subplot(4,3,8)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,2,1,:)));
    hold on;
end;
plot(t_aux,wsf(1,:),'k--')
xlabel('time')
text(0.1,0.8,'waveform no. 1','fontsize',20)
set(gca,'fontsize',20)

subplot(4,3,9)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,3,1,:)));
    hold on;
end;
plot(t_aux,wsf(1,:),'k--')
xlabel('time')
text(0.1,0.8,'waveform no. 1','fontsize',20)
set(gca,'fontsize',20)




figure;
subplot(4,3,[1 4])
boxplot(squeeze(error2(:,1,:))')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,14,'waveform no. 2','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,[2 5])
boxplot(squeeze(error2(:,2,:))')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,14,'waveform no. 2','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,[3 6])
boxplot(squeeze(error2(:,3,:))')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,17,'waveform no. 2','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,7)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,1,2,:)));
    hold on;
end;
plot(t_aux,wsf(2,:),'k--')
xlabel('time')
text(0.1,3.2,'waveform no. 2','fontsize',20)
set(gca,'fontsize',20)

subplot(4,3,8)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,2,2,:)));
    hold on;
end;
plot(t_aux,wsf(2,:),'k--')
xlabel('time')
text(0.1,3.2,'waveform no. 2','fontsize',20)
set(gca,'fontsize',20)

subplot(4,3,9)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,3,2,:)));
    hold on;
end;
plot(t_aux,wsf(2,:),'k--')
xlabel('time')
text(0.1,3.2,'waveform no. 2','fontsize',20)
set(gca,'fontsize',20)



figure;
subplot(4,3,[1 4])
boxplot(squeeze(error3(:,1,:))')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,5.5,'waveform no. 3','fontsize',20)
grid on
set(gca,'fontsize',20)

subplot(4,3,[2 5])
boxplot(squeeze(error3(:,2,:))')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,5.5,'waveform no. 3','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,[3 6])
boxplot(squeeze(error3(:,3,:))')
ylabel('waveform estimation error')
xlabel('iterations')
text(2.5,8.5,'waveform no. 3','fontsize',20)
set(gca,'fontsize',20)
grid on

subplot(4,3,7)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,1,3,:)));
    hold on;
end;
plot(t_aux,wsf(3,:),'k--')
xlabel('time')
text(0.1,1.6,'waveform no. 3','fontsize',20)
set(gca,'fontsize',20)

subplot(4,3,8)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,2,3,:)));
    hold on;
end;
plot(t_aux,wsf(3,:),'k--')
xlabel('time')
text(0.1,1.6,'waveform no. 3','fontsize',20)
set(gca,'fontsize',20)

subplot(4,3,9)
set(gca,'fontsize',9.5)
for iter = 1:3
    plot(t_aux,squeeze(WSF(iter,3,3,:)));
    hold on;
end;
plot(t_aux,wsf(3,:),'k--')
xlabel('time')
text(0.1,1.6,'waveform no. 3','fontsize',20)
set(gca,'fontsize',20)


figure;
subplot(4,3,[1 4])
boxplot(squeeze(svd_entropy(:,1,:))')
ylabel('matrix SVD entropy')
xlabel('iterations')
set(gca,'fontsize',20)
grid on

subplot(4,3,[2 5])
boxplot(squeeze(svd_entropy(:,2,:))')
ylabel('matrix SVD entropy')
xlabel('iterations')
set(gca,'fontsize',20)
grid on

subplot(4,3,[3 6])
boxplot(squeeze(svd_entropy(:,3,:))')
ylabel('matrix SVD entropy')
xlabel('iterations')
set(gca,'fontsize',20)
grid on
