load('213m.mat');
x = val(1,360*60*10+38:360*60*10.5);

%---high-pass filtering--
X = fft(x);
X(1:20) = 0;
X(end-19:end) = 0;
x = real(ifft(X));
%--------------------

N = length(x);
fs = 360;
redun = 4;
fmax = 0.025;
sigma = 0.05;%
bt = 1:N;
ft = 1/redun:1/redun:round(fmax*N)-1/redun;
gamma = 0; 
[Fx,~,~,~,~,~,omega2] = sstn_test_modL(x,gamma,sigma,ft,bt);
Fx(1:5*redun,:) = 0; 
c1 = exridge(Fx(1:320,:),0,0,10*redun);
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c1(i),i); 
    A(i) = abs(Fx(c1(i),i));
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;

%---
sobremuestreo = 1;
t_m1 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m1,'spline');
y1 = interp1(1:N,x./A,m,'spline');
[Fy1,~,~,~,~,~,~] = sstn_test_modL(y1,gamma,sigma,ft,bt);
%-----------------------

P = ((t_m1(end)-t_m1(1))/(2*pi));
T1 = floor(N/P); 
P = floor(P);
M1 = zeros(P,T1); 

for i = 1:P
    M1(i,:) = interp1(0.5*t_m1/pi,y1,i-1+[0:1/(T1):1-1/(T1)],'spline');
end;
%---------------------

c2 = exridge(Fx(1:520,:),0,0,10*redun);
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c2(i),i); 
    A(i) = abs(Fx(c2(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;

%---
sobremuestreo = 1;
t_m2 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m2,'spline');
y2 = interp1(1:N,x./A,m,'spline');
[Fy2,~,~,~,~,~,~] = sstn_test_modL(y2,gamma,sigma,ft,bt);
%-----------------------

%-------
c2 = exridge(Fy2,0,0,10*redun);
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c2(i),i); 
    A(i) = abs(Fx(c2(i),i)); 
end;
A = 2*pi*A;
phi_est = 2*pi*cumsum(frec)/N;

sobremuestreo = 1;
t_m2 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m2,'spline');
y2 = interp1(1:N,x./A,m,'spline');

[Fy2,~,~,~,~,~,~] = sstn_test_modL(y2,gamma,sigma,ft,bt);

%--------------------------

%---
t_m2 = t_m2/2;
P = ((t_m2(end)-t_m2(1))/(2*pi));
T2 = floor(N/P);
P = floor(P);
M2 = zeros(P,T2);

for i = 1:P
    M2(i,:) = interp1(0.5*t_m2/pi,y2,i-1+[0:1/(T2):1-1/(T2)],'spline');
end;
%---------------------

sobremuestreo = 1;%
t_m3 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m3,'spline');
y3 = interp1(1:N,x./A,m,'spline');
[Fy3,~,~,~,~,~,~] = sstn_test_modL(y3,gamma,sigma,ft,bt);
%-----------------------

%--
t_m3 = t_m3/1.991;%
P = ((t_m3(end)-t_m3(1))/(2*pi));%
T3 = floor(N/P); 
P = floor(P);
M3 = zeros(P,T3);
for i = 1:P
    M3(i,:) = interp1(0.5*t_m3/pi,y3,i-1+[0:1/(T3):1-1/(T3)],'spline');
end;
%---------------------



%-----figures----------------------------------
t = 0:1/fs:N/fs-1/fs;
f = fs/(N*redun):fs/(N*redun):fs*fmax;
aux = 0.5*diff(t_m1)/pi;
f_m1 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
figure;
subplot(10,4,1); plot(t,x,'k'); set(gca,'xticklabel',[]); xlim([t(1) t(end)]); title('original ECG signal');% text(1,-380,'ECG signal')
subplot(10,4,[5 9 13 17]); imagesc(t,f,abs(Fx)); colormap(1-gray); set(gca,'ydir','normal')
hold on; plot(t,fs*c1/(N*redun),'r--'); plot(t,fs*c2/(N*redun),'b--') 
xlabel('time [s]')
ylabel('frequency [Hz]')

subplot(10,4,2); plot(0.5*t_m1/pi,y1,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m1(1)/pi 0.5*t_m1(end)/pi]);
title('after 1 demodulation w/ 1st harm. (red)')% text(2,-380,'warped ECG signal')
subplot(10,4,[6 10 14 18]); imagesc(0.5*t_m1/pi,f_m1,abs(Fy1)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

aux = 0.5*diff(t_m2)/pi;
f_m2 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,4,3); plot(0.5*t_m2/pi,y2,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m2(1)/pi 0.5*t_m2(end)/pi]);
title('after 1 demodulation w/ 2nd harm. (blue)');% text(1,-380,'ECG signal')
subplot(10,4,[7 11 15 19]); imagesc(0.5*t_m2/pi,f_m2,abs(Fy2)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,4,4); plot(0.5*t_m2/pi,y2,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m2(1)/pi 0.5*t_m2(end)/pi]);
title('after 1 demodulation w/ 2nd harm. (blue)');% text(1,-380,'ECG signal')
subplot(10,4,[8 12 16 20]); imagesc(0.5*t_m2/pi,f_m2,abs(Fy2)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')


subplot(10,4,[26 30 34 38]); imagesc(M1); colormap(1-gray)
title({'Waveform Matrix','SVD entropy = 3.0230'})
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

subplot(10,4,[27 31 35 39]); imagesc(M2); colormap(1-gray)
title({'Waveform Matrix (dividing by 2)','SVD entropy = 2.7013'})
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

subplot(10,4,[28 32 36 40]); imagesc(M3); colormap(1-gray)
title({'Waveform Matrix (dividing by 1.991)','SVD entropy = 2.2299'})
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')



t = 0:1/fs:N/fs-1/fs;
f = fs/(N*redun):fs/(N*redun):fs*fmax;
aux = 0.5*diff(t_m1)/pi;
f_m1 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
figure;
subplot(10,15,1:5); plot(t,x,'k'); set(gca,'xticklabel',[]); xlim([t(1) t(end)]); title('original ECG signal');% text(1,-380,'ECG signal')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(t,f,abs(Fx)); colormap(1-gray); set(gca,'ydir','normal')
hold on; plot(t,fs*c1/(N*redun),'r--');
xlabel('time [s]')
ylabel('frequency [Hz]')
subplot(10,15,76:80); plot(0.5*t_m1/pi,y1,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m1(1)/pi 0.5*t_m1(end)/pi]);
title('after 1 demodulation w/ 1st harm.')% text(2,-380,'warped ECG signal')
subplot(10,15,[91:95 106:110 121:125 136:140]); imagesc(0.5*t_m1/pi,f_m1,abs(Fy1)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[6:15 21:30]); plot(t,x,'k'); xlim([0 15]); title('zoomed in original ECG signal'); xlabel('time [s]')
subplot(10,15,[36:45 51:60]); plot(t,x,'k'); xlim([15 30]); xlabel('time [s]')


subplot(10,15,[81:85 96:100 111:115 126:130 141:145]);
imagesc(M1); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M1))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%--------
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',50));
eva = evalclusters(M1,myfunc,'CalinskiHarabasz','KList',[1:6]);
K1 = eva.OptimalK;
I1 = kmeans(M1,K1,'replicates',50);
WSF = zeros(K1,T1); 
regWSF = zeros(K1,T1);
orden = 7; 
for i = 1:K1
    WSF(i,:) = median(M1(I1==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T1:1-1/T1,orden);
end;
%------------------------------------------------------

subplot(10,15,[86 101 116 131 146]); plot(I1,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K1,5)
    subplot(10,15,[15*(i+5)-3:15*(i+5)]); plot(0:2*pi/T1:2*pi-2*pi/T1,WSF(i,:),'k')
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;

figure;
aux = 0.5*diff(t_m2)/pi;
f_m2 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,15,1:5); plot(0.5*t_m2/pi,y2,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m2(1)/pi 0.5*t_m2(end)/pi]);
title('after 1 demodulation w/ 2nd harm.');% text(1,-380,'ECG signal')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(0.5*t_m2/pi,f_m2,abs(Fy2)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[6:10 21:25 36:40 51:55 66:70]); imagesc(M2); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M2))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%------
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',50));
eva = evalclusters(M2,myfunc,'CalinskiHarabasz','KList',[1:6]);
K2 = eva.OptimalK;
I2 = kmeans(M2,K2,'replicates',50);
WSF = zeros(K2,T2); 
regWSF = zeros(K2,T2);
orden = 7; 
for i = 1:K2
    WSF(i,:) = median(M2(I2==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T2:1-1/T2,orden);
end;
%------------------------------------------------------

subplot(10,15,[11 26 41 56 71]); plot(I2,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K2,5)
    subplot(10,15,[15*(i)-3:15*(i)]); plot(0:2*pi/T2:2*pi-2*pi/T2,WSF(i,:),'k')
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;

aux = 0.5*diff(t_m3)/pi;
f_m3 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,15,[76:80]); plot(0.5*t_m3/pi,y3,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m3(1)/pi 0.5*t_m3(end)/pi]);
title('after 1 demodulation w/ 2nd harm. (corrected)');% text(1,-380,'ECG signal')
subplot(10,15,[91:95 106:110 121:125 136:140]); imagesc(0.5*t_m3/pi,f_m3,abs(Fy3)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[81:85 96:100 111:115 126:130 141:145]); imagesc(M3); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M3))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%-----
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',50));
eva = evalclusters(M3,myfunc,'CalinskiHarabasz','KList',[1:6]);
K3 = eva.OptimalK;
I3 = kmeans(M3,K3,'replicates',50);
WSF = zeros(K3,T3); 
regWSF = zeros(K3,T3);
orden = 7; 
for i = 1:K3
    WSF(i,:) = median(M3(I3==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T3:1-1/T3,orden);
end;
%------------------------------------------------------

subplot(10,15,[86 101 116 131 146]); plot(I3,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K3,5)
    subplot(10,15,[15*(i+5)-3:15*(i+5)]); plot(0:2*pi/T3:2*pi-2*pi/T3,WSF(i,:),'k')
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;



figure;
M3s = synchronization(M3);
aux = 0.5*diff(t_m2)/pi;
f_m2 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,15,1:5); plot(0.5*t_m2/pi,y2,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m2(1)/pi 0.5*t_m2(end)/pi]);
title('after 1 demodulation w/ 2nd harm. synchronized');% text(1,-380,'ECG signal')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(0.5*t_m2/pi,f_m2,abs(Fy2)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[6:10 21:25 36:40 51:55 66:70]); imagesc(M3s); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M3s))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%-----
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',50));
eva = evalclusters(M3s,myfunc,'CalinskiHarabasz','KList',[1:6]);
K3s = eva.OptimalK;
% K = 3;
I3s = kmeans(M3s,K3s,'replicates',50);
WSF = zeros(K3s,T3); 
regWSF = zeros(K3s,T3);
orden = 7; 
for i = 1:K3s
    WSF(i,:) = median(M3s(I3s==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T3:1-1/T3,orden);
end;
%------------------------------------------------------

subplot(10,15,[11 26 41 56 71]); plot(I3s,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K3s,5)
    subplot(10,15,[15*(i)-3:15*(i)]);
    plot(0:2*pi/T3:2*pi-2*pi/T3,transpose(M3s(I3s==i,:)))
    hold on
    plot(0:2*pi/T3:2*pi-2*pi/T3,WSF(i,:),'k','linewidth',2)
    
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;


figure;
subplot(311)
for i = 1:length(I3s)
    plot(i-1+[0:1/T3:1-1/T3],M3(i,:),'k'); hold on
    if I3s(i)==2
        plot(i-1+[0:1/T3:1-1/T3],M3(i,:),'r')
    end;
end;
xlim([0 30])
title('warped ECG signal')
xlabel('warped time')

subplot(312)
for i = 1:length(I3s)
    plot(i-1+[0:1/T3:1-1/T3],M3(i,:),'k'); hold on
    if I3s(i)==2
        plot(i-1+[0:1/T3:1-1/T3],M3(i,:),'r')
    end;
end;
xlim([30 60])
xlabel('warped time')


