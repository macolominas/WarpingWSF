x = importdata('62-a_lhl.wav');
fs = x.fs;
x = x.data;
x = resample(x,1,10); %downsampling
x = x(1:end-1);
x = x(:)';
fs = fs/10;

%-----STFT
N = length(x);
redun = 1/4;
fmax = 0.5;
sigma = 0.005;
bt = 1:N; 
ft = 1/redun:1/redun:round(fmax*N);
gamma = 0;
[Fx,~,~,~,~,~,omega2] = sstn_test_modL(x,gamma,sigma,ft,bt);
%------------

c = exridge(Fx,0,0,16*redun);
c = round(c/3);
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fx(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;

sobremuestreo = 1;
t_m = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m,'spline');
y = interp1(1:N,x./A,m,'spline');
[Fy,~,~,~,~,~,omega2] = sstn_test_modL(y,gamma,sigma,ft,bt);

%-------------


P = ((t_m(end)-t_m(1)+1)/(2*pi));
T = floor(N/P); 
P = floor(P);
M = zeros(P,T); 

for i = 1:P
    M(i,:) = interp1(0.5*t_m/pi,y,i-1+[0:1/(T):1-1/(T)],'spline');
end;

%---------------------

eva = evalclusters(M,'kmeans','CalinskiHarabasz','KList',[1:6]);
K = eva.OptimalK;
% K = 3;
I = kmeans(M,K,'replicates',50);
figure;plot(I)
WSF = zeros(K,T); 
regWSF = zeros(K,T);
orden = 7; 
for i = 1:K
    WSF(i,:) = median(M(I==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T:1-1/T,orden);
end;
%------------------------------------------------------

%-----graficas----------------------------------

t = 0:1/fs:N/fs-1/fs;
f = fs/(N*redun):fs/(N*redun):fs*fmax;
aux = 0.5*diff(t_m)/pi;
f_m = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
figure;
subplot(5,2,1); plot(t,x,'k'); set(gca,'xticklabel',[]); xlim([t(1) t(end)]); text(2.3,0.75,'/a/')
subplot(5,2,[3 5 7 9]); imagesc(t,f,abs(Fx)); colormap(1-gray); set(gca,'ydir','normal')
hold on; plot(t,fs*c/(N*redun),'r--')
%hold on; plot(t,0.5*frec,'r--')
xlabel('time [s]')
ylabel('frequency [Hz]')
subplot(5,2,2); plot(0.5*t_m/pi,y,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m(1)/pi 0.5*t_m(end)/pi]); text(500,1.75,'warped /a/')
subplot(5,2,[4 6 8 10]); imagesc(0.5*t_m/pi,f_m,abs(Fy)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')


figure;subplot(2,10,[1 2 3 4 5 11 12 13 14 15])
imagesc(0:2*pi/T:2*pi-2*pi/T,1:P,M); colormap(1-gray)
title('Waveforms')
xlabel('phase')
ylabel('cylces')
subplot(2,10,[6 16]); plot(I,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
subplot(2,10,[7 8 9 10]); plot(0:2*pi/T:2*pi-2*pi/T,WSF(1,:),'k')
text(0.5,1.2,'waveform no. 1')
set(gca,'yticklabel',[])
xlabel('phase')
subplot(2,10,[17 18 19 20]); plot(0:2*pi/T:2*pi-2*pi/T,WSF(2,:),'k')
text(0.5,1.2,'waveform no. 2')
set(gca,'yticklabel',[])
xlabel('phase')



%-------------------------------
t = 0:1/fs:N/fs-1/fs;
f = fs/(N*redun):fs/(N*redun):fs*fmax;
aux = 0.5*diff(t_m)/pi;
f_m = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
figure;
subplot(10,20,[1:5]); plot(t,x,'k'); set(gca,'xticklabel',[]); xlim([t(1) t(end)]); title('original /a/ signal')%text(2.3,0.75,'/a/')
subplot(10,20,[21:25 41:45 61:65 81:85]); imagesc(t,f,abs(Fx)); colormap(1-gray); set(gca,'ydir','normal')
hold on; plot(t,fs*c/(N*redun),'r--')
%hold on; plot(t,0.5*frec,'r--')
xlabel('time [s]')
ylabel('frequency [Hz]')
subplot(10,20,[6:10]); plot(0.5*t_m/pi,y,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m(1)/pi 0.5*t_m(end)/pi]); title('warped signal')%text(500,1.75,'warped /a/')
subplot(10,20,[26:30 46:50 66:70 86:90]); imagesc(0.5*t_m/pi,f_m,abs(Fy)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,20,[11:15 31:35 51:55 71:75 91:95])
imagesc(0:2*pi/T:2*pi-2*pi/T,1:P,M); colormap(1-gray)
title('Waveforms')
xlabel('phase')
ylabel('cycles')
subplot(10,20,[16 36 56 76 96]); plot(I,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
subplot(10,20,[17:20 37:40]); plot(0:2*pi/T:2*pi-2*pi/T,WSF(1,:),'k')
title('waveform no. 1')
%text(0.5,1.2,'waveform no. 1')
set(gca,'yticklabel',[])
xlabel('phase')
subplot(10,20,[57:60 77:80]); plot(0:2*pi/T:2*pi-2*pi/T,WSF(2,:),'k')
title('waveform no. 2')
%text(0.5,1.2,'waveform no. 2')
set(gca,'yticklabel',[])
xlabel('phase')
%------------------------------

%-----------------------------
t = 0:1/fs:N/fs-1/fs;
f = fs/(N*redun):fs/(N*redun):fs*fmax;
aux = 0.5*diff(t_m)/pi;
f_m = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
figure;
subplot(10,15,[1:5]); plot(t,x,'k'); set(gca,'xticklabel',[]); xlim([t(1) t(end)]); title('original /a/ signal')%text(2.3,0.75,'/a/')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(t,f,abs(Fx)); colormap(1-gray); set(gca,'ydir','normal')
hold on; plot(t,fs*c/(N*redun),'r--')
xlabel('time [s]')
ylabel('frequency [Hz]')
subplot(10,15,[76:80]); plot(0.5*t_m/pi,y,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m(1)/pi 0.5*t_m(end)/pi]); title('warped signal')%text(500,1.75,'warped /a/')
subplot(10,15,[91:95 106:110 121:125 136:140]); imagesc(0.5*t_m/pi,f_m,abs(Fy)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[81:85 96:100 111:115 126:130 141:145])
imagesc(0:2*pi/T:2*pi-2*pi/T,1:P,M); colormap(1-gray)
title('Waveforms')
xlabel('phase')
ylabel('cycles')
subplot(10,15,[86 101 116 131 146]); plot(I,fliplr(1:P),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
subplot(10,15,[87:90 102:105]); plot(0:2*pi/T:2*pi-2*pi/T,WSF(1,:),'k')
title('waveform no. 1')
%text(0.5,1.2,'waveform no. 1')
set(gca,'yticklabel',[])
xlabel('phase')
subplot(10,15,[117:120 132:135]); plot(0:2*pi/T:2*pi-2*pi/T,WSF(2,:),'k')
title('waveform no. 2')
%text(0.5,1.2,'waveform no. 2')
set(gca,'yticklabel',[])
xlabel('phase')
%------------------------------


% subplot(10,15,[6:10 21:25]); plot(t,x,'k'); xlim([0.6 0.9])
% xlabel('time [s]')
% subplot(10,15,[11:15 26:30]); plot(t,x,'k'); xlim([1.6 1.9])
% xlabel('time [s]')
% 
% subplot(10,15,[36:40 51:55]); plot(.5*t_m/pi,y,'k'); xlim([135 160])
% xlabel('warped time')
% subplot(10,15,[41:45 56:60]); plot(.5*t_m/pi,y,'k'); xlim([508 532])
% xlabel('warped time')

subplot(10,15,[6:15 21:30]); plot(t,x,'k'); xlim([0.65 0.85]); ylim([-1 1])
hold on
plot([.661 .661],[.8 .9],'r')
plot([.6658 .6658],[.8 .9],'r')
plot([.661 .6658],[.85 .85],'r')
text(.668,.85,[num2str(.6658-.661) ' s'])

plot([.8334 .8334],[.8 .9],'r')
plot([.8368 .8368],[.8 .9],'r')
plot([.8334 .8368],[.85 .85],'r')
text(.817,.85,[num2str(.8368-.8334) ' s'])
title('zoomed in original signal')
xlabel('time [s]')

subplot(10,15,[36:45 51:60]); plot(.5*t_m/pi,y,'k'); xlim([125 170])
hold on
plot([126.7798 126.7798],[1.6 1.8],'r')
plot([127.7778 127.7778],[1.6 1.8],'r')
plot([126.7798 127.7778],[1.7 1.7],'r')
text(128,1.7,[num2str(127.7778-126.7798)])

plot([165.9384 165.9384],[1.6 1.8],'r')
plot([166.9363 166.9363],[1.6 1.8],'r')
plot([165.9384 166.9363],[1.7 1.7],'r')
text(167.1,1.7,num2str(166.9363 - 165.9384))
title('zoomed in warped signal')
xlabel('warped time')

