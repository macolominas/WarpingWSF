opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
% opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["activity", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "ra_x", "ra_y", "ra_z"];
opts.SelectedVariableNames = ["activity", "ra_x", "ra_y", "ra_z"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "double", "double", "double"];
opts = setvaropts(opts, [2, 3, 4, 5, 6, 7, 8, 9, 10, 11], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 3, 4, 5, 6, 7, 8, 9, 10, 11], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
% aux = readtable("C:\Marcelo\WarpingWSF\IU_walking_driving_climbing\raw_accelerometry_data\id00b70b13.csv", opts);
% aux = readtable("C:\Marcelo\WarpingWSF\IU_walking_driving_climbing\raw_accelerometry_data\id3e3e50c7.csv", opts);
% aux = readtable("C:\Marcelo\WarpingWSF\IU_walking_driving_climbing\raw_accelerometry_data\id7c20ee7a.csv", opts);
aux = readtable("id86237981.csv", opts);

A = table2array(aux);

fs = 100;


ankle = sqrt(A(:,2).^2+A(:,3).^2+A(:,4).^2);
x = ankle(78501:83500);
x = x(:)';
x = x - mean(x);



%-----STFT
N = length(x);
redun = 4;%
fmax = 0.05;
sigma = 0.1;
bt = 1:N; %
ft = 1/redun:1/redun:round(fmax*N); %
gamma = 0; %
[Fx,~,~,~,~,~,omega2] = sstn_test_modL(x,gamma,sigma,ft,bt);
%------------

c1 = exridge(Fx(100:end,:),0,0,5*redun);
c1 = c1+100;
frec = zeros(size(x));
A = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c1(i),i); 
    A(i) = abs(Fx(c1(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;

%--
sobremuestreo = 1;
t_m1 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m1,'spline');
y1 = interp1(1:N,x,m,'spline');
[Fy1,~,~,~,~,~,omega2] = sstn_test_modL(y1,gamma,sigma,ft,bt);
%-------------

%-------
t_m1 = t_m1/2;
P = ((t_m1(end)-t_m1(1)+1)/(2*pi));
T1 = floor(N/P); 
P = floor(P);
M1 = zeros(P,T1); 


for i = 1:P
    M1(i,:) = interp1(0.5*t_m1/pi,y1,i-1+[0:1/(T1):1-1/(T1)],'spline');
end;

%---------------------

%---------
c = exridge(Fy1(100:end,:),0,0,20);
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i);
    A(i) = abs(Fy1(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;
t_m2 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m2,'spline');
y2 = interp1(1:N,y1,m,'spline');
[Fy2,~,~,~,~,~,omega2] = sstn_test_modL(y2,gamma,sigma,ft,bt);
%---------------------------------------

%------------
c = exridge(Fy2(100:end,:),0,0,20);%
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); %
    A(i) = abs(Fy2(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;
t_m3 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m3,'spline');%
y3 = interp1(1:N,y2,m,'spline');
[Fy3,~,~,~,~,~,omega2] = sstn_test_modL(y3,gamma,sigma,ft,bt);
%---------------------------------------

%--------
t_m3 = t_m3/2;
P = ((t_m3(end)-t_m3(1)+1)/(2*pi));
T3 = floor(N/P); 
P = floor(P);
M3 = zeros(P,T3); 

for i = 1:P
    M3(i,:) = interp1(0.5*t_m3/pi,y3,i-1+[0:1/(T3):1-1/(T3)],'spline');
end;

%---------------------

%----------
c = exridge(Fy3(100:end,:),0,0,20);
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fy3(c(i),i));
end;
A = 2*pi*A;
phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;% 
t_m4 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m4,'spline');
y4 = interp1(1:N,y3,m,'spline');
[Fy4,~,~,~,~,~,omega2] = sstn_test_modL(y4,gamma,sigma,ft,bt);
%---------------------------------------
%--------------
c = exridge(Fy4(100:end,:),0,0,20);
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fy4(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;%
t_m5 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m5,'spline');
y5 = interp1(1:N,y4,m,'spline');
[Fy5,~,~,~,~,~,omega2] = sstn_test_modL(y5,gamma,sigma,ft,bt);
%---------------------------------------
% 
%-------------
c = exridge(Fy5(100:end,:),0,0,20);
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i);
    A(i) = abs(Fy5(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;%
t_m5 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m5,'spline');%
y5 = interp1(1:N,y5,m,'spline');
[Fy5,~,~,~,~,~,omega2] = sstn_test_modL(y5,gamma,sigma,ft,bt);
%---------------------------------------

%----------------
c = exridge(Fy5(100:end,:),0,0,20);
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fy5(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;
t_m5 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m5,'spline');%
y5 = interp1(1:N,y5,m,'spline');
[Fy5,~,~,~,~,~,omega2] = sstn_test_modL(y5,gamma,sigma,ft,bt);
%---------------------------------------
% 
%-------------
c = exridge(Fy5(100:end,:),0,0,20);
c = c+100;
frec = zeros(size(x));
for i = 1:N
    frec(i) = omega2(c(i),i); 
    A(i) = abs(Fy5(c(i),i)); 
end;
A = 2*pi*A;

phi_est = 2*pi*cumsum(frec)/N;
sobremuestreo = 1;
t_m5 = linspace(min(phi_est),max(phi_est),sobremuestreo*N);
m = interp1(phi_est,1:N,t_m5,'spline');
y5 = interp1(1:N,y5,m,'spline');
[Fy5,~,~,~,~,~,omega2] = sstn_test_modL(y5,gamma,sigma,ft,bt);
% %---------------------------------------


%------
t_m5 = t_m5/2;
P = ((t_m5(end)-t_m5(1)+1)/(2*pi));
T5 = floor(N/P); 
P = floor(P);
M5 = zeros(P,T5);

for i = 1:P
    M5(i,:) = interp1(0.5*t_m5/pi,y5,i-1+[0:1/(T5):1-1/(T5)],'spline');
end;

%---------------------


%figures------------------------
t = 0:1/fs:N/fs-1/fs;
f = fs/(N*redun):fs/(N*redun):fs*fmax;
aux = 0.5*diff(t_m1)/pi;
f_m1 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
figure;
subplot(10,15,1:5); plot(t,x,'k'); set(gca,'xticklabel',[]); xlim([t(1) t(end)]); title('original ABP signal');% text(1,-380,'ECG signal')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(t,f,abs(Fx)); colormap(1-gray); set(gca,'ydir','normal')
hold on; plot(t,fs*c1/(N*redun),'r--');
xlabel('time [s]')
ylabel('frequency [Hz]')
subplot(10,15,76:80); plot(0.5*t_m1/pi,y1,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m1(1)/pi 0.5*t_m1(end)/pi]);
title('after 1 demodulation w/ 2nd harm.')% text(2,-380,'warped ECG signal')
subplot(10,15,[91:95 106:110 121:125 136:140]); imagesc(0.5*t_m1/pi,f_m1,abs(Fy1)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[6:15]); plot(t,x,'k'); xlim([5 10]); title('zoomed in original ABP signal'); xlabel('time [s]')
subplot(10,15,[21:30]); plot(t,x,'k'); xlim([25 30]); title('zoomed in original ABP signal'); xlabel('time [s]')
subplot(10,15,[36:45]); plot(t,x,'k'); xlim([40 45]); title('zoomed in original ABP signal'); xlabel('time [s]')


subplot(10,15,[81:85 96:100 111:115 126:130 141:145]);
imagesc(M1); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M1))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%------
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicates',50));
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

subplot(10,15,[86 101 116 131 146]); plot(I1,fliplr(1:length(I1)),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
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
aux = 0.5*diff(t_m3)/pi;
f_m3 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,15,1:5); plot(0.5*t_m3/pi,y3,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m3(1)/pi 0.5*t_m3(end)/pi]);
title('after 3 demodulations w/ 2nd harm.');% text(1,-380,'ECG signal')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(0.5*t_m3/pi,f_m3,abs(Fy3)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[6:10 21:25 36:40 51:55 66:70]); imagesc(M3); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M3))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%--------
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicates',50));
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

subplot(10,15,[11 26 41 56 71]); plot(I3,fliplr(1:length(I3)),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K3,5)
    subplot(10,15,[15*(i)-3:15*(i)]); plot(0:2*pi/T3:2*pi-2*pi/T3,WSF(i,:),'k')
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;

aux = 0.5*diff(t_m5)/pi;
f_m5 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,15,[76:80]); plot(0.5*t_m5/pi,y5,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m5(1)/pi 0.5*t_m5(end)/pi]);
title('after 8 demodulations w/ 2nd harm.');% text(1,-380,'ECG signal')
subplot(10,15,[91:95 106:110 121:125 136:140]); imagesc(0.5*t_m5/pi,f_m5,abs(Fy5)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[81:85 96:100 111:115 126:130 141:145]); imagesc(M5); colormap(1-gray)
title(['Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M5))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%-------
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicates',50));
eva = evalclusters(M5,myfunc,'CalinskiHarabasz','KList',[1:6]);
K5 = eva.OptimalK;
I5 = kmeans(M5,K5,'replicates',50);
WSF = zeros(K5,T5); 
regWSF = zeros(K5,T5);
orden = 7; 
for i = 1:K5
    WSF(i,:) = median(M5(I5==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T5:1-1/T5,orden);
end;
%------------------------------------------------------

subplot(10,15,[86 101 116 131 146]); plot(I5,fliplr(1:length(I5)),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K5,5)
    subplot(10,15,[15*(i+5)-3:15*(i+5)]); plot(0:2*pi/T5:2*pi-2*pi/T5,WSF(i,:),'k')
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;


figure;
M5s = synchronization(M5);
aux = 0.5*diff(t_m5)/pi;
f_m5 = 1/(aux(1)*N*redun):1/(aux(1)*N*redun):fmax/aux(1);
subplot(10,15,1:5); plot(0.5*t_m5/pi,y5,'k'); set(gca,'xticklabel',[]); xlim([0.5*t_m5(1)/pi 0.5*t_m5(end)/pi]);
title('after 8 demodulations w/ 2nd harm.');% text(1,-380,'ECG signal')
subplot(10,15,[16:20 31:35 46:50 61:65]); imagesc(0.5*t_m5/pi,f_m5,abs(Fy5)); colormap(1-gray); set(gca,'ydir','normal')
xlabel('warped time')
ylabel('warped frequency')

subplot(10,15,[6:10 21:25 36:40 51:55 66:70]); imagesc(M5s); colormap(1-gray)
title(['Synchronized Waveform Matrix; SVD entropy = ' num2str(svd_entropy(M5s))])
set(gca,'xticklabel',[])
xlabel('phase')
ylabel('cycles')

%------
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',50));
eva = evalclusters(M5s,myfunc,'CalinskiHarabasz','KList',[1:6]);
K5s = eva.OptimalK;
I5s = kmeans(M5s,K5s,'replicates',50);
WSF = zeros(K5s,T5); 
regWSF = zeros(K5s,T5);
orden = 7; 
for i = 1:K5s
    WSF(i,:) = median(M5s(I5s==i,:),1);
    regWSF(i,:) = trigon_reg(WSF(i,:),0:1/T5:1-1/T5,orden);
end;
%------------------------------------------------------

subplot(10,15,[11 26 41 56 71]); plot(I5s,fliplr(1:length(I5s)),'ro-','markersize',3); ylim([1 P]);%camroll(-90); 
%set(gca,'xticklabel',[]); 
set(gca,'yticklabel',[])
title('Clustering')
for i = 1:min(K5s,5)
    subplot(10,15,[15*(i)-3:15*(i)]);
    plot(0:2*pi/T5:2*pi-2*pi/T5,transpose(M5s(I5s==i,:)))
    hold on
    plot(0:2*pi/T5:2*pi-2*pi/T5,WSF(i,:),'k','linewidth',2)
    
    title(['waveform no. ' num2str(i)])
    set(gca,'yticklabel',[])
    xlabel('phase')
    xlim([0 2*pi])
end;


colores = ['b','r','m','k','g','c'];
figure;
subplot(411)
for i = 1:length(M5s(:,1))
    plot(i-1:1/T5:i-1/T5,M5(i,:),colores(I5s(i)));
    hold on
end;
title('warped IMU signal')
xlim([0 50])
ylim([-1.5 7])
xlabel('warped time')
for i = 1:max(I5s)
    plot([8*(i-1)+1 8*(i-1)+2],[6 6],colores(i));
    text(8*(i-1)+2.3,6,['waveform no. ' num2str(i)]);
end;


subplot(412)
for i = 1:120
    plot(i-1:1/T5:i-1/T5,M5(i,:),colores(I5s(i)));
    hold on
end;
xlim([40 80])
ylim([-1.2 2.5])
xlabel('warped time')
subplot(413)
for i = 1:120
    plot(i-1:1/T5:i-1/T5,M5(i,:),colores(I5s(i)));
    hold on
end;
xlim([80 120])
ylim([-1.2 2.5])
xlabel('warped time')



figure;
subplot(411)
plot(t,x); title('original ABP signal');
xlabel('time (s)')

subplot(412);
plot(.5*t_m1/pi,y1); title('after 1 demodulation')
xlabel('warped time')

subplot(413);
plot(.5*t_m1/pi,y3); title('after 3 demodulations')
xlabel('warped time')

subplot(414);
plot(.5*t_m1/pi,y5); title('after 5 demodulations')
xlabel('warped time')

