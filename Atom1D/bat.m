% The following basic laser parameters are set to be transformed into atomic units
% See this in the Wikipedia deltailly
clear;clc;
addpath(genpath('D:/Program/lib'));
addpath(genpath('../lib'));
ConstantAtom;
param;
C=linspecer(nE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,V,u,E]=TIFFT;
hold on
subplot(211);
plot(x,V,'b','linewidth',2);
set(gca,'linewidth',1);
set(gca,'FontSize',14);
subplot(212);
plot(x,abs(u).^2,'r','linewidth',2);
hold off
box on
set(gca,'linewidth',1);
set(gca,'FontSize',14);
print(gcf,'-dtiff','-r500',strcat('./out/Density1D.tiff'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,dt,Ex,At,At0]=TDFFT(u);
hold on
plot(t*omegax1/2/pi,Ex,'b','linewidth',2);
hold off
box on
set(gca,'linewidth',1);
set(gca,'FontSize',14);
print(gcf,'-dtiff','-r500',strcat('./out/Ex.tiff'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HHG
% The x-coordinate of the spectrum, the unit is order
At=At'.*hann(length(At),'periodic');
At0=At0'.*hann(length(At0),'periodic');
f=(1:nt/2)/dt/nt;
Y=fft(At)*dt;
Y0=fft(At0)*dt;
Yabs=abs(Y).^2;
Y0abs=abs(Y0).^2;
figure
subplot(211);
semilogy(f*2*pi/omegax1,Yabs(1:nt/2),'r','linewidth',2);
xlim([0 150])
subplot(212);
semilogy(f*2*pi/omegax1,Y0abs(1:nt/2),'b','linewidth',2);
xlim([0 150])
print(gcf,'-dtiff','-r500',strcat('./out/I_omega.tiff'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier inverse transformation can obtain time domain pulse
% Attosecond
I1=zeros(1,nt);
for ii=round(5*omegax1*TF1/2/pi):round(105*omegax1*TF1/2/pi)
    I1=I1+Y(ii)*exp(i*2*pi*f(ii)*t)*2*pi/TF1;
end
I1=(abs(I1)).^2;
figure
plot(t*omegax1/2/pi,I1,'b','linewidth',2);
print(gcf,'-dtiff','-r500',strcat('./out/I_t.tiff'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time domain analysis
% Spectrogram
figure
specgram(At,f*2*pi/omegax1,2*pi/dt/omegax1,kaiser(200,13),190);
ylim([0 250])
box on
set(gca,'linewidth',1);
set(gca,'FontSize',14);
print(gcf,'-dtiff','-r500',strcat('./out/omega_time1.tiff'))

figure
spectrogram(At, 200, 90, f*2*pi/omegax1, max(f*2*pi/omega)*2, 'yaxis');
ylim([0 250])
box on
set(gca,'linewidth',1);
set(gca,'FontSize',14);
print(gcf,'-dtiff','-r500',strcat('./out/omega_time2.tiff'))
