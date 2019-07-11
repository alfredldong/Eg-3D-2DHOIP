% The following basic laser parameters are set to be transformed into atomic units
% See this in the Wikipedia deltailly
clear;clc;
addpath(genpath('D:/Program/lib'));
addpath(genpath('../lib'));
ConstantAtom;
paramn;
C=linspecer(nE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=potential('cos1',x,0.37,ax);
[state,value]=FEM('5',V,x);
energy=V0*real(value);
figure
scatter(1:nE,energy(1:nE))
box on
set(gca,'linewidth',2.5);
set(gca,'FontSize',24);
set(gca,'Fontname', 'Calibri')
set(gcf,'paperpositionmode','auto');
print(gcf,'-dtiff','-r300',strcat('./out/EigState1D.tiff'))
close all
hold on
subplot(211);
plot(x,V,'b','linewidth',2);
set(gca,'linewidth',1);
set(gca,'FontSize',14);
subplot(212);
u=state(:,121);
plot(x,abs(u).^2,'r','linewidth',2);
hold off
box on
set(gca,'linewidth',1);
set(gca,'FontSize',14);
print(gcf,'-dtiff','-r500',strcat('./out/Density1D.tiff'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,V,u,E]=TIFFTn;
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
u=state(:,122);
u=u/sqrt(sum(conj(u).*u)*dx);
[t,dt,Ex,At,At0]=TDFFTn(u');
f=(1:nt/2)/dt/nt;
At=At'.*hann(length(At),'periodic');
Y=fft(At)*dt;
Yabs=abs(Y).^2;
semilogy(f*2*pi/omegax1,Yabs(1:nt/2),'r','linewidth',2)
xlim([0 140])
box on
set(gca,'linewidth',1);
set(gca,'FontSize',14);
print(gcf,'-dtiff','-r500',strcat('./out/HHG1D.tiff'))