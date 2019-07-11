% space grid
nx=3e3;
xmax=244;
xmin=-xmax;
ax=8;
x=linspace(xmin,xmax,nx);
dx=(xmax-xmin)/(nx-1);
% nx is even
if mod(nx,2)==0
    xshift=[0:nx/2-1,-(nx/2):-1];
% nx is odd
elseif mod(nx,2)==1
    xshift=[0:(nx-1)/2,-(nx-1)/2:-1];
end
%%%%%%%%%%%%%%%%%%%%%
nE=300;
% which bands are used as the initial population bands
nE0= 1:2;
% the samples number in k space is always an odd
nkx=181;
nGx=11;
if nkx==1
    kx = 0;
else
    kx=linspace(-pi/ax,pi/ax,nkx);
end
dkx=2*pi/ax/(nkx-1);
Gx=linspace(-pi/ax*(nGx-1),pi/ax*(nGx-1),nGx);
dGx=2*pi/ax;
%%%%%%%%%%%%%%%%%%%%%
% the error in searching ground state
% energy convergence criterion
err0=1e-15;
%%%%%%%%%%%%%%%%%%%%%
% Iaser x 1 intensity
Ix1= 8.09e11;
% Iaser x 2 intensity
Ix2=0;
% wavelength
lambdax1=3200;
lambdax2=500;
% initial phase, arc
phix1 = 0;
phix2 = 0;
% The intense of the laser field
Ex01=2740*sqrt(Ix1)/E0;
Ex02=2740*sqrt(Ix2)/E0;
% The circle frequency of the laser field (nm->au)
omegax1=2*pi*c0/(lambdax1*nm)/f0;
omegax2=2*pi*c0/(lambdax2*nm)/f0;
% One photoperiod
T1=2*pi/omegax1;
T2=2*pi/omegax2;
% The number of the optical cylce
nTF1=10;nTF2=10;
nTE1=2;nTE2=2;
% Pulse width, full width of nTF1/nTE1 photoperiod
TF1=nTF1*T1;TE1=nTE1*T1;
tmax=TF1;
% Time grid
tmin=0;
%{
nt=1e4;
t=linspace(tmin,tmax,nt);
dt=(tmax-tmin)/(nt-1);
%}
dt=0.03;
nt=fix((tmax-tmin)/dt);
if mod(nt,2)~=0
    nt=nt+1;
end
t=linspace(tmin,tmax,nt);
% Absorber
Axmax=40;