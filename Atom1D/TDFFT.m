function [t,dt,Ex,At,At0]=TDFFT(u)
ConstantAtom;
param;
% Momentum
Pexp=exp(-1i*(2*pi)^2*(xshift/(xmax-xmin)).^2*dt/4);
% Coulomb potential, use "soft core potential"
V=potential('softcore1',x,1,0.4713);
Vxg=gradient(V,dx);
Ex=Ex01*laser('trapezoid',t,TF1,TE1).*cos(omegax1*t+phix1);
% Initial dipole acceleration
At=zeros(1,nt);
At0=zeros(1,nt);
% Absorber
absorber=sin((xmax-abs(x))/Axmax/2*pi).^(1/8).*(abs(x)>=xmax-Axmax)+(abs(x)<xmax-Axmax);
% The exited function is invoked to solve the exited state with time dependent evolution method 
for ii=1:1:nt
    Vexp=exp(-1i*(V+x.*Ex(ii))*dt);
    u=SplitOperator('1d',u,Pexp,Vexp);
    u=u.*absorber;
    u=u/sqrt(sum(conj(u).*u)*dx);
    At(ii)=-sum(conj(u).*(Vxg+Ex(ii)).*u)*dx;
    At0(ii)=-sum(conj(u).*(Vxg).*u)*dx;
    if mod(ii,1e3)==0
        fprintf('loop step ii = %d\n',ii);   
    end
end
end