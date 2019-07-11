function [x,V,u,E]=TIFFT
ConstantAtom;
param;
% Momentum
Pexp=exp(-(2*pi)^2*(xshift/(xmax-xmin)).^2*dt/4);
% Coulomb potential, use "soft core potential"
V=potential('softcore1',x,1,0.4713);
Vexp=exp(-V*dt);
% Initial Value
u=-V;
u=u/sqrt(sum(conj(u).*u)*dx);
E=sum(conj(u).*(-2*del2(u,dx)+V.*u))*dx;
% The ground function is invoked to solve the groud state with imaginary time evolution method
% Note that the frequency of FFT in MATLAB is 0-->max, and we need -max/2-->0-->max/2 in our calculation
err=1;
ii=0;
while err>err0
    ii=ii+1;
    E0=E;
    u=SplitOperator('1d',u,Pexp,Vexp);
    u=u/sqrt(sum(conj(u).*u)*dx);
    E=sum(conj(u).*(-2*del2(u,dx)+V.*u))*dx;
    err=abs(E-E0);
    if mod(ii,1e3)==0
        fprintf('loop step ii = %d, err = %e\n',ii,err);
    end
end
fprintf('loop end ii = %d, err = %e\n',ii,err);
end
