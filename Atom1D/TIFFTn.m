function [x,V,state,value]=TIFFTn
ConstantAtom;
param;
% Momentum
Pexp=exp(-(2*pi)^2*(xshift/(xmax-xmin)).^2*dt/4);
% Coulomb potential, use "soft core potential"
V=potential('cos1',x,0.37,ax);
Vexp=exp(-V*dt);
% Initial Value
u=-V;
u=u/sqrt(sum(conj(u).*u)*dx);
E=sum(conj(u).*(-2*del2(u,dx)+V.*u))*dx;
state=rand(nx,nE);state(:,1)=u;value=zeros(1,nE);
% The ground function is invoked to solve the groud state with imaginary time evolution method
% Note that the frequency of FFT in MATLAB is 0-->max, and we need -max/2-->0-->max/2 in our calculation
err=zeros(1,nE)+100;
% find the large cycle of the excited state
for jj=1:nE
    u=state(:,jj);
    u=u/sqrt(sum(conj(u).*u)*dx);
    ii=0;
    while err(jj)>err0
        ii=ii+1;
        E0=E;
        du=excitestate('1d',u,state,jj,dx);
        u=u-du;
        u=SplitOperator('1d',u,Pexp.',Vexp);
        u=u/sqrt(sum(conj(u).*u)*dx);
        %E=sum(conj(u).*(H*u))*dx;
        E=sum(conj(u).*(-2*del2(u,dx)+V'.*u))*dx;
        err(jj)=abs(E-E0);
        if mod(ii,1e2)==0
            fprintf('loop step ii = %d, err = %e\n',ii,err(jj));
        end
        ii
    end
    value(jj)=E;
    state(:,jj)=u;
    fprintf('state %d : %f  %e\n',jj,value(jj),err(jj));
end
end
