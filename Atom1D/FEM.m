function [state,value]=FEM(flag,V,x)
% matrix element H=-1/2*Lapace+V
% differential operator without periodic boundary condition is easy to represent by diag function
% but we write the assignment format just for a better understand
xmin=x(1);xmax=x(end);
nx=length(x);
dx=(xmax-xmin)/(nx-1);
H=zeros(nx,nx);
switch flag
    case '3'
        % three point difference
        for i=2:nx-1
            % The first term is the second order difference, and the second term is the first order difference
            H(i,i+1)=-1/2/dx/dx;
            H(i,i-1)=-1/2/dx/dx;
        end
        % start and end item
        H(1,2)=-1/2/dx/dx;
        H(nx,nx-1)=-1/2/dx/dx;
        %{ 
        periodic boundary condition
        H(1,end)=H(2,1);
        H(end,1)=H(end-1,end);
        %}
        H=H+eye(nx)/dx/dx;
    case '5'
        % five point difference
        for i=2:nx-1
            % The first term is the second order difference, and the second term is the first order difference
            H(i,i+1)=-2/3/dx/dx;
            H(i,i-1)=-2/3/dx/dx;
        end
        % start and end item
        H(1,2)=-2/3/dx/dx;
        H(nx,nx-1)=-2/3/dx/dx;
        for i=3:nx-2
            % The first term is the second order difference, and the second term is the first order difference
            H(i,i+2)=1/24/dx/dx;
            H(i,i-2)=1/24/dx/dx;
        end
        H(1,3)=1/24/dx/dx;
        H(2,4)=1/24/dx/dx;
        H(end-1,end-3)=1/24/dx/dx;
        H(end,end-2)=1/24/dx/dx;
        %{ 
        periodic boundary condition
        H(1,end)=H(2,1);
        H(end,1)=H(end-1,end);
        H(1,end-1)=H(3,1);
        H(2,end)=H(3,1);
        H(end-1,1)=H(end-2,end);
        H(end,2)=H(end-2,end);
        %}
        H=H+eye(nx)*5/dx/dx/4;
    case '7'
        % seven point difference
        for i=2:nx-1
            % The first term is the second order difference, and the second term is the first order difference
            H(i,i+1)=-3/4/dx/dx-1i*k*3/4/dx;
            H(i,i-1)=-3/4/dx/dx+1i*k*3/4/dx;
        end
        % start and end item
        H(1,2)=-3/4/dx/dx-1i*k*3/4/dx;
        H(nx,nx-1)=-3/4/dx/dx+1i*k*3/4/dx;
        for i=3:nx-2
            % The first term is the second order difference, and the second term is the first order difference
            H(i,i+2)=3/40/dx/dx+1i*k*3/20/dx;
            H(i,i-2)=3/40/dx/dx-1i*k*3/20/dx;
        end
        H(1,3)=3/40/dx/dx+1i*k*3/20/dx;
        H(2,4)=3/40/dx/dx+1i*k*3/20/dx;
        H(end-1,end-3)=3/40/dx/dx-1i*k*3/20/dx;
        H(end,end-2)=3/40/dx/dx-1i*k*3/20/dx;
        for i=4:nx-3
            % The first term is the second order difference, and the second term is the first order difference
            H(i,i+3)=-1/180/dx/dx-1i*k/60/dx;
            H(i,i-3)=-1/180/dx/dx+1i*k/60/dx;
        end
        H(1,4)=-1/180/dx/dx-1i*k/60/dx;
        H(2,5)=-1/180/dx/dx-1i*k/60/dx;
        H(3,6)=-1/180/dx/dx-1i*k/60/dx;
        H(end-2,end-5)=-1/180/dx/dx+1i*k/60/dx;
        H(end-1,end-4)=-1/180/dx/dx+1i*k/60/dx;
        H(end,end-3)=-1/180/dx/dx+1i*k/60/dx;
        %{ 
        periodic boundary condition
        H(1,end-2)=H(4,1);
        H(1,end-1)=H(3,1);
        H(1,end)=H(2,1);
        H(2,end-1)=H(4,1);
        H(2,end)=H(3,1);
        H(3,end)=H(4,1);
        H(end-2,1)=H(end-3,end);
        H(end-1,1)=H(end-2,end);
        H(end-1,2)=H(end-3,end);
        H(end,1)=H(end-1,end);
        H(end,2)=H(end-2,end);
        H(end,3)=H(end-3,end);
        %}
        H=H+eye(nx)*49/dx/dx/36;
end
H=H+diag(V);
% solve eigen value and eigen vector
[state,value]=eig(H);
state=state./sqrt(sum(conj(state).*state)*dx);
value=diag(value);
end