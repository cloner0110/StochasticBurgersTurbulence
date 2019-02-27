du(1)=(u(2,1)-u(1,1))/dx;
FilteredU=FilterFunction(u,Nx);
FilteredU2=FilterFunction(u.^2,Nx);
du(Nx)=(u(Nx,1)-u(Nx-1,1))/dx;
for i=2:Nx-1
    du(i)=(u(i+1,1)-u(i-1,1))/(2*dx);
end
filteereddu=FilterFunction(du,Nx);
FilteredS=FilterFunction(abs(du).*du,Nx);
for i=1:Nx
    N(i)=-2*(delta^(4/3))*filteereddu(i)*(2^(4/3)-1);
    M(i)=2*(delta^2)*FilteredS(i)-...
           2*(delta2 ^2)*abs(filteereddu(i))*filteereddu(i);
    L11(i)=FilteredU2(i)-FilteredU(i)*FilteredU(i);
end
Cs=sqrt((AveragingFunc(L.*M,Nx))/(AveragingFunc(M.*M,Nx)));
Cw=(AveragingFunc(L.*N,Nx))/(AveragingFunc(N.*N,Nx));
if (Cs^2)<0
    Cs=0;
end
if (Cw)<0
    Cw=0;
end