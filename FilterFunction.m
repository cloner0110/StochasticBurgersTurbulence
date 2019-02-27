%RUN MAIN.m
function result=FilterFunction(u,Nx)
result=zeros(Nx,1);
for i=2:Nx-1
     result(i,1)=(u(i-1)+2*u(i)+u(i+1))/4; 
end
  result(1,1)=u(1);
  result(Nx,1)=u(Nx); 
end