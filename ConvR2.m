function result=ConvR2(u,u_temp,i,dx)
result= (1/(16*dx))*(u(i+1,1).*u(i+1,1)+(2*u(i+1,1)).*u(i+1,1)+...
          (2*u(i+1,1)).*u(i,1)-2*u(i,1).*u(i-1,1)+2*u(i+1,1).*u_temp(i,1)-...
          (2*u(i,1).*u(i-1,1))+(u(i+1,1).*u(i+1,1))+(2*u(i+1,1).*u(i,1))-...
          (2*u_temp(i,1).*u(i-1,1))+(2*u(i+1,1).*u_temp(i,1))-...
          (2*u_temp(i,1).*u(i-1,1))-...
          u(i-1,1).*u(i-1,1)-2*u(i-1,1).*u(i-1,1)-u(i-1,1).*u(i-1,1));
end