function result=ConvR4(u,u_temp,i,dx)
  result=((4*ConvR2(u,u_temp,i,dx))/3)-...
         ((1/3)*((1/(32*dx))*(u(i+2,1).*u(i+2,1)+(2*u(i+2,1).*u(i+2,1))+...
          (2*u(i+2,1).*u(i,1))-(2*u(i,1).*u(i-2,1))+...
          (2*u(i+2,1).*u_temp(i,1))-...
          (2*u(i,1).*u_temp(i-2,1))+(u_temp(i+2,1).*u_temp(i+2,1))+(2*...
          u_temp(i+2,1).*u(i,1))-...
          (2*u_temp(i,1).*u(i-2,1))+(2*u_temp(i+2,1)...
          .*u_temp(i,1))-(2*u_temp(i,1).*u_temp(i-2,1))-...
          (u(i-2,1).*u(i-2,1))-(2*u(i-2,1).*u_temp(i-2,1))-...
          (u_temp(i-2,1).*u_temp(i-2,1)))));
end