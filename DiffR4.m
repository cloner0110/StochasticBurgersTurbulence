function result=DiffR4(u,u_temp,i,visc,dx)
  result=((4/3)*DiffR2(u,u_temp,i,visc,dx))-...
         ((1/3)*((visc/(8*(dx^2)))*(u(i+2,1)+u(i+2,1)-(2*u(i,1))-...
         (2*u_temp(i,1))+u(i-2,1)+u(i-2,1))));
end