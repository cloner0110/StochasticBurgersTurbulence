function result=RSWONG(u,i,dx,delta,Cs)
  result=-1*2*(Cs*delta^(4/3))*((u(i+1,1)-u(i-1,1))/(2*dx));
end