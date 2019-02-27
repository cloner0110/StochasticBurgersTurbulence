function result=RSSMAG(u,i,dx,delta,Cs)
  result=-1*2*((Cs*delta)^2)*(abs((u(i+1,1)-u(i-1,1))/(2*dx)))*...
      ((u(i+1,1)-u(i-1,1))/(2*dx));
end