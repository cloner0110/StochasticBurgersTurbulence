%RUN MAIN.m
function result=AveragingFunc(u,nx)
  sigma=0;
  for i=1:nx
      sigma=sigma+u(i,1);
  end
  result=sigma/nx;
end