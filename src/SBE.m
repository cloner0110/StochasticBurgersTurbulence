clear all; clc; close all;
visc = 1e-5; D = 1e-6; beta = -0.75;
Nx  = 8192; dx = 2*pi/(Nx-1); NSteps = 2e6; dt = 1e-4; NDNS = 8192;
delta=dx;m=Nx;delta2=delta*2;Cs=0.15;Cw=0.009;
e=0.000001;u=zeros(Nx,1);u_temp=zeros(Nx,1);M=zeros(Nx,1);L=zeros(Nx,1);
N=zeros(Nx,1);x=0:dx:2*pi;Ninfo = 100; Nstat = 1000;
%Initial velocity field******************************************************
randn('state',0);                   %Initialize Random viscmber Generator
%%chosing the solution Type
Solutionmodel=questdlg('Which method you want to model with ?',...
    'Question',...
    'DNS','Smagorinsky','Wong','DNS');
modelorder=questdlg('which order of accuracy you want to solve with ?'...
    ,'Question',...
    '2nd order RK','4th order RK','2nd order RK');
%%% Time Marching
for s = 1:200%NSteps
message=msgbox('Your computer is now solving the problem, Please wait..... ');
switch Solutionmodel
case 'DNS'
        if modelorder=='2nd order RK'
            fprintf('%d\n',s);
            f= sqrt(2*D/dt)*SBE_FBM(beta,m)'; 
            ff=FilterFunction(f,Nx);
            for i=3:Nx-2
                  RK1=DiffR2(u,u_temp,i,visc,dx)+ff(i)...
                      -ConvR2(u,u_temp,i,dx);
                  u_temp(i,1)=u_temp(i,1)+0.5*RK1*dt;
                  RK2=DiffR2(u,u_temp,i,visc,dx)+ff(i)...
                      -ConvR2(u,u_temp,i,dx);
                  u(i,1)=u_temp(i,1)+(0.5*RK1+0.5*RK2)*dt;
            end
            end
            
        if modelorder=='4th order RK'
            fprintf('%d\n',s);
            f= sqrt(2*D/dt)*SBE_FBM(beta,m)'; 
            ff=FilterFunction(f,Nx);
            for i=3:Nx-2
                  RK1=DiffR4(u,u_temp,i,visc,dx)+ff(i)...
                      -ConvR4(u,u_temp,i,dx);
                  u_temp(i,1)=u_temp(i,1)+0.5*RK1*dt;
                  RK2=DiffR4(u,u_temp,i,visc,dx)+ff(i)...
                      -ConvR4(u,u_temp,i,dx);
                   u_temp(i,1)=u_temp(i,1)+0.5*RK2*dt;
                  RK3=DiffR4(u,u_temp,i,visc,dx)+ff(i)...
                      -ConvR4(u,u_temp,i,dx);
                   u_temp(i,1)=u_temp(i,1)+RK3*dt;
                  RK4=DiffR4(u,u_temp,i,visc,dx)+ff(i)...
                      -ConvR4(u,u_temp,i,dx);
                  u(i,1)=u_temp(i,1)+(1/6)*(RK1+2*RK2+2*RK3+RK4)*dt;
            end
        end
case 'Smagorinsky'
    if modelorder=='2nd order RK'
        fprintf('%d\n',s);
        f= sqrt(2*D/dt)*SBE_FBM(beta,m)'; 
        ff=FilterFunction(f,Nx);
        for i=3:Nx-2
                  RK1=DiffR2(u,u_temp,i,visc,dx)+ff(i)-0.5*...
                  (RSSMAG(u,i+1,dx,delta,Cs)-...
                  RSSMAG(u,i-1,dx,delta,Cs))/(2*dx)...
                  -ConvR2(u,u_temp,i,dx);
                  u_temp(i,1)=u_temp(i,1)+0.5*RK1*dt;
                  RK2=DiffR2(u,u_temp,i,visc,dx)+ff(i)-...
                  0.5*(RSSMAG(u,i+1,dx,delta,Cs)-...
                  RSSMAG(u,i-1,dx,delta,Cs))/(2*dx)...
                  -ConvR2(u,u_temp,i,dx);
            u(i,1)=u_temp(i,1)+(0.5*RK1+0.5*RK2)*dt;
        end
    end
    if modelorder=='4th order RK'
        fprintf('%d\n',s);
        f= sqrt(2*D/dt)*SBE_FBM(beta,m)'; 
        ff=FilterFunction(f,Nx);
        for i=4:Nx-3
            RK1=DiffR4(u,u_temp,i,visc,dx)+ff(i)-0.5*...
            (RSSMAG(u,i+1,dx,delta,Cs)-...
            RSSMAG(u,i-1,dx,delta,Cs))/(2*dx)...
             -ConvR4(u,u_temp,i,dx);
            u_temp(i,1)=u_temp(i,1)+0.5*RK1*dt;
            RK2=DiffR4(u,u_temp,i,visc,dx)+ff(i)-0.5*...
            (RSSMAG(u,i+1,dx,delta,Cs)-...
            RSSMAG(u,i-1,dx,delta,Cs))/(2*dx)...
            -ConvR4(u,u_temp,i,dx);
            u_temp(i,1)=u_temp(i,1)+0.5*RK2*dt;
            RK3=DiffR4(u,u_temp,i,visc,dx)+ff(i)-0.5*...
            (RSSMAG(u,i+1,dx,delta,Cs)-...
            RSSMAG(u,i-1,dx,delta,Cs))/(2*dx)...
            -ConvR4(u,u_temp,i,dx);
            u_temp(i,1)=u_temp(i,1)+RK3*dt;
            RK4=DiffR4(u,u_temp,i,visc,dx)+ff(i)-0.5*...
            (RSSMAG(u,i+1,dx,delta,Cs)-...
            RSSMAG(u,i-1,dx,delta,Cs))/(2*dx)...
            -ConvR4(u,u_temp,i,dx);
            u(i,1)=u_temp(i,1)+(1/6)*(RK1+2*RK2+2*RK3+RK4)*dt;
        end
        end
case 'Wong'
    if modelorder=='2nd order RK'
            fprintf('%d\n',s); 
            f= sqrt(2*D/dt)*SBE_FBM(beta,m)'; 
            ff=FilterFunction(f,Nx);
            for i=3:Nx-2
                  RK1=DiffR2(u,u_temp,i,visc,dx)+ff(i)-0.5*...
                      (RSWONG(u,i+1,dx,delta,Cs)-...
                      RSWONG(u,i-1,dx,delta,Cs))/(2*dx)...
                      -ConvR2(u,u_temp,i,dx);
                  u_temp(i,1)=u_temp(i,1)+0.5*RK1*dt;
                  RK2=DiffR2(u,u_temp,i,visc,dx)+ff(i)-0.5*...
                      (RSWONG(u,i+1,dx,delta,Cs)-...
                      RSWONG(u,i-1,dx,delta,Cs))/(2*dx)...
                      -ConvR2(u,u_temp,i,dx);
                  u(i,1)=u_temp(i,1)+(0.5*RK1+0.5*RK2)*dt;
            end
    end
    if modelorder=='4th order RK'
        fprintf('%d\n',s); 
            f= sqrt(2*D/dt)*SBE_FBM(beta,m)'; 
            ff=FilterFunction(f,Nx);
            for i=3:Nx-2
              RK1=DiffR4(u,u_temp,i,visc,dx)+ff(i)-...
                  0.5*(RSWONG(u,i+1,dx,delta,Cs)-...
                  RSWONG(u,i-1,dx,delta,Cs))/(2*dx)...
                  -ConvR4(u,u_temp,i,dx);
              u_temp(i,1)=u_temp(i,1)+0.5*RK1*dt;
              RK2=DiffR4(u,u_temp,i,visc,dx)+ff(i)-...
                  0.5*(RSWONG(u,i+1,dx,delta,Cs)-...
                  RSWONG(u,i-1,dx,delta,Cs))/(2*dx)...
                  -ConvR4(u,u_temp,i,dx);
               u_temp(i,1)=u_temp(i,1)+0.5*RK2*dt;
              RK3=DiffR4(u,u_temp,i,visc,dx)+ff(i)-...
                  0.5*(RSWONG(u,i+1,dx,delta,Cs)-...
                  RSWONG(u,i-1,dx,delta,Cs))/(2*dx)...
                  -ConvR4(u,u_temp,i,dx);
               u_temp(i,1)=u_temp(i,1)+RK3*dt;
              RK4=DiffR4(u,u_temp,i,visc,dx)+ff(i)-...
                  0.5*(RSWONG(u,i+1,dx,delta,Cs)-...
                  RSWONG(u,i-1,dx,delta,Cs))/(2*dx)...
                  -ConvR4(u,u_temp,i,dx);
              u(i,1)=u_temp(i,1)+(1/6)*(RK1+2*RK2+2*RK3+RK4)*dt;
            end
    end
    if mod (s,1000)==0
    resolved_KE=0.5*AveragingFunc(u.*u,Nx);
    if t>2000 && t<200000
        fmodes=0:1:nx/2;
        fu=fft(u_temp(:,1))/nx;
        fustar=conj(fu(1:nx/2+1));
        E(:,p)=0.5*(fustar.*fu);
      end
      if t==20000
         for i=1:nx/2+1
            Emean(i)=mean(E(i,:)); 
         end
      end
    end    
end
CSCW
u_temp(:,1)=u(:,1);
end
fclose all;
