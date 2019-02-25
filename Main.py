import numpy as np
import matplotlib.pyplot as plt
import random
import math
def SBE_FBM(beta,m):
    x=math.sqrt(m)*random.Random(1,m)
    L=m/2;
    freqM=np.zeros(2*l)
    freqM[0]=1
    freqM[(2*l)-1]=1
    for i in range (1,2*l):
        if i<=l:
            freqM[i]=i
        else:
            freqM[i]=(2*l)-i
    xft=np.fft.fft(x)
    xft[0]=1;xft[L]=0
    x1ft=np.multiply(xft,np.power(freqM,beta/2))
    x1=np.real(np.fft.ifft(x1ft))
    return x1
def Averaging(u,nx):
    sigma=0
    for i in range(nx):
        sigma+=u[i]
    return sigma/nx
def FilterFunction(u,nx):
    result=np.zeros(nx)
    for i in range(1,nx):
        result[i]=(u[i-1]+2*u[i]+u[i+1])/4
    result[0]=u[0]
    result[nx-1]=u[nx-1]
    return result
def convR2(u,u_temp,i,dx):
    result=(1/16*dx)
