import numpy as np
import matplotlib.pyplot as plt
import scipy

Nbridge=100000000
Nsteps=50000
Bigdt = 6.25e-2 
dt=Bigdt/Nsteps
X_final=np.zeros(Nbridge)
t=np.arange(1,Nsteps+1)*1.0/Nsteps*Bigdt
for k in range(0,Nbridge):
    W_t=np.cumsum(np.sqrt(dt)*scipy.randn(Nsteps))

    #subtrack linear piece
    #delta=W_t[Nsteps-1]-Bigdt/2
    #W_t=W_t-delta*t/Bigdt
    #plt.figure(1)
    #plt.plot(t,W_t)
    #  Trapezoid rule
    int_term=dt/2*np.exp(0.5)+dt*np.cumsum(np.exp(0.5*t+W_t)) - dt/2*np.exp(0.5*t[-1]+W_t[-1])
    X_t=np.exp(0.5*t+W_t)/(2+int_term)
    #plt.figure(2)
    #plt.plot(t,X_t)
    X_final[k]=X_t[-1]

ave_end=np.mean(X_final)
var_end=np.var(X_final)

print('ave=',ave_end)
print('var=',var_end)
print('sd=',np.sqrt(var_end))
plt.figure(3)
plt.hist(X_final,20)
plt.title('solutions at t_end')

#plt.figure(2)
#plt.title('solutions')
#plt.figure(1)
#plt.title('random walks')
plt.show() 
