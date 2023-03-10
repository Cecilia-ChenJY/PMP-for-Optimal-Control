# -*- coding: utf-8 -*-
"""
Jianyu Chen

This is a project for optimal control.
Python code for 1-dim double well system
"""
#%%
import torch
from torch import nn
from torch.nn import functional as F
import math
import numpy as np
import matplotlib.pyplot as plt
#%%
t = np.arange(100) * 0.01
print(t)
sigma = 1
#%%
u= np.random.normal(0, 0.008, 100)
#u = np.linspace(0, 1, 100)
print(u)
sigma = 1 
def runge_kutta(x, t, dt, f, u, i):
    """   x is the initial value for x   
          t is the initial value for t
          dt is the time step in t
          f is derivative of function x(t)

    """

    k1 = dt * (f(x, t)+u[i])
    k2 = dt * (f(x + k1, t + dt)+u[i+1])
    return x + (k1 + k2) / 2.

def get_x(u):
    t = 0.  
    x = -1.  
    T = 1.0 
    dt = 0.01
    xs, ts = [x], [0.0]
    def func(x, t):
        return x-x**3
    i = 0
    while i < 99:
        #print (i, t)
        x = runge_kutta(x, t, dt, func, u, i)
        i += 1
        t += dt
        xs.append(x)
        ts.append(t)
    return np.array([np.float32(xx) for xx in xs])
#%%
import math
import numpy as np
import matplotlib.pyplot as plt 
def runge_kutta(p, t, dt, f, xs, i):
    """   p is the terminal value for p   
          t is the initial value for t
          dt is the time step in t
          f is derivative of function p(t)

    """

    k1 = -dt * (f(p,t)*(3*xs[i]**2-1)-3*xs[i])
    k2 = -dt * (f(p + k1, t + dt)*(3*xs[i-1]**2-1)-3*xs[i-1])                      
    return p + (k1 + k2) / 2.


def get_p(x):
    t = 1.0 
    p = -2*(x[99]-1)/((x[99]-1)**2+1)**2
    dt = 0.01
    ps, ts = [p], [0.0]
    def func(p, t):
        return p
    i = 99
    while i > 0:
        #print (i, t)
        p = runge_kutta(p, t, dt, func, x, i)
        i -= 1
        t -= dt
        ps.append(p)
        ts.append(t)

    return np.array([np.float32(pp) for pp in ps])
#%%
class MLP(nn.Module):
     def __init__(self):#call the constractor class"Module" to perform necessary initialization.
        super().__init__()
        self.hidden1 = nn.Linear(1, 8)#hidden
        self.hidden2 = nn.Linear(8, 16)#hidden
        self.out = nn.Linear(16,1)#out put layer
 #define the forward propagation of the model,how to return the required model output 
 #based on the input "t"（input t is "u"） 
     def forward(self, u):
        u = self.hidden1(u)
        u = F.relu(u)
        u = self.hidden2(u)
        u = F.relu(u)
        out = self.out(u)
        return out
net = MLP()
#%%
t = torch.tensor(t)
t = t.to(torch.float32)
x = get_x(u)
p = get_p(x)
reversed_p = []
for i in reversed(p):
    reversed_p.append(i)
    
print(reversed_p)
p = reversed_p
x = torch.tensor(x).to(torch.float32)
p = torch.tensor(p).to(torch.float32)
#%%
print(x)
print(p)
#%%
import numpy as np
import matplotlib.pylab as plt
#from gradient import numerical_gradient
Loss = []
lr =0.05
epoch = 1000
N = 20#iteration times
Q = 2*(x[99]-1)**2
print(Q)
dt = 0.01
dx= (x[1:]-x[0:-1])/dt
dp = (p[1:]-p[0:-1])/dt
rho = 1 #修正系数
updater = torch.optim.SGD(net.parameters(), lr)  # change to adaptive
print(dx.shape)
#%%
def train():

    for ii in range(epoch):
        H = 0
        for j in range(99): 
            t_j = t[j].reshape(1) 
            u_hat_j = net(torch.tensor(t_j))
            x_j = x[j].reshape(1)
            p_j = p[j].reshape(1)
            H = H + p_j* (x_j-x_j**3+sigma*u_hat_j) - 0.5*(u_hat_j**2+1-3*x_j**2)-0.5*rho*(abs(dx[j]-x_j+x_j**3-sigma*u_hat_j)**2)+abs((dp[j]+p_j*(1-3*x_j**2)+3*x_j)**2)
            

    
        l = -(H-Q)/100 # l>0

        updater.zero_grad()
        l.backward()
        updater.step()
        Loss.append(l.detach().item())
    
    
    return
#%%
for k in range(N):
    print(k)
    U = [net(torch.tensor(t_j.reshape(1))) for t_j in t]
    u_op = [u.detach().numpy() for u in U]
    #print(u_op)
    flatten_u_op = [element for sublist in u_op for element in sublist]
    print(flatten_u_op)
    x = get_x(flatten_u_op)
    p = get_p(x)
    
    reversed_p = []
    for i in reversed(p):
        reversed_p.append(i)
    

    p = reversed_p

    x = torch.tensor(x).to(torch.float32)
    p = torch.tensor(p).to(torch.float32)
    train()
    print(x)
#%%
x = t
y = flatten_u_op
#y = list(map(lambda x:x.detach().numpy(), y))
plt.title("optimal control")
plt.xlabel("Time")
plt.ylabel("u")


plt.plot(x, y)
plt.show()
#%%
x = get_x(y)
print(x)
#%%
import scipy.io as scio
scio.savemat('pathway_data.mat',{'x':x})




