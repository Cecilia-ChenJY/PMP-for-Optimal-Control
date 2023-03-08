% Generate sample paths for asymmetric alpha-stable Levy process
% Generate sample paths for the SDE by MC
% dX_t = f(X_t)dt +d1*dB_t+ h1(X_t)*dL_t where L_t is asymmetric alpha-stable levy process
% dY_t = g(Y_t)dt +d2*dB_t+ h2(Y_t)*dL_t 
% L(t_{n+1})-L(t_n) = (t_{n+1}-t_n)^\frac{1}{alpha}Z_n  where Z_n come

clear, clc;
randn('state',100);
% rand('state',100);
format long;

tic;
%alpha1 = 0.5;
%alpha2 = 0.7;
%beta1 = 0.5;
%beta2 = 0;
gamma = 1;
eps = 0.2;

T = 1; 
f = @(x,y)      x-x^3-gamma*x*y;  
g = @(x,y)     -(1+x^2)*y;
d1 = @(x,y)       eps*(x^2);
d2 = @(x,y)       1;


M = 20000; dt=T/M;
t=0:dt:T;
N=5000;
X = zeros(M+1,1);
Y = zeros(M+1,1);
X(1) = -1;             %initial point
Y(1) = 0;
n=1;
while n<N
for i = 1:M
     W = sqrt(dt)*randn;
     %L1 = stblrnd(alpha1,beta1,gamma,delta);
     %L2 = stblrnd(alpha2,beta2,gamma,delta);
     X(i+1) = X(i) + f(X(i),Y(i))*dt + d1(X(i),Y(i))*W; 
     Y(i+1) = Y(i) + g(X(i),Y(i))*dt + d2(X(i),Y(i))*W;
end
a=max(abs(X));
b=max(abs(Y));
if abs(X(end)-1)<1.5 && abs(Y(end)-0)<0.3 && abs(a)<3.5 && abs(b)<5
plot(t,X);
%plot(t,Y);
%plot3(t,X,Y);hold on
%axis([-1.5 1.5 -1.5 1.5]);
end
n=n+1;
end

toc
% save  data t X