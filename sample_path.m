% Generate sample paths for Gaussian process
% Generate sample paths for the SDE by MC
% dX_t = f(X_t)dt +d1*dB_t where B_t is Gaussian process
clear, clc;
randn('state',100);
% rand('state',100);
format long;

tic;

T = 1; 
f = @(x)      x-x^3;  
d1 = @(x)       1;%º”–‘‘Î…˘

M = 20000; dt=T/M;
t=0:dt:T;
N=3000;
X = zeros(M+1,1);
X(1) = -1;             %initial point
n=1;
P = [];
while n<N
for i = 1:M
     W = sqrt(dt)*randn;
     X(i+1) = X(i) + f(X(i))*dt + d1(X(i))*W; 
end
a=max(abs(X));
if abs(X(end)-1)<0.1
   P = [P X];
%plot(t,X);hold on
%axis([-1.5 1.5 -1.5 1.5]);
end
n=n+1;
end

toc

%plot(t',P)

x_op = zeros(1,M+1);
for i = 1:M+1
    [F,xi] = ksdensity(P(i,:));
    [m,I] = max(F);
    x_op(i) = xi(I); 
end
figure
plot(t,x_op)

% save  data t X
%t=linspace(0,1,1000);
%plot(t,x,'*','LineWidth',3)