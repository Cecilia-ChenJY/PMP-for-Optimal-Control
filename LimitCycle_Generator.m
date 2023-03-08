%% compute the limit cycle in three dimension

clear,clc;
tic

%% initialization parameters
D = 0.1;
D1 = 0.2;
D2 = 2.1;
a = 1;
b = 1;
c = 5;
d = 0.1;
alpha = 1;
beta = 0.5;
N0 = 9.96;
%g(y) = c*y/(1+d*y);
%f1(x) = (b/a)*x;(0<N<a)
%f2(x) = b;（N>a)

%% N-P-Z system: vector filed (f,g,h)
%N 营养物质； P 浮游植物； Z 浮游生物
F1 = @(x,y,z)    D*(N0-x)-(b/a)*x*y;
F2 = @(x,y,z)    D*(N0-x)-b*y;
G1 = @(x,y,z)    alpha*(b/a)*x*y-(c*y/(1+d*y))*z-D1*y;
G2 = @(x,y,z)    alpha*b*y-(c*y/(1+d*y))*z-D1*y;
H = @(x,y,z)    beta*(c*y/(1+d*y))*z-D2*z;
%% Compute the trajectory using Euler scheme
T = 500;
N = 5*1e5;
dt = T/N;
t = 0 : dt: T;
x(1) =1;
y(1) =1;
z(1) =1;

for i = 1 : N
    if x(i)> a    %use f2 and g2
        x(i+1) = x(i) + F2(x(i),y(i),z(i))*dt;
        y(i+1) = y(i) + G2(x(i),y(i),z(i))*dt;
        z(i+1) = z(i) + H(x(i),y(i),z(i))*dt;
    else
        x(i+1) = x(i) + F1(x(i),y(i),z(i))*dt;
        y(i+1) = y(i) + G1(x(i),y(i),z(i))*dt;
        z(i+1) = z(i) + H(x(i),y(i),z(i))*dt;
    end
        
end
  plot3(x,y,z,'b')
  hold on
  stable=[x(1:5000);y(1:5000);z(1:5000)];
p = 1;
LCV(1,1) = x(end);
LCV(2,1) = y(end);
LCV(3,1) = z(end);
while p > 0 || i==0
    i = i - 1;
    LCV(1,N-i) = x(i+1);
    LCV(2,N-i) = y(i+1);
    LCV(3,N-i) = z(i+1);
    if norm(LCV(:,1)-LCV(:,N-i))<0.1 && t(end)-t(i)>5
       p = 0;
    end
end
LCV(:,N-i+1)=[x(end);y(end);z(end)];
plot3(LCV(1,:),LCV(2,:),LCV(3,:),'r','LineWidth',3)
xlabel('nutrient N');ylabel('phytoplankton P');zlabel('zooplankton Z');
grid on;
plot3(1.01868,0.91703,0.170157,'g.','MarkerSize',20)
save('LimitCycle_N0','LCV','-ascii','-double');
hold on;
toc