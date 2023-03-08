%画三维系统的向量场
[x,y,z] = meshgrid(linspace(0,3,10),linspace(0,3,10),linspace(0,3,10));%网格,范围，100个点
%参数设置
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
%主函数方程
u1 = D.*(N0-x)-(b/a).*x.*y;
v1 = alpha.*(b./a).*x.*y-(c.*y./(1+d.*y)).*z-D1.*y;
w = beta.*(c.*y./(1+d.*y)).*z-D2.*z;
figure(1)
q = quiver3(x,y,z,u1,v1,w,2);
xlim([0,1]);
%axis equal
%axis([-0.01,0.01,-0.01,0.02,-0.01,0.01])
%text(0.7123,1.9898,0.3627,'SN1','color','r');
q.ShowArrowHead = 'on';
xlabel('x');ylabel('y');zlabel('z')
title(sprintf('The Vector Field of Synthetic Genetic Oscillator'))
figure(2)
u2 = D.*(N0-x)-b.*y;
v2 = alpha.*b.*y-(c.*y./(1+d.*y)).*z-D1.*y;
w = beta.*(c.*y./(1+d.*y)).*z-D2.*z;
q = quiver3(x,y,z,u2,v2,w,2.5);
xlabel('x');ylabel('y');zlabel('z')
xlim([1,3]);


