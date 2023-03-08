t=0:0.01:0.99;%
length(t)
ux = normrnd(8,3,[1 100]);
uy = normrnd(0,0.08,[1 100]);
plot3(t, ux, uy,'-.g')
grid on
hold on;
%load the optimal control u*,and plot them in one figure
u_op = load('u_op.mat');
ux_op = u_op.y(:,1);
uy_op = u_op.y(:,2);
plot3(t, ux_op, uy_op)