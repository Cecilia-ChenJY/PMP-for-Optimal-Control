x=[0.0025,0.005,0.01,0.02];%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
%epoch from 40000 to 41000
epoch1=[0.0816427171,0.154962778,0.32683205604,0.4471741616,];%epoch10000   
epoch2=[0.0816822425,0.155119508,0.3268330097,0.4471741914,];%epoch20000
epoch3=[0.0817213729,0.155270114541,0.3268337249,0.4471741914,];%epoch30000
epoch4=[0.0817601383,0.1554147303,0.3268342316,0.4471741914,];%epoch40000
epoch5=[0.0817985460,0.155553475,0.326834559,0.4471741914,];%epoch44000
 %dt=0.02 a=[-0.4471741616,-0.4471741914,-0.4471741914,-0.4471741914,-0.4471741914]; %dt=0.02(一共每隔1000步取一次，取5个点）
 b=[]; %dt =0.01
 %dt = 0.005 c=[-0.154962778,-0.155119508,-0.155270114541,-0.1554147303,-0.155553475]; 
 %dt = 0.0025 d=[-0.0816427171,-0.0816822425,-0.0817213729,-0.0817601383,-0.0817985460]; 
 plot(x,epoch1,'-*r','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','r');
 hold on;
 plot(x,epoch2,'-dc','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','c');
 hold on;
 plot(x,epoch3,'--pb','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','b')
 hold on;
 plot(x,epoch4,'-.og','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','g')
 hold on;
 plot(x,epoch5,'--dy','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','y')
axis([0,0.021,0.05,0.45])  %确定x轴与y轴框图大小
set(gca,'XTick',[0.0025:0.0025:0.02]) %x轴范围1-6，间隔1
%set(gca,'YTick',[0:100:700]) %y轴范围0-700，间隔100
%legend('epoch1','epoch2','epoch3','epoch4');   %右上角标注
xlabel('time interval')  %x轴坐标描述
ylabel('errors ') %y轴坐标描述