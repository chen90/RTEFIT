function exl50flux(shotnum,t1,t2)
step=0.01;
datatime=[num2str(t1),':',num2str(t2),':',num2str(step)];
[flux,~]=exl50dbN(shotnum,'flux\d*',datatime);
for i=1:39
    m(i)=mean(flux(1:end,i));
end
%---------设置积分时间常数，转换为物理单位 Wb
m(1:14)=m(1:14)/500;
m(15:21)=m(15:21)/200;
m(22:32)=m(22:32)/100;
m(33:39)=m(33:39)/200;
if shotnum<1210
    m(33:39)=rot90(m(33:39)); %为了便于画图，更换顺序
    m(1)=m(1)/5;
    m(32)=m(32)/5;
    x=m(15);
    m(15)=m(16);
    m(16)=x;    
end    
m=abs(m);
num=1:39;
%-------------有问题数据点单独画图--------------
% m_err=[0,0,0];
% num_err=num([21,27,37]);
% m([21,27,37])=[];
% num([21,27,37])=[];
%-----------------画图------------------------------
figure;plot(num,m*1000,':o','Color','k','LineWidth',1.5,'MarkerSize',6,'MarkerFace','b','MarkerEdgeColor','b');
% hold on;plot(num_err,m_err*1000,'o','MarkerSize',8,'MarkerFace','k','MarkerEdgeColor','k');
fillall(1,14,'[1,0,0]');
fillall(15,21,'[0,1,0]');
fillall(22,32,'[0,0,1]');
xlabel('$Loop-No$','interpreter','Latex');
ylabel('$\rm\phi (mWb)$','interpreter','Latex');
title('Flux Loop Measurement')
set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');



% pf12=[ 0.0006    0.0004    0.0004    0.0003    0.0002    0.0002    0.0002    0.0002    0.0002    0.0002    0.0003    0.0004    0.0004    0.0006    0.0056    0.0103    0.0147...
%     0.0197    0.0233    0.0250    0.0246    0.0219    0.0188    0.0160    0.0139    0.0125    0.0115    0.0126    0.0141    0.0162    0.0190    0.0218    0.0245    0.0248    0.0232... ]
%     0.0196    0.0146    0.0103    0.0056];
% 
% pf34=[0.0007    0.0007    0.0007    0.0006    0.0006    0.0006    0.0006    0.0006    0.0006    0.0006    0.0006    0.0007    0.0007    0.0007    0.0045    0.0085    0.0127...
%     0.0190    0.0269    0.0365    0.0574    0.0891    0.0862    0.0776    0.0692    0.0631    0.0589    0.0636    0.0699    0.0783    0.0867    0.0890    0.0574    0.0365 ...
%      0.0269   0.0190    0.0127    0.0085    0.0045];
%  
% pf56=[0.0007    0.0008    0.0008    0.0009    0.0009    0.0010    0.0010    0.0010    0.0010    0.0009    0.0009    0.0008    0.0008    0.0007    0.0042    0.0080    0.0117 ...
%     0.0172    0.0237    0.0310    0.0448    0.0632    0.0759    0.0918    0.1089    0.1215    0.1252    0.1206    0.1075    0.1075    0.0746    0.0633    0.0449    0.0311   ...
%     0.0237    0.0173    0.0118    0.0080    0.0042 


end