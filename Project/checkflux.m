function checkflux(shotnum,t1,t2)
step=0.01;
datatime=[num2str(t1),':',num2str(t2),':',num2str(step)];

[flux,t]=exl50dbN(shotnum,'flux\d*',datatime);
%%
num=1;RC_inner=100/1000;  %积分时间/1000单位转换为mWb
figure;stackplot({{t,smooth(flux(:,1)/RC_inner,num)/1,'$\phi(mWb)$'}}, [], [], [], ['Inner Flux'], 'Time (s)');%figure 1
hold on;plot(t,smooth(flux(:,2)/RC_inner,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,3)/RC_inner,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,4)/RC_inner,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,5)/RC_inner,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,6)/RC_inner,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,7)/RC_inner,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,8)/RC_inner,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,9)/RC_inner,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,10)/RC_inner,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,11)/RC_inner,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,12)/RC_inner,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,13)/RC_inner,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,14)/RC_inner,num),'linewidth',2.5);
l2=legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14');
set(l2,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);
xlim([-2,10]);

RC_ud=200/1000;
figure;stackplot({{t,smooth(flux(:,15)/RC_ud,num)/1,'$\phi(mWb)$'}}, [], [], [], ['Up & Down Flux'], 'Time (s)');%figure 2
hold on;plot(t,smooth(flux(:,16)/RC_ud,num)/1,'linewidth',2.5);
hold on;plot(t,-smooth(flux(:,17)/RC_ud,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,18)/RC_ud,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,19)/RC_ud,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,20)/RC_ud,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,21)/RC_ud,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,33)/RC_ud,num),'--','linewidth',2.5);
hold on;plot(t,-smooth(flux(:,34)/RC_ud,num),'--','linewidth',2.5);
hold on;plot(t,smooth(flux(:,35)/RC_ud,num),'--','linewidth',2.5);
hold on;plot(t,smooth(flux(:,36)/RC_ud,num),'--','linewidth',2.5);
hold on;plot(t,smooth(flux(:,37)/RC_ud,num),'--','linewidth',2.5);
hold on;plot(t,smooth(flux(:,38)/RC_ud,num),'--','linewidth',2.5);
hold on;plot(t,smooth(flux(:,39)/RC_ud,num),'--','linewidth',2.5);
l3=legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21','Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39');
set(l3,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

RC_outside=100/1000;
figure;stackplot({{t,smooth(flux(:,22)/RC_outside,num)/1,'$\phi(mWb)$'}}, [], [], [], ['OutSide Flux'], 'Time (s)');%figure 3
hold on;plot(t,smooth(flux(:,23)/RC_outside,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,24)/RC_outside,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,25)/RC_outside,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,26)/RC_outside,num)/1,'linewidth',2.5);
hold on;plot(t,-smooth(flux(:,28)/RC_outside,num)/1,'linewidth',2.5);
hold on;plot(t,smooth(flux(:,29)/RC_outside,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,30)/RC_outside,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,31)/RC_outside,num),'linewidth',2.5);
hold on;plot(t,smooth(flux(:,32)/RC_outside,num),'linewidth',2.5);
l4=legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32');
set(l4,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

end