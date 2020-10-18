shotnum=331;
[ipf1,~]=exl50db(shotnum,'IPF1','-15:10:0.01');
[ipf2,t]=exl50db(shotnum,'IPF2','-15:10:0.01');
[itf,~]=exl50db(shotnum,'IPF2','-15:10:0.01'); %'ITF02'
[ip315,~]=exl50db(shotnum,'ip','-15:10:0.01');
figure;hg=stackplot({{t/1e3,itf,'TF(kA)'},{t/1e3,ipf1,'PF1-4(A)'},{t/1e3,ipf2,'PF5-6(A)'},{t/1e3,ip315,'Ip(kA)'}}, [], [], [], ['Shot',num2str(shotnum)], 'Time (s)');
xlim([0,10]);
%%
[flux,t_flux]=exl50dbN(shotnum,'flux\d*','0:10:0.01');
%%
t=t_flux(:,1);
figure;stackplot({{t/1e3,flux(:,1),'$\phi$'}}, [], [], [], 'Flux-inner', 'Time (s)');
hold on;plot(t/1e3,flux(:,2),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,3),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,4),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,5),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,6),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,7),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,8),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,9),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,10),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,11),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,12),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,13),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,14),'linewidth',2.5);
legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14')

figure;stackplot({{t/1e3,flux(:,15),'$\phi$'}}, [], [], [], 'Flux_Top', 'Time (s)');
hold on;plot(t/1e3,flux(:,16),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,17),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,18),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,19),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,20),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,21),'linewidth',2.5);
legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21')

figure;stackplot({{t/1e3,flux(:,33),'$\phi$'}}, [], [], [], 'Flux-down', 'Time (s)');
hold on;plot(t/1e3,-flux(:,34),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,35),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,36),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,37),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,38),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,39),'linewidth',2.5);
legend('Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39')


figure;stackplot({{t/1e3,flux(:,22),'$\phi$'}}, [], [], [], 'Flux-Out', 'Time (s)');
hold on;plot(t/1e3,-flux(:,23),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,24),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,25),'linewidth',2.5);
hold on;plot(t/1e3,flux(:,26),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,28),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,29),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,30),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,31),'linewidth',2.5);
hold on;plot(t/1e3,-flux(:,32),'linewidth',2.5);
legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32')
%%