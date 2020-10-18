%%plotting countour 
contour(xx,yy,flux_efit',80,'w');

hold on;
flux_surface=flux_efit'-sibry;
surf(xx,yy,min(min(flux_efit))*ones(size(xx))*1.05,ones(size(xx))*1.0);
surf(xx,yy,flux_surface);
% colormap(flipud(jet));
colormap(flipud(hot));
% colormap(hot);
% colorbar;
% colorbar('Ticks',[],'TickLabels',[]);
shading interp;
xlabel('R(m)','fontsize',18,'fontname','times');
ylabel('Z(m)','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times')
set(gca,'linewidth',2)
grid off
box off
% axis equal;
axis([0 2.2 -2.05 2.05]);
view([0,0,1]);
zlim([min(min(flux_efit))*1.1,0]);
caxis([min(min(flux_efit))*1.1,0]);

%%%plot outline/ flux surface/ coils/ probes after Cal_density.m is called
hold on; 
d=(sibry-simag)/20;
v=-sibry:d:-simag;
v(length(v)+1)=-sibry;
for i=1:length(v)
    [cc,hh]=contour(xx,yy,-flux_efit',[v(i),v(i)]);
    
    if max(size(cc)<1)
        continue
    end
    
    lentmp=size(cc);
    count=1:3;
    sumtmp=1;
    for j=1:count(end)
        if sumtmp>lentmp(2)
              continue
        end 
         count(j)=cc(2,sumtmp);
         tmpx=cc(1,sumtmp+1:sumtmp+count(j));
         tmpy=cc(2,sumtmp+1:sumtmp+count(j));
         tmpz=tmpx*0.0+max(max(Nxy));
         if i==length(v)
            plot3(tmpx,tmpy,tmpz,'r','linewidth',3.0);
         else
            plot3(tmpx,tmpy,tmpz,'color',[0.5 0.5 0.5],'linewidth',0.01);
         end
         sumtmp=sumtmp+count(j)+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%divertor configure%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [xxtmp,yytmp]=meshgrid(x,y);
% % [cc,hh]=contour(xxtmp,yytmp,-flux_efit',[-sibry,-sibry],'w');    
% %     lentmp=size(cc);
% %     count=1:3;
% %     sumtmp=1;
% %     for j=1:count(end)
% %         if sumtmp>lentmp(2)
% %               continue
% %         end 
% %          count(j)=cc(2,sumtmp);
% %          tmpx=cc(1,sumtmp+1:sumtmp+count(j));
% %          tmpy=cc(2,sumtmp+1:sumtmp+count(j));
% %          tmpz=tmpx*0.0+max(max(Nxy));
% %          plot3(tmpx,tmpy,tmpz,'r','linewidth',3.0);
% %          sumtmp=sumtmp+count(j)+1
% %     end
% %   clear  xxtmp yytmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%limiter configure%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot3(rbbbs,zbbbs,rbbbs*0.0+max(max(Nxy)),'b','linewidth',4.0);
plot3(rbbbs,zbbbs,rbbbs*0.0+max(max(Nxy)),'color',[0.8 0 0],'linewidth',3.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vv=[0.158 1.655 -1.405 1.405];  %%%the dimensions of the VV
plot3([vv(1) vv(2)],[vv(4) vv(4)],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0.667 0.667 1],'linewidth',2.5);
plot3([vv(2) vv(2)],[vv(4) vv(3)],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0.667 0.667 1],'linewidth',2.5)
plot3([vv(2) vv(1)],[vv(3) vv(3)],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0.667 0.667 1],'linewidth',2.5)
plot3([vv(1) vv(1)],[vv(3) vv(4)],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0.667 0.667 1],'linewidth',2.5)

limz(1:length(rlim))=tmpz(1);
 plot3(rlim,zlimit,limz,'color',[0.5 0.5 0.5],'linewidth',2);

position_coil=[
    0.963	1.908	0.144	0.215	0 4 6 24 738
    2.1065	1.335	0.144	0.215	0 4 6 24 738
    2.1065	0.445	0.144	0.215	0 4 6 24 738
    0.963	-1.908	0.144	0.215	0 4 6 24 738
    2.1065	-1.335	0.144	0.215	0 4 6 24 738
    2.1065	-0.445	0.144	0.215	0 4 6 24 738
    ];
for i=1:length(position_coil(:,1))
    plot3([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)-position_coil(i,4)/2],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0 0 1],'linewidth',1.5)
    hold on
    plot3([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)+position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0 0 1],'linewidth',1.5)
    plot3([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)-position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0 0 1],'linewidth',1.5)
    plot3([position_coil(i,1)+position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0 0 1],'linewidth',1.5)
    plot3([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)-position_coil(i,4)/2 position_coil(i,2)+position_coil(i,4)/2],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0 0 1],'linewidth',1.5)
    plot3([position_coil(i,1)-position_coil(i,3)/2 position_coil(i,1)+position_coil(i,3)/2],[position_coil(i,2)+position_coil(i,4)/2 position_coil(i,2)-position_coil(i,4)/2],[max(max(-flux_efit))+0.001 max(max(-flux_efit))+0.001],'color',[0 0 1],'linewidth',1.5)
end
