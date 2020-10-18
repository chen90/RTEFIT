function movie(path_gfiles,path_movie,ISHOT,TIME_t_ms,FrameRate,Systems,Play,ax)
%% Movie program of the real time EFIT
% -inputs:
% path_gfiles                       string,  path where the gfiles are saved
% path_movie                       string,  path where the movies are saved
% ISHOT[000000]                    int,  shot number
% TIME_t_ms[00000]                 int,  time slice for the discharge of ISHOT (unit:ms)
% FrameRate                        int,   number of frames per second in avi format
% DelayTime                        Float, Delay time for the movie in GIF format
% Systems                          logical, Operation systems for the local computer:1 or 2，which implies Windows or Unix, respectively
% ax                               cell, %the two axeses from the GUI fig,1/2, 1:outline, 2:surface and parameteres

% -outputs:
% 2 movies which are saved in the path of path_movie: one is in format avi, the other is in format GIF

% Edited by Bin Chen in 2019/06/03
% Contact: chenbino@enn.cn
% ENN Group 1989-2019, all rights reserved.

set(0,'defaultfigurecolor','w');
%% make 2 movies: one is in format avi, the other is in format GIF
%创建avi文件对象
if Systems==1 %Windows operation system
   aviobj = VideoWriter([path_movie,'\EquFitting_shot',num2str(ISHOT,'%06d'),'.avi']); % movie in format avi
elseif Systems==2 %Unix operation system
   aviobj = VideoWriter([path_movie,'/EquFitting_shot',num2str(ISHOT,'%06d'),'.avi']); % movie in format avi
end

aviobj.FrameRate=FrameRate; %number of frames per second
open(aviobj); 

%动画部分代码
% tshow=vpa(TIME_t_s*1000,6);    %change ploting time unit from s to ms
TIME_t_s=roundn(TIME_t_ms*1E-3,-3); %change ploting time unit from ms to s

if Play==0
    fig=figure('Visible','off');
%     fig=figure();
    set(fig,'units','normalized','position',[0.1 0.05 0.6 0.88]);
    ax{1}=gca;
    outline; % plot the outline for the device
    ax1_position=get(ax{1},'position');
    axes('units','normalized','position',ax1_position);
    outline; % plot the outline for the device
    ax{2}=gca;
    info_equ('',0,ax{1});    
else 
    axes(ax{2});
    fig=gcf; %define the plot area
    outline; % plot the outline for the device
end

for k=1:length(TIME_t_ms)
    axes(ax{2});
    if Systems==1 %Windows operation system
           path_gfile=[path_gfiles,'\',num2str(ISHOT,'%06d'),'.',num2str(k,'%05d'),'\g',num2str(ISHOT,'%06d'),'.',num2str(floor(str2num(num2str(TIME_t_ms(k),'%12e'))),'%05i')];
           path_afile=[path_gfiles,'\',num2str(ISHOT,'%06d'),'.',num2str(k,'%05d'),'\a',num2str(ISHOT,'%06d'),'.',num2str(floor(str2num(num2str(TIME_t_ms(k),'%12e'))),'%05i')];
    elseif Systems==2 %Unix operation system
           path_gfile=[path_gfiles,'/',num2str(ISHOT,'%06d'),'.',num2str(k,'%05d'),'/g',num2str(ISHOT,'%06d'),'.',num2str(floor(str2num(num2str(TIME_t_ms(k),'%12e'))),'%05i')];
           path_afile=[path_gfiles,'/',num2str(ISHOT,'%06d'),'.',num2str(k,'%05d'),'/a',num2str(ISHOT,'%06d'),'.',num2str(floor(str2num(num2str(TIME_t_ms(k),'%12e'))),'%05i')];
    end
    
    if ~exist(path_gfile,'file') || ~exist(path_afile,'file') 
        hold off;
        plot(rlimter,zlimter,'color',[0.5 0.5 0.5],'linewidth',2); %if no gfile, then just plot the limiter
        set(gca,'fontweight','bold','fontsize',16,'fontname','times','XTick',[],'YTick',[],'Color', 'None');
        axis equal;
        axis([0 2.2 -2.0 2.0]);
        xdim=sum(get(ax{1},'xlim'));
        ydim=get(ax{1},'ylim');
        title(['SHOT: ',num2str(ISHOT,'%06d'), ', t=', num2str(TIME_t_ms(k),'%05d') '(ms)']);
        %title(['SHOT: ',num2str(ISHOT,'%06d'), ', t=', num2str(TIME_t_ms(k)/1000.0-0.6,'%6.4f') '(s)']);
        text(xdim*1.1,ydim(2)*0.8,'No Reconstruction','horiz','left','color','r','fontsize',16);
        box on;
    else
%         gfile=fluxsurface(path_gfile); %make the movie with lines
        fluxsurface(path_gfile,'contour'); %make the movie with contour colors
        info_equ(path_afile,1,ax{1});
        %title(['SHOT: ',num2str(ISHOT,'%06d'), ', t=' num2str(TIME_t_s(k),'%6.3f') '(s)']);
        title(['SHOT: ',num2str(ISHOT,'%06d'), ', t=' num2str(TIME_t_ms(k),'%05d') '(ms)']);
        %title(['SHOT: ',num2str(ISHOT,'%06d'), ', t=', num2str(TIME_t_ms(k)/1000.0-0.6,'%6.4f') '(s)']);
        axes(ax{1});
    end
   
    %获取当前画面
    F = getframe(fig);
    %加入avi对象中
    writeVideo(aviobj,F);
    %转成gif图片,只能用256色
    im = frame2im(F);
    [I,map] = rgb2ind(im,256);
    %写入 GIF89a 格式文件    
    if k == 1;
        if Systems==1 %Windows operation system
           imwrite(I,map,[path_movie,'\EquFitting_shot',num2str(ISHOT,'%06d'),'.gif'],'GIF', 'Loopcount',inf,'DelayTime',0.1); % movie in format gif
        elseif Systems==2 %Unix operation system
           imwrite(I,map,[path_movie,'/EquFitting_shot',num2str(ISHOT,'%06d'),'.gif'],'GIF', 'Loopcount',inf,'DelayTime',0.1); % movie in format gif
        end
    else
        if Systems==1 %Windows operation system
           imwrite(I,map,[path_movie,'\EquFitting_shot',num2str(ISHOT,'%06d'),'.gif'],'GIF','WriteMode','append','DelayTime',0.1);
        elseif Systems==2 %Unix operation system
           imwrite(I,map,[path_movie,'/EquFitting_shot',num2str(ISHOT,'%06d'),'.gif'],'GIF','WriteMode','append','DelayTime',0.1);
        end
    end

%         clf(fig); %clear the figure window
%         cla(gca); %reset the axes
end
%关闭avi对象
close(aviobj);

end