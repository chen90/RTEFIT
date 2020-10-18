%% routine for the flux loop data check
%figure 1-3 flux vs time
%figure 4 flux vs flux loop No

%% constant and control flags
symmetry_left=1; %0/1,if symmetry_left=0, the data on the right are used (FLUX028-39) and the data on the left are the copy of the right. Vice versa (FLUX015-026).
Ipdome_smooth=0; %0/1,if Ipdome_smooth=0, the flux loop data FLUX0022-32 will keep the symmetric raw data; otherwise 1, smooth the data by hand.
Opt_flux=1; %0/1, Optimize the flux loop data near the central rod, if Opt_flux=1, the data are optimized near the central rod; otherwise 0, use the raw data.
COILS_t_tmp=COILS_t*2*pi; %unit Wb
t=t_opt/1000.0; %time seires of the signals, unit: s
fluxloop_No=1:1:size(COILS_t_tmp,2);

%% remove the drift and check the flux loop data
if handles.Rdrift_on
    if handles.ntime==1 %only one time slice ntime=1
        figure('Name','Remove drift from Flux loop signal');
        set(gcf,'units','normalized','position',[0.05 0.1 0.9 0.6]);
        subplot(1,2,1); % figure flux vs flux loop No
        j=1;
        for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
            hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
            l{j}=['t=',num2str(handles.tstart/1000.0),' s'];
            j=j+1;
        end
        xlabel('$Loop-No$','interpreter','Latex');
        ylabel('$\rm\phi (\textit{Wb})$','interpreter','Latex'); %unit Wb
        title('Flux Raw Data')
        set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
        legend(l,'location','SouthWest','fontsize',16);
        
       %% Optimize the flux loop data near the central rod
        if Opt_flux
            % interpolate the data for the flux loop 2-5
            for i=1:length(t)
                COILS_t_tmp(i,2:5) = interp1(fluxloop_No([1,6:end]),COILS_t_tmp([1,6:end]),2:5);
            end 
            COILS_t_tmp(:,8:1:14)=COILS_t_tmp(:,7:-1:1);

            % interpolate the data for the flux loop 7 and 8    
            for i=1:length(t)
                COILS_t_tmp(i,7:8) = interp1(fluxloop_No([1:6,9:end]),COILS_t_tmp([1:6,9:end]),7:8,'spline');
                COILS_t_tmp(i,2) =  COILS_t_tmp(i,1)*0.5;
                COILS_t_tmp(i,3) =  COILS_t_tmp(i,3)*0.2;
                COILS_t_tmp(i,4) =  COILS_t_tmp(i,4)*0.7;
                COILS_t_tmp(i,5) =  COILS_t_tmp(i,5)*0.9;
                coef_tmp=polyfit(1:6,COILS_t_tmp(i,1:6),2);
                COILS_t_tmp(i,2:5)=polyval(coef_tmp,2:5);
            end 
        end %Opt_flux
        
       %% linear fitting of the data to remove the drift from the raw data
        if symmetry_left==1
            if ISHOT<2190 %when the shot number below 2190, remove the offset of some unreasonable signals
                for i=1:length(t)
                    flux_offset=abs(COILS_t_tmp(i,15)-(COILS_t_tmp(i,14)+COILS_t_tmp(i,16))/2);
                    flux_offset = flux_offset+abs(COILS_t_tmp(i,17)-(COILS_t_tmp(i,16)+COILS_t_tmp(i,18))/2);
                    COILS_t_tmp(i,17:end)=COILS_t_tmp(i,17:end)-flux_offset;
                end 

                % interpolate the data for the flux loop 15-21  
                for i=1:length(t)
                    COILS_t_tmp(i,15:21) = interp1(fluxloop_No([1:14,22]),COILS_t_tmp([1:14,22]),15:21,'spline');
                end
            end

            % set the flux loop data completely symmetric
            COILS_t_tmp(:,8:1:14)=COILS_t_tmp(:,7:-1:1);
            COILS_t_tmp(:,39:-1:33)=COILS_t_tmp(:,15:1:21);
            COILS_t_tmp(:,32:-1:28)=COILS_t_tmp(:,22:1:26);

        else %symmetry_left
            if ISHOT<2190 %when the shot number below 2190, remove the offset of some unreasonable signals
                % interpolate the data for the flux loop 33-37  
                for i=1:length(t)
                    COILS_t_tmp(i,33:37) = interp1(fluxloop_No([32,38:39]),COILS_t_tmp([32,38:39]),33:37,'spline');
                end 
            end

                % set the flux loop data completely symmetric
                COILS_t_tmp(:,8:1:14)=COILS_t_tmp(:,7:-1:1);
                COILS_t_tmp(:,15:1:21)=COILS_t_tmp(:,39:-1:33);
                COILS_t_tmp(:,22:1:26)=COILS_t_tmp(:,32:-1:28);

        end %symmetry_left

        % interpolate the data for the flux loop 27    
        for i=1:length(t)
            COILS_t_tmp(i,27) = interp1(fluxloop_No([1:26,28:end]),COILS_t_tmp([1:26,28:end]),27,'spline');
        end

        if handles.Smooth==1 %smooth the symmetric data for the flux loops
            for i=1:length(t)
                COILS_t_tmp(i,22:32)=smooth(COILS_t_tmp(i,22:32));
            end
        end

        if Ipdome_smooth==1  %0/1,if Ipdome_smooth=0, the flux loop data FLUX0022-32 will keep the symmetric raw data; otherwise 1, smooth the data by hand.
            x=0:1:10;
            y=gaussmf(x,[2 5]);
            y(2)=y(2)*1.4;
            y(3)=y(3)*1.35;
            y(4)=y(4)*1.1;
            y(5)=y(5)*0.98;
            y(7:11)=y(5:-1:1);
            y=y-min(y);
            y=y/max(y);
            COILS_t_tmp(:,22:32)=repmat(COILS_t_tmp(:,22),1,11).*(-1.0/1.5.*repmat(y,length(t),1)+1);
        end

        % figure flux vs flux loop No
        subplot(1,2,2);
        j=1;
        for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
            hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
            l{j}=['t=',num2str(handles.tstart/1000.0),' s'];
            j=j+1;
        end
        xlabel('$Loop-No$','interpreter','Latex');
        ylabel('$\rm\phi (\textit{Wb})$','interpreter','Latex'); %unit Wb
        title('Flux Data Processed')
        set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
        % legend(l,'location','northeastOutside','fontsize',10);
        legend(l,'location','SouthWest','fontsize',16);
        
    else %more than one time slice
      
        %% plot flux loops raw data
        figure('Name','Remove drift from Flux loop signal');
        set(gcf,'units','normalized','position',[0.01 0.04 0.98 0.88]);
        subplot(2,4,1);
        stackplot({{t, COILS_t_tmp(:,1),'$\phi(Wb)$'}}, [], [], [], ['Inner Flux'], 'Time (s)');%figure 1
        hold on;plot(t, COILS_t_tmp(:,2),'linewidth',2.5);
        hold on;plot(t, COILS_t_tmp(:,3),'linewidth',2.5);
        hold on;plot(t, COILS_t_tmp(:,4),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,5),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,6),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,7),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,8),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,9),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,10),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,11),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,12),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,13),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,14),'linewidth',2.5);
        l2=legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14');
        set(l2,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);
        % xlim([-2,10]);

        subplot(2,4,2);stackplot({{t,COILS_t_tmp(:,15),''}}, [], [], [], ['Up & Down Flux'], 'Time (s)');%figure 2
        hold on;plot(t,COILS_t_tmp(:,16),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,17),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,18),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,19),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,20),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,21),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,33),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,34),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,35),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,36),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,37),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,38),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,39),'--','linewidth',2.5);
        l3=legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21','Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39');
        set(l3,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

        subplot(2,4,3);stackplot({{t,COILS_t_tmp(:,22),''}}, [], [], [], ['OutSide Flux'], 'Time (s)');%figure 3
        hold on;plot(t,COILS_t_tmp(:,23),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,24),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,25),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,26),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,28),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,29),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,30),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,31),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,32),'linewidth',2.5);
        l4=legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32');
        set(l4,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

        % figure 4 flux vs flux loop No
        subplot(2,4,4);
        j=1;
        for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
            hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
            l{j}=['t=',num2str(t(i)),' s'];
            j=j+1;
        end
        xlabel('$Loop-No$','interpreter','Latex');
        % ylabel('$\rm\phi (Wb)$','interpreter','Latex'); %unit Wb
        title('Flux Raw Data')
        set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
        legend(l,'location','SouthWest','fontsize',10);
        
   %% Optimize the flux loop data near the central rod
        if Opt_flux
            % interpolate the data for the flux loop 2-5
            for i=1:length(t)
                COILS_t_tmp(i,2:5) = interp1(fluxloop_No([1,6:end]),COILS_t_tmp(i,[1,6:end]),2:5);
            end 
            COILS_t_tmp(:,8:1:14)=COILS_t_tmp(:,7:-1:1);

            % interpolate the data for the flux loop 7 and 8    
            for i=1:length(t)
                COILS_t_tmp(i,7:8) = interp1(fluxloop_No([1:6,9:end]),COILS_t_tmp(i,[1:6,9:end]),7:8,'spline');
                COILS_t_tmp(i,2) =  COILS_t_tmp(i,1)*0.5;
                COILS_t_tmp(i,3) =  COILS_t_tmp(i,3)*0.2;
                COILS_t_tmp(i,4) =  COILS_t_tmp(i,4)*0.7;
                COILS_t_tmp(i,5) =  COILS_t_tmp(i,5)*0.9;
                coef_tmp=polyfit(1:6,COILS_t_tmp(i,1:6),2);
                COILS_t_tmp(i,2:5)=polyval(coef_tmp,2:5);
            end 
        end %Opt_flux
  
    %% linear fiting of the data to remove the drift from the raw data
        %after the slop was removed, the data can't be used as well, the
        %distribution of the flux loops is wrong compared to the data from the
        %equlibrium calculation. How to deal with the drift is not clear yet!
        % %     if t(1)<handles.t_start_Rdrift && t(end)>handles.t_end_Rdrift
        % %         t_start_Rdrift=find(t==max(t(t<handles.t_start_Rdrift)));
        % %         t_end_Rdrift=find(t==min(t(t>handles.t_end_Rdrift)));
        % %     elseif t(1)>handles.t_start_Rdrift && t(1)<handles.t_end_Rdrift && t(end)>handles.t_end_Rdrift
        % %         t_start_Rdrift=1;
        % %         t_end_Rdrift=find(t==min(t(t>handles.t_end_Rdrift)));
        % %     elseif t(1)<handles.t_start_Rdrift && t(end)<handles.t_end_Rdrift && t(end)>handles.t_start_Rdrift
        % %         t_start_Rdrift=find(t==max(t(t<handles.t_start_Rdrift)));
        % %         t_end_Rdrift=length(t);
        % %     else
        % %         t_start_Rdrift=1;
        % %         t_end_Rdrift=length(t);
        % %     end
        % %     
        % %     for i=1:39
        % %         coef_tmp=polyfit(t(t_start_Rdrift:t_end_Rdrift),COILS_t_tmp(t_start_Rdrift:t_end_Rdrift,i),1);
        % %         COILS_t_tmp(:,i)=COILS_t_tmp(:,i)-coef_tmp(1).*t;
        % % %         if t(1)<=0 && t(end)>0
        % % %             flux_offset = interp1(t,COILS_t_tmp(:,i),0);
        % % %         else
        % % %             flux_offset=0.0;
        % % %         end
        % % %         COILS_t_tmp(:,i)=COILS_t_tmp(:,i)-flux_offset;
        % %     end
                
        if symmetry_left==1
            if ISHOT<2190 %when the shot number below 2190, remove the offset of some unreasonable signals
                for i=1:length(t)
                    flux_offset=abs(COILS_t_tmp(i,15)-(COILS_t_tmp(i,14)+COILS_t_tmp(i,16))/2);
                    flux_offset = flux_offset+abs(COILS_t_tmp(i,17)-(COILS_t_tmp(i,16)+COILS_t_tmp(i,18))/2);
                    COILS_t_tmp(i,17:end)=COILS_t_tmp(i,17:end)-flux_offset;
                end 

                % interpolate the data for the flux loop 15-21  
                for i=1:length(t)
                    COILS_t_tmp(i,15:21) = interp1(fluxloop_No([1:14,22]),COILS_t_tmp(i,[1:14,22]),15:21,'spline');
                end
            end

            % set the flux loop data completely symmetric
            COILS_t_tmp(:,8:1:14)=COILS_t_tmp(:,7:-1:1);
            COILS_t_tmp(:,39:-1:33)=COILS_t_tmp(:,15:1:21);
            COILS_t_tmp(:,32:-1:28)=COILS_t_tmp(:,22:1:26);

       else %symmetry_left
            if ISHOT<2190 %when the shot number below 2190, remove the offset of some unreasonable signals
                % interpolate the data for the flux loop 33-37  
                for i=1:length(t)
                    COILS_t_tmp(i,33:37) = interp1(fluxloop_No([32,38:39]),COILS_t_tmp(i,[32,38:39]),33:37,'spline');
                end 
            end

                % set the flux loop data completely symmetric
                COILS_t_tmp(:,8:1:14)=COILS_t_tmp(:,7:-1:1);
                COILS_t_tmp(:,15:1:21)=COILS_t_tmp(:,39:-1:33);
                COILS_t_tmp(:,22:1:26)=COILS_t_tmp(:,32:-1:28);
        end %symmetry_left
      
        % interpolate the data for the flux loop 27    
        for i=1:length(t)
            COILS_t_tmp(i,27) = interp1(fluxloop_No([1:26,28:end]),COILS_t_tmp(i,[1:26,28:end]),27,'spline');
        end

        if handles.Smooth==1 %smooth the symmetric data for the flux loops
            for i=1:length(t)
                COILS_t_tmp(i,22:32)=smooth(COILS_t_tmp(i,22:32));
            end
        end

        if Ipdome_smooth==1  %0/1,if Ipdome_smooth=0, the flux loop data FLUX0022-32 will keep the symmetric raw data; otherwise 1, smooth the data by hand.
            x=0:1:10;
            y=gaussmf(x,[2 5]);
            y(2)=y(2)*1.4;
            y(3)=y(3)*1.35;
            y(4)=y(4)*1.1;
            y(5)=y(5)*0.98;
            y(7:11)=y(5:-1:1);
            y=y-min(y);
            y=y/max(y);
            COILS_t_tmp(:,22:32)=repmat(COILS_t_tmp(:,22),1,11).*(-1.0/1.5.*repmat(y,length(t),1)+1);
        end
     
       %% figure 1-3 flux vs time
        subplot(2,4,5);
        stackplot({{t, COILS_t_tmp(:,1),'$\phi(Wb)$'}}, [], [], [], ['Inner Flux'], 'Time (s)');%figure 1
        hold on;plot(t, COILS_t_tmp(:,2),'linewidth',2.5);
        hold on;plot(t, COILS_t_tmp(:,3),'linewidth',2.5);
        hold on;plot(t, COILS_t_tmp(:,4),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,5),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,6),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,7),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,8),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,9),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,10),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,11),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,12),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,13),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,14),'linewidth',2.5);
        l2=legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14');
        set(l2,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);
        % xlim([-2,10]);

        subplot(2,4,6);stackplot({{t,COILS_t_tmp(:,15),''}}, [], [], [], ['Up & Down Flux'], 'Time (s)');%figure 2
        hold on;plot(t,COILS_t_tmp(:,16),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,17),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,18),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,19),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,20),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,21),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,33),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,34),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,35),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,36),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,37),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,38),'--','linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,39),'--','linewidth',2.5);
        l3=legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21','Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39');
        set(l3,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

        subplot(2,4,7);stackplot({{t,COILS_t_tmp(:,22),''}}, [], [], [], ['OutSide Flux'], 'Time (s)');%figure 3
        hold on;plot(t,COILS_t_tmp(:,23),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,24),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,25),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,26),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,28),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,29),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,30),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,31),'linewidth',2.5);
        hold on;plot(t,COILS_t_tmp(:,32),'linewidth',2.5);
        l4=legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32');
        set(l4,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

        % figure flux vs flux loop No
        subplot(2,4,8);
        j=1;
        for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
            hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
            l{j}=['t=',num2str(t(i)),' s'];
            j=j+1;
        end
        xlabel('$Loop-No$','interpreter','Latex');
        % ylabel('$\rm\phi (Wb)$','interpreter','Latex'); %unit Wb
        title('Flux data processed')
        set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
        legend(l,'location','SouthWest','fontsize',10);
    end %handles.ntime==1
    
else %Rdrift_on=false
    
    %% Not remove the drift and check the flux loop data, if only one time slice
     if handles.ntime==1 %only one time slice ntime=1
         if handles.Smooth==1 %smooth the raw data for all the flux loops
            figure('Name','Smooth Flux loop signals');
            set(gcf,'units','normalized','position',[0.05 0.1 0.9 0.6]);
            subplot(1,2,1);% figure flux vs flux loop No
            j=1;
            for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
                hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
                l{j}=['t=',num2str(handles.tstart/1000.0),' s'];
                j=j+1;
            end
            xlabel('$Loop-No$','interpreter','Latex');
            ylabel('$\rm\phi (\textit{Wb})$','interpreter','Latex'); %unit Wb
            title('Flux Raw Data')
            set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
            legend(l,'location','SouthWest','fontsize',16);
            
            COILS_t_tmp(:,13:1:14)=COILS_t_tmp(:,2:-1:1);
            COILS_t_tmp=smooth(COILS_t_tmp)'; %unit Wb %smooth the raw data for all the flux loops

            subplot(1,2,2);% figure flux vs flux loop No
            j=1;
            for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
                hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
                l{j}=['t=',num2str(handles.tstart/1000.0),' s'];
                j=j+1;
            end
            xlabel('$Loop-No$','interpreter','Latex');
            ylabel('$\rm\phi (\textit{Wb})$','interpreter','Latex'); %unit Wb
            title('Flux Data Processed')
            set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
            legend(l,'location','SouthWest','fontsize',16);
         
         else %handles.Smooth==1
         
            %plot the flux loops raw data    
            figure('Name','Original flux loop signal');% figure flux vs flux loop No
            set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6]);
            j=1;
            for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
                hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
                l{j}=['t=',num2str(t(i)),' s'];
                j=j+1;
            end
            xlabel('$Loop-No$','interpreter','Latex');
            ylabel('$\rm\phi (\textit{Wb})$','interpreter','Latex');
            title('Flux Raw Data')
            set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 18, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
            legend(l,'location','SouthWest','fontsize',18);
         end %handles.Smooth==1
         
     else %handles.ntime==1, more than one time slice
         
     %% more than one time slice, smooth the raw data for all the flux loops
        if handles.Smooth==1 
            figure('Name','Smooth Flux loop signals');
            set(gcf,'units','normalized','position',[0.01 0.04 0.98 0.88]);
            subplot(2,4,1);
            stackplot({{t, COILS_t_tmp(:,1),'$\phi(Wb)$'}}, [], [], [], ['Inner Flux'], 'Time (s)');%figure 1
            hold on;plot(t, COILS_t_tmp(:,2),'linewidth',2.5);
            hold on;plot(t, COILS_t_tmp(:,3),'linewidth',2.5);
            hold on;plot(t, COILS_t_tmp(:,4),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,5),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,6),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,7),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,8),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,9),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,10),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,11),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,12),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,13),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,14),'linewidth',2.5);
            l2=legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14');
            set(l2,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);
            % xlim([-2,10]);

            subplot(2,4,2);stackplot({{t,COILS_t_tmp(:,15),''}}, [], [], [], ['Up & Down Flux'], 'Time (s)');%figure 2
            hold on;plot(t,COILS_t_tmp(:,16),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,17),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,18),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,19),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,20),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,21),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,33),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,34),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,35),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,36),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,37),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,38),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,39),'--','linewidth',2.5);
            l3=legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21','Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39');
            set(l3,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

            subplot(2,4,3);stackplot({{t,COILS_t_tmp(:,22),''}}, [], [], [], ['OutSide Flux'], 'Time (s)');%figure 3
            hold on;plot(t,COILS_t_tmp(:,23),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,24),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,25),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,26),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,28),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,29),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,30),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,31),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,32),'linewidth',2.5);
            l4=legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32');
            set(l4,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

            % figure 4 flux vs flux loop No
            subplot(2,4,4);
            j=1;
            for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
                hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
                l{j}=['t=',num2str(t(i)),' s'];
                j=j+1;
            end
            xlabel('$Loop-No$','interpreter','Latex');
            % ylabel('$\rm\phi (Wb)$','interpreter','Latex'); %unit Wb
            title('Flux Raw Data')
            set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
            legend(l,'location','SouthWest','fontsize',10);

            COILS_t_tmp(:,13:1:14)=COILS_t_tmp(:,2:-1:1);
            for i=1:length(t)
                  COILS_t_tmp(i,:)=smooth(COILS_t_tmp(i,:)); %smooth the raw data for all the flux loops
            end
            
          %% figure 1-3 flux vs time
            subplot(2,4,5);
            stackplot({{t, COILS_t_tmp(:,1),'$\phi(Wb)$'}}, [], [], [], ['Inner Flux'], 'Time (s)');%figure 1
            hold on;plot(t, COILS_t_tmp(:,2),'linewidth',2.5);
            hold on;plot(t, COILS_t_tmp(:,3),'linewidth',2.5);
            hold on;plot(t, COILS_t_tmp(:,4),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,5),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,6),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,7),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,8),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,9),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,10),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,11),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,12),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,13),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,14),'linewidth',2.5);
            l2=legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14');
            set(l2,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);
            % xlim([-2,10]);

            subplot(2,4,6);stackplot({{t,COILS_t_tmp(:,15),''}}, [], [], [], ['Up & Down Flux'], 'Time (s)');%figure 2
            hold on;plot(t,COILS_t_tmp(:,16),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,17),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,18),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,19),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,20),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,21),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,33),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,34),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,35),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,36),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,37),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,38),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,39),'--','linewidth',2.5);
            l3=legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21','Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39');
            set(l3,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

            subplot(2,4,7);stackplot({{t,COILS_t_tmp(:,22),''}}, [], [], [], ['OutSide Flux'], 'Time (s)');%figure 3
            hold on;plot(t,COILS_t_tmp(:,23),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,24),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,25),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,26),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,28),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,29),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,30),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,31),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,32),'linewidth',2.5);
            l4=legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32');
            set(l4,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

          % figure flux vs flux loop No
            subplot(2,4,8);
            j=1;
            for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
                hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
                l{j}=['t=',num2str(t(i)),' s'];
                j=j+1;
            end
            xlabel('$Loop-No$','interpreter','Latex');
            % ylabel('$\rm\phi (Wb)$','interpreter','Latex'); %unit Wb
            title('Flux data processed')
            set(gca,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
            legend(l,'location','SouthWest','fontsize',10);
            
        else %handles.Smooth==1
            
           %% Not smooth the raw data for all the flux loops
            %plot the flux loops raw data
            figure('Name','Original flux loop signal');
            set(gcf,'units','normalized','position',[0.01 0.04 0.98 0.88]);
            subplot(2,2,1);
            stackplot({{t, COILS_t_tmp(:,1),'$\phi(Wb)$'}}, [], [], [], ['Inner Flux'], 'Time (s)');%figure 1
            hold on;plot(t, COILS_t_tmp(:,2),'linewidth',2.5);
            hold on;plot(t, COILS_t_tmp(:,3),'linewidth',2.5);
            hold on;plot(t, COILS_t_tmp(:,4),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,5),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,6),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,7),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,8),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,9),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,10),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,11),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,12),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,13),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,14),'linewidth',2.5);
            l2=legend('Flux01','Flux02','Flux03','Flux04','Flux05','Flux06','Flux07','Flux08','Flux09','Flux10','Flux11','Flux12','Flux13','Flux14');
            set(l2,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);
            % xlim([-2,10]);

            subplot(2,2,2);stackplot({{t,COILS_t_tmp(:,15),'$\phi(Wb)$'}}, [], [], [], ['Up & Down Flux'], 'Time (s)');%figure 2
            hold on;plot(t,COILS_t_tmp(:,16),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,17),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,18),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,19),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,20),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,21),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,33),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,34),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,35),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,36),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,37),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,38),'--','linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,39),'--','linewidth',2.5);
            l3=legend('Flux15','Flux16','Flux17','Flux18','Flux19','Flux20','Flux21','Flux33','Flux34','Flux35','Flux36','Flux37','Flux38','Flux39');
            set(l3,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

            subplot(2,2,3);stackplot({{t,COILS_t_tmp(:,22),'$\phi(Wb)$'}}, [], [], [], ['OutSide Flux'], 'Time (s)');%figure 3
            hold on;plot(t,COILS_t_tmp(:,23),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,24),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,25),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,26),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,28),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,29),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,30),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,31),'linewidth',2.5);
            hold on;plot(t,COILS_t_tmp(:,32),'linewidth',2.5);
            l4=legend('Flux22','Flux23','Flux24','Flux25','Flux26','Flux28','Flux29','Flux30','Flux31','Flux32');
            set(l4,'fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 10);

            % figure flux vs flux loop No
            subplot(2,2,4);
            j=1;
            for i=1:ceil(size(COILS_t_tmp,1)/10):size(COILS_t_tmp,1);
                hold on;plot(fluxloop_No,COILS_t_tmp(i,:),':o','LineWidth',1.5,'MarkerSize',6); %unit Wb
                l{j}=['t=',num2str(t(i)),' s'];
                j=j+1;
            end
            xlabel('$Loop-No$','interpreter','Latex');
            ylabel('$\rm\phi (\textit{Wb})$','interpreter','Latex');
            title('Flux Raw Data')
            set(gca,'FontAngle',  'normal','fontname', 'Times New Roman', 'FontWeight', 'normal', 'FontSize', 16, 'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on','ticklength',[0.02 0.02],'Xgrid','off');
            legend(l,'location','SouthWest','fontsize',10);
        end %handles.Smooth==1
     end %handles.ntime==1
end %handles.Rdrift_on

if handles.Systems==1 %Windows operation system
   saveas(gcf,[path_movie,'\shot',num2str(ISHOT,'%06d'),'fluxloop.png']);
elseif handles.Systems==2 %Unix operation system
   saveas(gcf,[path_movie,'/shot',num2str(ISHOT,'%06d'),'fluxloop.png']);
end

% interpolate the data for the flux loop 27
if find(isnan(COILS_t_tmp(1,:))==1)==27
    for i=1:length(t)
        COILS_t_tmp(i,27) = interp1(fluxloop_No([1:26,28:end]),COILS_t_tmp(i,[1:26,28:end]),27,'spline');
    end
end

% save in txt
    %COILS_t_save=COILS_t'; 
    %save([path_movie,'\shot',num2str(ISHOT,'%06d'),'fluxloop_raw.txt'],'COILS_t_save','-ascii','-append'); %unit Wb/rad

    COILS_t=COILS_t_tmp/2/pi; %unit Wb/rad
    
    %COILS_t_save=COILS_t';
    %save([path_movie,'\shot',num2str(ISHOT,'%06d'),'fluxloop_smoothed.txt'],'COILS_t_save','-ascii','-append');