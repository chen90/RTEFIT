%% create fake data by adding random error to the simulated data from a EFIT equilibrim simulation
% only use, if no data from the MDSplus available

fs=handles.fs; % sample rate, unit:Hz
t_length=handles.t_length; %length of the signal, unit: s
t_signal=0:(t_length*1000.0)/(fs*t_length-1):t_length*1000.0; %time seires of the signals, unit: ms, linspace(0,t_length*1000.0,fs*t_length);
MP_num=length(EXPMP2); %number of magnetic probe  unit T
COILS_num=length(COILS); %number of flux loops  unit v.s/rad(wb/rad)
PF_num=handles.PF_num; %number of PF coils
ISHOT=handles.ISHOT;
ntime = handles.ntime; %number of time slice choosed from the raw data
Error_Relative=handles.Error_Relative;
Error_Relative_MP = Error_Relative; %relative error for MP
Error_Relative_COIL = Error_Relative; %relative error for flux loops
Error_Relative_PF = Error_Relative; %relative error for PF coils

%create the raw data: add the error to data with standard normal distribution
%interpolate raw data to the selected time slice
t_opt = round(linspace(handles.tstart,handles.tend,ntime)); %select time slice: start from tstart, end at the tend and divided by number of ntime

PLASMA_t_fake=ones(length(t_signal),1)*PLASMA;
PLASMA_t_fake=PLASMA_t_fake+PLASMA_t_fake.*randn(length(t_signal),length(PLASMA)).*Error_Relative_PF;
PLASMA_t = interp1(t_signal,PLASMA_t_fake',t_opt);
clear varname_IN1 varvalue_IN1 PLASMA_t_fake;

EXPMP2_t_fake=ones(length(t_signal),1)*EXPMP2;
EXPMP2_t_fake=EXPMP2_t_fake+EXPMP2_t_fake.*randn(length(t_signal),MP_num).*Error_Relative_MP;
[x_EXPMP2,y_EXPMP2]=meshgrid(t_opt,1:MP_num);
[x_EXPMP2_raw,y_EXPMP2_raw]=meshgrid(t_signal,1:MP_num);
EXPMP2_t = interp2(x_EXPMP2_raw,y_EXPMP2_raw,EXPMP2_t_fake',x_EXPMP2,y_EXPMP2)'; %for the interp2, the value of x coordinate should change for each row, the value of y coordinate should change for each column
clear EXPMP2_t_fake x_EXPMP2_raw y_EXPMP2_raw x_EXPMP2 y_EXPMP2;

COILS_t_fake=ones(length(t_signal),1)*COILS;
COILS_t_fake=COILS_t_fake+COILS_t_fake.*randn(length(t_signal),COILS_num).*Error_Relative_COIL;
[x_COILS,y_COILS]=meshgrid(t_opt,1:COILS_num);
[x_COILS_raw,y_COILS_raw]=meshgrid(t_signal,1:COILS_num);
COILS_t = interp2(x_COILS_raw,y_COILS_raw,COILS_t_fake',x_COILS,y_COILS)';
clear COILS_t_fake x_COILS_raw y_COILS_raw x_COILS y_COILS;

BRSP_length=20;
BRSP_t_fake_tmp=ones(length(t_signal),1)*BRSP(1:BRSP_length);
BRSP_t_fake=BRSP_t_fake_tmp+BRSP_t_fake_tmp.*randn(length(t_signal),BRSP_length).*Error_Relative_PF;
BRSP_t_fake(:,PF_num+1:end)=BRSP_t_fake_tmp(:,PF_num+1:end);
[x_BRSP,y_BRSP]=meshgrid(t_opt,1:BRSP_length);
[x_BRSP_raw,y_BRSP_raw]=meshgrid(t_signal,1:BRSP_length);
BRSP_t = interp2(x_BRSP_raw,y_BRSP_raw,BRSP_t_fake',x_BRSP,y_BRSP)';
clear BRSP_t_fake_tmp BRSP_t_fake x_BRSP_raw y_BRSP_raw Error_Relative_PF Error_Relative_COIL Error_Relative_MP x_BRSP y_BRSP t_signal BRSP_length;

% % % % %create the raw data: no longer used, due to the slow running speed
% % % % %add the error to data with standard normal distribution
% % % % TIME_t_fake = 1:1:ntime;
% % % % for i=1:ntime
% % % %     EXPMP2_t_fake(i,:)=EXPMP2+EXPMP2.*randn(MP_num,1)'.*Error_Relative_MP;
% % % %     COILS_t_fake(i,:)=COILS+COILS.*randn(COILS_num,1)'.*Error_Relative_COIL;
% % % %     BRSP_t_fake(i,:)=BRSP+BRSP.*randn(length(BRSP),1)'.*Error_Relative_PF;
% % % %     BRSP_t_fake(i,PF_num+1:end)=BRSP(PF_num+1:end);
% % % %     PLASMA_t_fake(i)=PLASMA+PLASMA.*randn(length(PLASMA),1)'.*Error_Relative_PF;
% % % % end