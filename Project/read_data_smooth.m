%% note: The currents of PF coils from the database have positive signs. So in the simulations the signs have to be changed according to the sign of the plasma current (always positive in EFIT).
% In this script the currents in the PF coils are set to negative, according to the experiments have been done for EXL-50
% **Check the value of the experiment data before use: COILS_t**
% read data from the MDSplus server, t_raw(s), PLASMA_t_MDSplus(kA), EXPMP2_t_MDSplus(T), COILS_t_MDSplus(v.s(wb)), BRSP_t_MDSplus(A) 

fs=handles.fs; % sample rate, unit:Hz
MP_num=length(EXPMP2); %number of magnetic probe
MP_T_num=38; %number of poloidal magnetic probe
MP_N_num=25; %number of magnetic probe in Normal direction
COILS_num=length(COILS); %number of flux loops  unit v.s/rad(wb/rad)
PF_num=handles.PF_num; %number of PF coils
PF_turn=22; %turns of PF coils
ISHOT=handles.ISHOT;
ntime = handles.ntime; %number of time slice choosed from the raw data
tstart = handles.tstart;
tend = handles.tend; %unit:ms
% setappdata(0,'verbose',1) % automatic plotting
span=10; %sets the span of the moving average to span for smooth, span is a percentage of the total number of data points, less than or equal to 1; or the sample rate is 100kHz, 10 means doing the moving average for every 0.1ms
% offset=3.0e2; %Now which is turned off. number of points to calclulate the offset of the signals, default: the average of the data in the first 3 seconds
sct('exl50',ISHOT); %refresh the shots

%%%read EXL50 magnetic positions from excel
    %A1 = xlsread('E:\Engineering Simulation\MM_calibration\EXL50_1A_symmetry.xls','MPT');
    %A2 = xlsread('E:\Engineering Simulation\MM_calibration\EXL50_1A_symmetry.xls','MPN');
    %A3 = xlsread('E:\RTEFIT\RTEFIT-developer\Project\EXL50_1A_symmetry.xls','FLUX');
    A3 = xlsread('EXL50_1A_symmetry.xls','FLUX');

%read the experimental data from the MDSplus server
%smooth then interpolate the smoothed data to the selected time slice
t_opt = round(linspace(tstart,tend,ntime)); %select time slice: start from tstart, end at the tend and divided by number of ntime
% time_range=['-15:',num2str(tend/1000.0),':0.01']; %read raw data every 10 millisecond
time_range=[num2str(tstart/1000.0*0.9),':',num2str(tend/1000.0*1.01),':0.05']; %read raw data every 20 millisecond

%read and smooth the plasma current
try
    [PLASMA_t_MDSplus,t_current]=exl50db(ISHOT,'IP',time_range); %plasma current in unit PLASMA_unit(kA), t_signal in unit s, channels from ip01 to ip04 in different toroidal locations
catch
    errordlg('Error: Can not connect to the EXL50 data base','File Error');
    delete(h);
    error('Error: Can not connect to the EXL50 data base');
end
t_current=t_current*1000.0; %unit s to ms

 if ~isempty(find(isnan(PLASMA_t_MDSplus)==1,1)) || ~isempty(find(isempty(PLASMA_t_MDSplus)==1,1))
        error('Error: The plasma current is NaN or empty');
 end

PLASMA_t_smoothed=smooth(t_current,PLASMA_t_MDSplus);
% PLASMA_t_smoothed=smooth(t_current,PLASMA_t_MDSplus,span,'rloess');
% PLASMA_t = (interp1(t_signal,PLASMA_t_smoothed,t_opt')-repmat(sum(PLASMA_t_smoothed(1:offset))/offset,ntime,1))*1000.0;
if handles.Method==1
    PLASMA_t=interp1(t_current,PLASMA_t_smoothed,t_opt')*1000.0*1.13; %plasma current in unit A, the plasma current multiplys by 1.35 (1.13 by Guo)
elseif handles.Method==2
    PLASMA_t=interp1(t_current,PLASMA_t_smoothed,t_opt')*1000.0*1.13; %plasma current in unit A, the plasma current multiplys by 1.35 (1.13 by Guo)
end
clear varname_IN1 varvalue_IN1 PLASMA_t_MDSplus PLASMA_t_smoothed;

%read and smooth the flux loop data
[COILS_t_MDSplus,t_COILS]=exl50dbN(ISHOT,'flux\d*',time_range); %read the experimental data of flux loops from the MDSplus server, unit v.s(wb)
t_COILS=t_COILS*1000.0;%unit s to ms
for i=1:COILS_num
%     COILS_t_smoothed(:,i)=smooth(COILS_t_MDSplus(:,i),span,'rloess');
    COILS_t_smoothed(:,i)=smooth(COILS_t_MDSplus(:,i));
end
[x_COILS,y_COILS]=meshgrid(t_opt,1:COILS_num);
[x_COILS_raw,y_COILS_raw]=meshgrid(t_COILS',1:COILS_num);
% COILS_t = interp2(x_COILS_raw,y_COILS_raw,COILS_t_smoothed',x_COILS,y_COILS)'-repmat(sum(COILS_t_smoothed(1:offset,:))/offset,ntime,1);
COILS_t = interp2(x_COILS_raw,y_COILS_raw,COILS_t_smoothed',x_COILS,y_COILS)';
clear COILS_t_MDSplus x_COILS_raw y_COILS_raw x_COILS y_COILS COILS_t_smoothed t_COILS;

for i=1:ntime
    COILS_t(i,:)=COILS_t(i,:)./A3(:,6)'.*A3(:,7)'./2./pi;
end
COILS_t((isnan(COILS_t)==1)==1)=0.0;

%modify the data by hand (the flux loop data are not accurate enough to do the fitting, so some of the data are adjusted by hand)
COILS_t(:,35)=COILS_t(:,19);
COILS_t(:,21)=COILS_t(:,33);
COILS_t(:,10:1:14)=COILS_t(:,5:-1:1);
% COILS_t(:,39:-1:33)=COILS_t(:,15:1:21);
% COILS_t(:,15:1:21)=COILS_t(:,39:-1:33);
% % % % % COILS_t(:,9)=COILS_t(:,6);

%check the diagnostic data of flux loop and remove the drift
checkflux_routine;
if handles.Rdrift_on
    COILS_t(:,8:1:14)=COILS_t(:,7:-1:1);
    COILS_t(:,39:-1:33)=COILS_t(:,15:1:21);
    COILS_t(:,32:-1:28)=COILS_t(:,22:1:26);
end

%read and smooth the magnetic probe data
EXPMP2_t=zeros(length(t_opt),MP_num);
% EXPMP2_t_smoothed=zeros(length(t_signal),MP_num);
% for i=14:MP_T_num
%     EXPMP2_t_smoothed(:,i)=smooth(exl50db(ISHOT,['MP',num2str(i,'%03d'),'T'],time_range));  %read the experimental data of poloidal magnetic probes from the MDSplus server
% end
% for i=MP_T_num+1:MP_num
%     EXPMP2_t_smoothed(:,i)=smooth(exl50db(ISHOT,['MP',num2str(i-MP_T_num,'%03d'),'N'],time_range));  %read the experimental data of Normal direction magnetic probe from the MDSplus server
% end
% [x_EXPMP2,y_EXPMP2]=meshgrid(t_opt,1:MP_num);
% [x_EXPMP2_raw,y_EXPMP2_raw]=meshgrid(t_signal,1:MP_num);
% % EXPMP2_t = interp2(x_EXPMP2_raw,y_EXPMP2_raw,EXPMP2_t_smoothed',x_EXPMP2,y_EXPMP2)'-repmat(sum(EXPMP2_t_smoothed(1:offset,:))/offset,ntime,1); %for the interp2, the value of x coordinate should change for each row, the value of y coordinate should change for each column
% EXPMP2_t = interp2(x_EXPMP2_raw,y_EXPMP2_raw,EXPMP2_t_smoothed',x_EXPMP2,y_EXPMP2)'; %for the interp2, the value of x coordinate should change for each row, the value of y coordinate should change for each column
% EXPMP2_t(find(isnan(EXPMP2_t)==1))=0;
% clear EXPMP2_t_smoothed x_EXPMP2_raw y_EXPMP2_raw x_EXPMP2 y_EXPMP2;

%read and smooth the PF current
if handles.Method==1
    BRSP_length=6; %for the option of free boundary EFIT
elseif handles.Method==2
    BRSP_length=20; %for the option of Fitting
end
BRSP_t_smoothed=ones(length(t_current),1)*BRSP(1:BRSP_length);
% for i=1:2, the currents in the PF1-2, PF3-4 and PF5-6 are in series connection, respectively.
% IPF1-6 is the current in PF1-6, respectively, For example, the PF1 is the same as the PF2.
BRSP_t_smoothed(:,1:2)=-repmat(smooth(exl50db(ISHOT,'IPF1',time_range)),1,2); %read the experimental data of PF current from the MDSplus server
BRSP_t_smoothed(:,3:4)=-repmat(smooth(exl50db(ISHOT,'IPF3',time_range)),1,2); %read the experimental data of PF current from the MDSplus server
% BRSP_t_smoothed(:,5:6)=-repmat(smooth(exl50db(ISHOT,'IPF5',time_range)),1,2); %read the experimental data of PF current from the MDSplus server
BRSP_t_smoothed(:,5:6)=0.0; %read the experimental data of PF current from the MDSplus server
[x_BRSP,y_BRSP]=meshgrid(t_opt,1:BRSP_length);
[x_BRSP_raw,y_BRSP_raw]=meshgrid(t_current',1:BRSP_length);
BRSP_t = interp2(x_BRSP_raw,y_BRSP_raw,BRSP_t_smoothed',x_BRSP,y_BRSP)';
% BRSP_t(:,1:PF_num) =-abs(BRSP_t(:,1:PF_num)-repmat(sum(BRSP_t_smoothed(1:offset,1:PF_num))/offset,ntime,1));
BRSP_t((isnan(BRSP_t)==1)==1)=0.0;
clear BRSP_t_smoothed x_BRSP_raw y_BRSP_raw x_BRSP y_BRSP t_signal BRSP_length;


save('data.mat','t_opt','PLASMA_t','COILS_t','EXPMP2_t','BRSP_t');%%Error of the Magnetic probe
% load('data.mat');