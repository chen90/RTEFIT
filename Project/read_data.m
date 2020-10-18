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
% offset=3.0e2; %Now which is turned off. number of points to calclulate the offset of the signals, default: the average of the data in the first 3 seconds
% sct('exl50',ISHOT); %refresh the shots

%%%read EXL50 magnetic positions from excel
    %A1 = xlsread('E:\Engineering Simulation\MM_calibration\EXL50_1A_symmetry.xls','MPT');
    %A2 = xlsread('E:\Engineering Simulation\MM_calibration\EXL50_1A_symmetry.xls','MPN');
    %A3 = xlsread('E:\RTEFIT\RTEFIT-developer\Project\EXL50_1A_symmetry.xls','FLUX');
    A3 = xlsread('EXL50_1A_asymmetry.xls','FLUX');
    %fprintf('%10.4f,',A3(:,1)/1000.0);
    
%read the experimental data from the MDSplus server, the raw data will not be smoothed
if handles.ntime==1 %only one time slice ntime=1
    handles.tend=handles.tstart;
    tend=handles.tend;
    set(handles.tend_edit, 'string', handles.tstart/1000.0);
    guidata(hObject, handles);
    ntime=5;
    time_range=[num2str(tstart/1000.0-0.1),':',num2str(tend/1000.0+0.1),':',num2str(((tend-tstart)/1000.0+0.2)/ntime)];
else
    time_range=[num2str(tstart/1000.0),':',num2str(tend/1000.0),':',num2str((tend-tstart)/1000.0/ntime)];
end

%read the plasma current
try
    [PLASMA_t_MDSplus,t_signal]=exl50db(ISHOT,'IP',time_range); %plasma current in unit PLASMA_unit(kA), t_signal in unit s, channels from ip01 to ip04 in different toroidal locations
catch
    errordlg('Error: Can not connect to the EXL50 data base','File Error');
    delete(h);
    error('Error: Can not connect to the EXL50 data base');
end

 if ~isempty(find(isnan(PLASMA_t_MDSplus)==1,1)) || ~isempty(find(isempty(PLASMA_t_MDSplus)==1,1))
      error('Error: The plasma current is NaN or empty');
 end

if handles.Method==1
    PLASMA_t=PLASMA_t_MDSplus(1:ntime)*1000.0; %plasma current in unit A, if want to solve the vacuum field *0
    if ISHOT<2190 %when the shot number below 2190, the plasma current data should multiplys by 1.13 (by Guo)
        PLASMA_t=PLASMA_t*1.13;
    end
elseif handles.Method==2
    PLASMA_t=PLASMA_t_MDSplus(1:ntime)*1000.0; %plasma current in unit A
    if ISHOT<2190 %when the shot number below 2190, the plasma current data should multiplys by 1.13 (by Guo)
        PLASMA_t=PLASMA_t*1.13;
    end
end

t_opt=t_signal(1:ntime)*1000.0; %time seires of the signals, unit: ms
clear PLASMA_t_MDSplus varname_IN1 varvalue_IN1;

%read the flux loop data
COILS_t_MDSplus=exl50dbN(ISHOT,'flux\d*',time_range); %read the experimental data of flux loops from the MDSplus server, unit v.s(wb)
COILS_t = COILS_t_MDSplus(1:ntime,1:COILS_num);
clear COILS_t_MDSplus;

if ISHOT<4816
    rc=[0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.001	0.002	0.002	0.002	0.002	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.01	0.01	0.01	0.01	0.01	0	0.01	0.01	0.01	0.01	0.01	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0	0 ];
    rc(28)=-rc(28);
    rc(1:10)=-rc(1:10);
    rc(22:32)=-rc(22:32);
    rc=-rc;
elseif ISHOT<5018
    rc=[0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.005	0.02	0.02	0.02	0.02	0.02	0.05	0.05	0.05	0.05	0.05	0.05	0.05	0	0.05	0.05	0.05	0.05	0.05	0.05	0.05	0.02	0.02	0.02	0.02	0.02	0.05	0.05];   
    rc(28)=-rc(28);
    rc(17)=-rc(17);
    rc(34)=-rc(34);
elseif ISHOT<5572
    rc=[0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01	0.01];   
    rc(28)=-rc(28);
    rc(17)=-rc(17);
    rc(34)=-rc(34);
elseif ISHOT<5790
    rc=[0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02	0.02];   
    rc(28)=-rc(28);
    rc(17)=-rc(17);
    rc(34)=-rc(34);
else
    rc=[0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.002	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047	0.047];   
    rc(28)=-rc(28);
    rc(17)=-rc(17);
    rc(34)=-rc(34);
end

for i=1:ntime
    COILS_t(i,:)=COILS_t(i,:).*rc(1:39)./2./pi;
end

if find(isnan(COILS_t)==1)
    warning('Warning: Checking the fluxloop data, some of the values are NaN');
end
% % COILS_t(find(isnan(COILS_t)==1))=0.0;
% COILS_t((isnan(COILS_t)==1)==1)=0.0;

%modify the data by hand (the flux loop data are not accurate enough to do the fitting, so some of the data are adjusted by hand)
if ISHOT<2190 %when the shot number below 2190, the flux loops data are not ready, some of the bad signals are adjusted here
    COILS_t(:,35)=COILS_t(:,19);
    COILS_t(:,21)=COILS_t(:,33);
    COILS_t(:,10:1:14)=COILS_t(:,5:-1:1);
end

% %read the magnetic probe data
EXPMP2_t_MDSplus=zeros(length(t_signal),MP_num);
% for i=14:MP_T_num
%     EXPMP2_t_MDSplus(:,i)=exl50db(ISHOT,['MP',num2str(i,'%03d'),'T'],time_range);  %read the experimental data of poloidal magnetic probes from the MDSplus server
% end
% for i=MP_T_num+1:MP_num
%     EXPMP2_t_MDSplus(:,i)=exl50db(ISHOT,['MP',num2str(i-MP_T_num,'%03d'),'N'],time_range);  %read the experimental data of Normal direction magnetic probe from the MDSplus server
% end
EXPMP2_t = EXPMP2_t_MDSplus(1:ntime,:);
% EXPMP2_t((isnan(EXPMP2_t)==1)==1)=0.0;
clear EXPMP2_t_MDSplus;

%read the PF current
if handles.Method==1
    BRSP_length=8; %for the option of free boundary EFIT
elseif handles.Method==2
    BRSP_length=20; %for the option of Fitting
end
BRSP_t_MDSplus=ones(length(t_signal),1)*BRSP(1:BRSP_length);
% for i=1:2, the currents in the PF1-2, PF3-4 and PF5-6 are in series connection, respectively.
% IPF1-6 is the current in PF1-6, respectively, For example, the PF1 is the same as the PF2.
BRSP_t_MDSplus(:,1:2)=-repmat(exl50db(ISHOT,'IPF1',time_range),1,2); %read the experimental data of PF current from the MDSplus server
if ISHOT<2190 %when the shot number below 2190, not use the PF5-6 signals
    BRSP_t_MDSplus(:,3:4)=-repmat(exl50db(ISHOT,'IPF3',time_range),1,2); %read the experimental data of PF current from the MDSplus server
    BRSP_t_MDSplus(:,5:6)=-repmat(exl50db(ISHOT,'IPF6',time_range),1,2); %read the experimental data of PF current from the MDSplus server
    BRSP_t_MDSplus(:,8)=0.0; %read the experimental data of TF current from the MDSplus server
else
    BRSP_t_MDSplus(:,3)=-exl50db(ISHOT,'IPF3',time_range); %read the experimental data of PF current from the MDSplus server
    BRSP_t_MDSplus(:,4)=-exl50db(ISHOT,'IPF4',time_range); %read the experimental data of PF current from the MDSplus server
    BRSP_t_MDSplus(:,5:6)=-repmat(exl50db(ISHOT,'IPF6',time_range),1,2); %read the experimental data of PF current from the MDSplus server
    BRSP_t_MDSplus(:,8)=smooth(exl50db(ISHOT,'ITF',time_range))*1000.0*11/12; %*11/12,1 turn, incomplete ring, open angle=?2*pi/12 rad. read the experimental data of TF current from the MDSplus server
end
BRSP_t_MDSplus(:,7)=-BRSP_t_MDSplus(:,8)/(11/12)*(12*0.08/1.945+2*pi*11/12)/2/pi; %*1.6(or 4.81/2.96) %read the experimental data of TF current from the MDSplus server
BRSP_t =BRSP_t_MDSplus(1:ntime,:);
% BRSP_t((isnan(BRSP_t)==1)==1)=0.0;
clear BRSP_t_MDSplus;
cd('..\Project');

if handles.ntime==1 %only one time slice ntime=1
    PLASMA_t=interp1(t_opt,PLASMA_t,tstart); 
    
    [x_COILS,y_COILS]=meshgrid(tstart,1:COILS_num);
    [x_COILS_raw,y_COILS_raw]=meshgrid(t_opt',1:COILS_num);
    COILS_t = interp2(x_COILS_raw,y_COILS_raw,COILS_t',x_COILS,y_COILS)';
    clear x_COILS_raw y_COILS_raw x_COILS y_COILS;
    
    [x_EXPMP2,y_EXPMP2]=meshgrid(tstart,1:MP_num);
    [x_EXPMP2_raw,y_EXPMP2_raw]=meshgrid(t_opt',1:MP_num);
    EXPMP2_t = interp2(x_EXPMP2_raw,y_EXPMP2_raw,EXPMP2_t',x_EXPMP2,y_EXPMP2)';
    clear x_EXPMP2_raw y_EXPMP2_raw x_EXPMP2 y_EXPMP2;
    
    [x_BRSP,y_BRSP]=meshgrid(tstart,1:BRSP_length);
    [x_BRSP_raw,y_BRSP_raw]=meshgrid(t_opt',1:BRSP_length);
    BRSP_t = interp2(x_BRSP_raw,y_BRSP_raw,BRSP_t',x_BRSP,y_BRSP)';
    clear x_BRSP_raw y_BRSP_raw x_BRSP y_BRSP;
    
    ntime=1;
    t_opt=tstart;
end

%check the diagnostic data of flux loop and remove the drift
checkflux_routine;
cd(handles.pwd_path);

save([path_movie,'\shot',num2str(ISHOT,'%06d'),'data.mat'],'t_opt','PLASMA_t','COILS_t','EXPMP2_t','BRSP_t');%%Error of the Magnetic probe
% load([path_movie,'\shot',num2str(ISHOT,'%06d'),'data.mat']);