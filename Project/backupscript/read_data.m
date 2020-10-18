% read data from the MDSplus server, t_raw(s), PLASMA_t_MDSplus(A), EXPMP2_t_MDSplus(T), COILS_t_MDSplus(v.s/rad(wb/rad)), BRSP_t_MDSplus(A)     
fs=handles.fs; % sample rate, unit:Hz
MP_num=length(EXPMP2); %number of magnetic probe
MP_T_num=38; %number of poloidal magnetic probe
MP_N_num=25; %number of magnetic probe in Normal direction
COILS_num=length(COILS); %number of flux loops  unit v.s/rad(wb/rad)
PF_num=handles.PF_num; %number of PF coils
PF_turn=22; %turns of PF coils
ISHOT=handles.ISHOT;
ntime = handles.ntime; %number of time slice choosed from the raw data
% setappdata(0,'verbose',1) % automatic plotting
try
    [PLASMA_t_MDSplus,t_raw,PLASMA_unit]=exl50db(ISHOT,'ip02'); %plasma current in unit PLASMA_unit(A), t_raw in unit s, channels from ip01 to ip04 in different toroidal locations
catch
    errordlg('Error: Can not connect to the EXL50 data base','File Error');
    delete(h);
    error('Error: Can not connect to the EXL50 data base');
end
t_signal=t_raw*1000.0; %time seires of the signals, unit: ms

%read the experimental data from the MDSplus server
%interpolate raw data to the selected time slice
t_opt = round(linspace(handles.tstart,handles.tend,ntime)); %select time slice: start from tstart, end at the tend and divided by number of ntime

PLASMA_t = interp1(t_signal,PLASMA_t_MDSplus,t_opt'); 
clear varname_IN1 varvalue_IN1 PLASMA_t_MDSplus;

COILS_t_MDSplus=zeros(length(t_signal),COILS_num);
for i=1:COILS_num
    COILS_t_MDSplus(:,i)=exl50db(ISHOT,['FLUX',num2str(i,'%02d')]); %read the experimental data of flux loops from the MDSplus server, unit v.s/rad(wb/rad)
end
[x_COILS,y_COILS]=meshgrid(t_opt,1:COILS_num);
[x_COILS_raw,y_COILS_raw]=meshgrid(t_signal,1:COILS_num);
COILS_t = interp2(x_COILS_raw,y_COILS_raw,COILS_t_MDSplus',x_COILS,y_COILS)';
clear COILS_t_MDSplus x_COILS_raw y_COILS_raw x_COILS y_COILS;

EXPMP2_t_MDSplus=zeros(length(t_signal),MP_num);
for i=14:MP_T_num
    EXPMP2_t_MDSplus(:,i)=exl50db(ISHOT,['MP',num2str(i,'%03d'),'T']);  %read the experimental data of poloidal magnetic probes from the MDSplus server
end
for i=MP_T_num+1:MP_num
    EXPMP2_t_MDSplus(:,i)=exl50db(ISHOT,['MP',num2str(i-MP_T_num,'%03d'),'N']);  %read the experimental data of Normal direction magnetic probe from the MDSplus server
end
[x_EXPMP2,y_EXPMP2]=meshgrid(t_opt,1:MP_num);
[x_EXPMP2_raw,y_EXPMP2_raw]=meshgrid(t_signal,1:MP_num);
EXPMP2_t = interp2(x_EXPMP2_raw,y_EXPMP2_raw,EXPMP2_t_MDSplus',x_EXPMP2,y_EXPMP2)'; %for the interp2, the value of x coordinate should change for each row, the value of y coordinate should change for each column
clear EXPMP2_t_MDSplus x_EXPMP2_raw y_EXPMP2_raw x_EXPMP2 y_EXPMP2;

BRSP_length=20;
BRSP_t_MDSplus=ones(length(t_signal),1)*BRSP(1:BRSP_length);
for i=1:PF_num
    BRSP_t_MDSplus(:,i)=exl50db(ISHOT,'ipf01')./PF_turn; %read the experimental data of PF current from the MDSplus server
end
[x_BRSP,y_BRSP]=meshgrid(t_opt,1:BRSP_length);
[x_BRSP_raw,y_BRSP_raw]=meshgrid(t_signal,1:BRSP_length);
BRSP_t = interp2(x_BRSP_raw,y_BRSP_raw,BRSP_t_MDSplus',x_BRSP,y_BRSP)';
clear BRSP_t_MDSplus x_BRSP_raw y_BRSP_raw x_BRSP y_BRSP t_signal BRSP_length;