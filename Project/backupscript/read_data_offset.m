%% note: This version is only useful for the EXL50 discharges with negative PF currents: In this script the currents in the PF coils are set with same signs (negative)
% read data from the MDSplus server, t_raw(ms), PLASMA_t_MDSplus(kA), EXPMP2_t_MDSplus(T), COILS_t_MDSplus(v.s(wb)), BRSP_t_MDSplus(A) 

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

%%%read EXL50 magnetic positions from excel
    %A1 = xlsread('E:\Engineering Simulation\MM_calibration\EXL50_1A_BeforeHeat.xls','MPT');
    %A2 = xlsread('E:\Engineering Simulation\MM_calibration\EXL50_1A_BeforeHeat.xls','MPN');
    A3 = xlsread('E:\RTEFIT\RTEFIT-developer\Project\EXL50_1A_BeforeHeat.xls','FLUX');

%read the experimental data from the MDSplus server
time_range=[num2str(tstart/1000.0),':',num2str(tend/1000.0),':',num2str((tend-tstart)/1000.0/(ntime+1))];
try
    PLASMA_offset=exl50db(ISHOT,'IP','-5:-2:0.03'); %offset
    [PLASMA_t_raw,t_raw]=exl50db(ISHOT,'IP',time_range); %plasma current in unit PLASMA_unit(kA), t_raw in unit s, channels from ip01 to ip04 in different toroidal locations
catch
    errordlg('Error: Can not connect to the EXL50 data base','File Error');
    delete(h);
    error('Error: Can not connect to the EXL50 data base');
end
PLASMA_t=(PLASMA_t_raw(1:ntime)-sum(PLASMA_offset)/length(PLASMA_offset))*1000.0;
t_opt=t_raw(1:ntime); %time seires of the signals, unit: ms
clear PLASMA_offset t_raw;

COILS_t=zeros(length(t_opt),COILS_num);
for i=1:COILS_num
    COIL_offset=exl50db(ISHOT,['FLUX',num2str(i,'%03d')],'-5:-2:0.03'); %offset
    COILS_t_tmp=exl50db(ISHOT,['FLUX',num2str(i,'%03d')],time_range)-sum(COIL_offset)/length(COIL_offset); %read the experimental data of flux loops from the MDSplus server, unit v.s(wb)
    COILS_t(:,i)=COILS_t_tmp(1:ntime);
end
clear COIL_offset COIL_offset;

for i=1:ntime
    COILS_t(i,:)=COILS_t(i,:)./A3(:,6)'.*A3(:,7)'./2./pi;
end

EXPMP2_t=zeros(length(t_opt),MP_num);
% % % % for i=14:MP_T_num
% % % %     EXPMP2_offset=exl50db(ISHOT,['MP',num2str(i,'%03d'),'T'],'-5:-2:0.03'); %offset
% % % %     EXPMP2_t_tmp=exl50db(ISHOT,['MP',num2str(i,'%03d'),'T'],time_range)-sum(EXPMP2_offset)/length(EXPMP2_offset);  %read the experimental data of poloidal magnetic probes from the MDSplus server
% % % %     EXPMP2_t(:,i)=EXPMP2_t_tmp(1:ntime);
% % % % end
% % % % for i=MP_T_num+1:MP_num
% % % %     EXPMP2_offset=exl50db(ISHOT,['MP',num2str(i-MP_T_num,'%03d'),'N'],'-5:-2:0.03'); %offset
% % % %     EXPMP2_t_tmp=exl50db(ISHOT,['MP',num2str(i-MP_T_num,'%03d'),'N'],time_range)-sum(EXPMP2_offset)/length(EXPMP2_offset);  %read the experimental data of Normal direction magnetic probe from the MDSplus server
% % % %     EXPMP2_t(:,i)=EXPMP2_t_tmp(1:ntime);
% % % % end
% % % % clear EXPMP2_t_tmp EXPMP2_offset;

BRSP_length=20;
BRSP_t=ones(length(t_opt),1)*BRSP(1:BRSP_length);
% for i=1:4, the currents in the PF1-4 is in series connection, the PF5-6 is in series connection.
% IPF2 is the current in PF2, IPF1 is the current in PF5-6.
BRSP_offset=exl50db(ISHOT,'IPF2','-5:-2:0.03'); %offset
BRSP_t_tmp=exl50db(ISHOT,'IPF2',time_range)-sum(BRSP_offset)/length(BRSP_offset); %read the experimental data of PF current from the MDSplus server
BRSP_t(:,1:4)=repmat(BRSP_t_tmp(1:ntime),1,4);
% for i=5:6, the PF5-6
% BRSP_t(:,5:6)=0.0;
BRSP_offset=exl50db(ISHOT,'IPF1','-5:-2:0.03'); %offset
BRSP_t_tmp=exl50db(ISHOT,'IPF1',time_range)-sum(BRSP_offset)/length(BRSP_offset); %read the experimental data of PF current from the MDSplus server
BRSP_t(:,5:6)=repmat(BRSP_t_tmp(1:ntime),1,2); 
BRSP_t(:,1:6)=abs(BRSP_t(:,1:6)).*(-1);
clear BRSP_offset BRSP_t_tmp;

save('data.mat','t_opt','PLASMA_t','COILS_t','EXPMP2_t','BRSP_t');%%Error of the Magnetic probe

load('data.mat');
% BRSP_t=-abs(BRSP_t(:,1:6));
BRSP_t(1:4)=BRSP_t(1:4)*4.55/2.3; %A resistor is connected in parallel with PF2 coil (2.3/2.25 R/PF2).