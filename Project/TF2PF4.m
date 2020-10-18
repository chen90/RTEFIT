%% note: The currents of PF coils from the database have positive signs. So in the simulations the signs have to be changed according to the sign of the plasma current (always positive in EFIT).
% In this script the currents in the PF coils are set to negative, according to the experiments have been done for EXL-50
% **Check the value of the experiment data before use: COILS_t**
% read data from the MDSplus server, t_raw(s), PLASMA_t_MDSplus(kA), EXPMP2_t_MDSplus(T), COILS_t_MDSplus(v.s(wb)), BRSP_t_MDSplus(A) 

clc;
path(path,strcat(pwd,'\..\DataProc'));  %   add mdsplus path
pmds; %set the path for the mdsplus
sct('exl50',-1) %update the data base
% DP; % open the Data view GUI by Dr. SONG Xianming
cd('..\Project');

COILS_num=39; %number of flux loops  unit v.s/rad(wb/rad)
PF_turn=22; %turns of PF coils
ISHOT=[3573, 3574, 3575];
ntime = 5; %number of time slice choosed from the raw data
tstart = 5000;
tend = 10000; %unit:ms
% sct('exl50',ISHOT); %refresh the shots

%%%read EXL50 magnetic positions from excel
A3 = xlsread('EXL50_1A_asymmetry.xls','FLUX');
        
%read the experimental data from the MDSplus server, the raw data will not be smoothed
time_range=[num2str(tstart/1000.0),':',num2str(tend/1000.0),':',num2str((tend-tstart)/1000.0/ntime)];

for i=1:length(ISHOT)
    [ITF(i), t_signal(i)]=-smooth(exl50db(ISHOT,'ITF',time_range))*1000.0; %read the experimental data of TF current from the MDSplus server
end
    ITF =ITF(1:ntime,:);

%read the flux loop data
for i=1:length(ISHOT)
    COILS_t_MDSplus(i)=exl50dbN(ISHOT,'flux\d*',time_range); %read the experimental data of flux loops from the MDSplus server, unit v.s(wb)
end
COILS_t = COILS_t_MDSplus(1:ntime,1:COILS_num,:);

for i=1:ntime
    for j=1:length(ISHOT)
        COILS_t(i,:,j)=COILS_t(i,:,j)./A3(:,6)'.*A3(:,7)'./2./pi;
    end
end


BRSP_t_MDSplus(:,1)=-exl50db(ISHOT,'IPF4',time_range); %read the experimental data of PF current from the MDSplus server
BRSP_t =BRSP_t_MDSplus(1:ntime,:);
clear BRSP_t_MDSplus;

%check the diagnostic data of flux loop and remove the drift
checkflux_routine;

% save([path_movie,'\shot',num2str(ISHOT,'%06d'),'data.mat'],'t_opt','PLASMA_t','COILS_t','EXPMP2_t','BRSP_t');%%Error of the Magnetic probe
% load([path_movie,'\shot',num2str(ISHOT,'%06d'),'data.mat']);