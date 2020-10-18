%% Main program of the real time EFIT
% -inputs:
% m_EXL.ini                namelist, the initial input file for the...
                           %equilibrium reconstruction(Fitting), where based on which the series of m[ISHOT].[1:1:ntime] files are built
% ISHOT[000000]            int,  shot number
% ITIME[00000]             int,  time slice for the discharge of ISHOT (unit:ms)
% PF_num                   number of PF coils for EXL PF_num=6 (default)
% path_gfiles              string,  path for the gfiles on the local machine
% path_mfile               string,  path for the mfiles on the local machine
% path_movie               tring,  path where the movies are saved on the local machine
% RTEFIT_path              tring,  path where the RTEFIT code are located on the Debye cluster
% input_file_path          tring,  path where the input files are located on the Debye cluster
% output_shot_path         tring,  path where the output files are located on the Debye cluster
% FrameRate                int,   number of frames per second
% configurations connecting to the server (default: ENN/Debye)
% Note: Don't change the relative path of the directories and files inside the EFIT code on Debye cluster and in the RTEFIT directory in the local machine
% Note: The time scale for the equilibrium reconstruction is 1ms, the ntime should not be too big, in case of delta(t_opt)<1ms
% Note: The CUTIP must be set CUTIP<PLASMA for the input files

% -outputs:
% m[ISHOT].[1:1:ntime]      namelist, the series of input files for the EFIT equlibrium reconstruction
% g[ISHOT].[ITIME]      TXT, the series of output files for the EFIT equlibrium reconstruction
% a[ISHOT].[ITIME]      TXT, the series of output files for the EFIT equlibrium reconstruction
% 2 movies which are saved in the path of path_movie: one is in format avi, the other is in format GIF

%This script is the interface program between the diagnostic data of
%Magnetic probes, flux loops, Ip, PF coil currents and the EFIT
%reconstruction code. Aimed at simulating the equalibrium reconstruction based
%on the data and show it on the SCREEN, which is similar to the funcion of
%RT-EFIT.

% Edited by Bin Chen in 2019/06/03
% Contact: chenbino@enn.cn (Bin Chen), sunshuyingc@enn.cn (Shuying Sun)
% ENN Group 1989-2019, all rights reserved.

clear; close all; clc;
%% Prepare the inputs and initialisation
% create useful variable according to the diagnostic data
% create a structure varialbe mfile with all input parameters from the m_EXL.ini namelist file
path_mfile='..\RT_input\m_EXL.ini'; %path on local machine
mfile=read_mfile(path_mfile);
listname_mfile=fieldnames(mfile); %fields name of the structure of mfile
varname_IN1=mfile.IN1.varname; %variable names of the IN1 namelist of mfile
varvalue_IN1=mfile.IN1.varvalue; %variable values of the IN1 namelist of mfile
format long;

for i=1:length(varname_IN1)
     eval([varname_IN1{i}, '= varvalue_IN1{i};']); % create varialbes from the IN1 namelist
end

%% read data from the MDSplus

%% create fake data by adding random error to the simulated data from a EFIT equilibrim simulation
% only use, if no data from the MDSplus available
fs=1.0e5; % sample rate, unit:Hz
t_length=5.0; %length of the signal, unit: s
t_signal=0:(t_length*1000.0)/(fs*t_length-1):t_length*1000.0; %time seires of the signals, unit: ms, linspace(0,t_length*1000.0,fs*t_length);
MP_num=length(EXPMP2); %number of magnetic probe  unit T
COILS_num=length(COILS); %number of flux loops  unit v.s/rad(wb/rad)
PF_num=6; %number of PF coils unit A
ISHOT=85;
ntime = 3; %number of time slice choosed from the raw data, (ms)
Error_Relative_MP = 0.03; %relative error for MP
Error_Relative_COIL = 0.03; %relative error for flux loops
Error_Relative_PF = 0.15; %relative error for PF coils
% Error_Relative_MP = 0.03; %relative error for MP
% Error_Relative_COIL = 0.03; %relative error for flux loops
% Error_Relative_PF = 0.15; %relative error for PF coils

%create the raw data: add the error to data with standard normal distribution
EXPMP2_t_fake=ones(length(t_signal),1)*EXPMP2;
EXPMP2_t_fake=EXPMP2_t_fake+EXPMP2_t_fake.*randn(length(t_signal),MP_num).*Error_Relative_MP;
COILS_t_fake=ones(length(t_signal),1)*COILS;
COILS_t_fake=COILS_t_fake+COILS_t_fake.*randn(length(t_signal),COILS_num).*Error_Relative_COIL;
BRSP_t_fake_tmp=ones(length(t_signal),1)*BRSP;
BRSP_t_fake=BRSP_t_fake_tmp+BRSP_t_fake_tmp.*randn(length(t_signal),length(BRSP)).*Error_Relative_PF;
BRSP_t_fake(:,PF_num+1:end)=BRSP_t_fake_tmp(:,PF_num+1:end);
PLASMA_t_fake=ones(length(t_signal),1)*PLASMA;
PLASMA_t_fake=PLASMA_t_fake+PLASMA_t_fake.*randn(length(t_signal),length(PLASMA)).*Error_Relative_PF;

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

%% choose time range to do the simulations and generate simulation inputs for the choosed time slices
figure(1);
clf;
plot(t_signal,abs(PLASMA_t_fake/std(PLASMA_t_fake)),'k');
hold on;
plot(t_signal,abs(EXPMP2_t_fake(:,1)/std(EXPMP2_t_fake(:,1))),'r');
plot(t_signal,abs(COILS_t_fake(:,1)/std(COILS_t_fake(:,1))),'b');
plot(t_signal,abs(BRSP_t_fake(:,3)/std(BRSP_t_fake(:,3))),'g');
legend('location','NorthWest', 'Ip','Probe1','Flux1','PF5');
xlabel('t ( ms )','fontsize',18,'fontname','times');
ylabel('Signals Amplitude (a.u.)','fontsize',18,'fontname','times');
set(gca,'fontweight','bold','fontsize',18,'fontname','times');
set(gca,'linewidth',2);
title('click twice on the plot to choose time range: first click the starting time, second click the ending time')

[tstart,y1]=ginput(1);
[tend,y2]=ginput(1);
close(1);
clear fs listname_mfile y1 y2 Error_Relative_PF Error_Relative_COIL Error_Relative_MP;

if tend<tstart
    error('error: tend < tstart!');
end

if tstart<=t_signal(1)
    tstart=t_signal(1);
end

if tend>=t_signal(end)
    tend=t_signal(end);
end

if abs(tend-tstart)/ntime<=1
    error('error: The ntime is too big, delta(t_opt)<1ms !');
end

%interpolate raw data to the selected time slice
t_opt = round(linspace(tstart,tend,ntime)); %select time slice: start from tstart, end at the tend and divided by number of ntime
[x_EXPMP2,y_EXPMP2]=meshgrid(t_opt,1:MP_num);
[x_EXPMP2_raw,y_EXPMP2_raw]=meshgrid(t_signal,1:MP_num);
EXPMP2_t = interp2(x_EXPMP2_raw,y_EXPMP2_raw,EXPMP2_t_fake',x_EXPMP2,y_EXPMP2)'; %for the interp2, the value of x coordinate should change for each row, the value of y coordinate should change for each column
[x_COILS,y_COILS]=meshgrid(t_opt,1:COILS_num);
[x_COILS_raw,y_COILS_raw]=meshgrid(t_signal,1:COILS_num);
COILS_t = interp2(x_COILS_raw,y_COILS_raw,COILS_t_fake',x_COILS,y_COILS)';
[x_BRSP,y_BRSP]=meshgrid(t_opt,1:length(BRSP));
[x_BRSP_raw,y_BRSP_raw]=meshgrid(t_signal,1:length(BRSP));
BRSP_t = interp2(x_BRSP_raw,y_BRSP_raw,BRSP_t_fake',x_BRSP,y_BRSP)';
PLASMA_t = interp1(t_signal,PLASMA_t_fake',t_opt);

clear varname_IN1 varvalue_IN1 t_signal BRSP_t_fake_tmp BRSP_t_fake EXPMP2_t_fake COILS_t_fake PLASMA_t_fake;

%% read data from the MDSplus, ISHOT_MDSplus, TIME_t_MDSplus(ntime), PLASMA_t_MDSplus(ntime), EXPMP2_t_MDSplus(ntime,:), COILS_t_MDSplus(ntime,:), BRSP_t_MDSplus(ntime,:)
% to do...
% ntime = length(TIME_t_MDSplus); %number of time slice

%% Create the series of m[ISHOT].[1:1:ntime] files for the EFIT equlibrium reconstruction, update the mfile according to the data
% update the mfile according to the fake data
path_inputshot = ['..\RT_input\input_file\shot',num2str(ISHOT,'%06d')]; %path on local machine

if ~exist(path_inputshot,'dir')
   mkdir(path_inputshot);
end

fileinfo = dir(path_inputshot);
if length(fileinfo) ~= 2 %判断是否为空，因为matlab有.和..，所以空文件夹的信息长度为2
    n=length(fileinfo);
    for i=1:n
        if (~strcmp(fileinfo(i).name,'.') && ~strcmp(fileinfo(i).name,'..') )
            %disp(fileinfo(i).name);
            %exist([path_inputshot,'\',fileinfo(i).name])
            delete([path_inputshot,'\',fileinfo(i).name]); %delete all the files in the shot[ISHOT] dir.
        end
    end
end

for j=1:ntime
    for i=1:length(mfile.IN1.varname)
        if j == 1
            if strcmp(mfile.IN1.varname{i},'PLASMA')
                mfile.IN1.varvalue{i}=PLASMA_t(j);
            elseif strcmp(mfile.IN1.varname{i},'EXPMP2')
                mfile.IN1.varvalue{i}=EXPMP2_t(j,:);
            elseif strcmp(mfile.IN1.varname{i},'COILS')
                mfile.IN1.varvalue{i}=COILS_t(j,:);
            elseif strcmp(mfile.IN1.varname{i},'BRSP')
                mfile.IN1.varvalue{i}=BRSP_t(j,:);
            elseif strcmp(mfile.IN1.varname{i},'ITIME')
                 mfile.IN1.varvalue{i}=t_opt(j);
            elseif strcmp(mfile.IN1.varname{i},'ISHOT')
                mfile.IN1.varvalue{i}=ISHOT;
            end
        else
            if strcmp(mfile.IN1.varname{i},'PLASMA')
                mfile.IN1.varvalue{i}=PLASMA_t(j);
            elseif strcmp(mfile.IN1.varname{i},'EXPMP2')
                mfile.IN1.varvalue{i}=EXPMP2_t(j,:);
            elseif strcmp(mfile.IN1.varname{i},'COILS')
                mfile.IN1.varvalue{i}=COILS_t(j,:);
            elseif strcmp(mfile.IN1.varname{i},'BRSP')
                mfile.IN1.varvalue{i}=BRSP_t(j,:);
            elseif strcmp(mfile.IN1.varname{i},'ITIME')
                 mfile.IN1.varvalue{i}=t_opt(j);
            end
        end
    end
    outputstr=[path_inputshot,'\m',num2str(ISHOT),'.',num2str(j)];
    write_mfile(mfile,outputstr);
end

%% connecting to the server (default: ENN/Debye)
ftpobj = ftp('10.1.141.210','public','public','System','UNIX');
% dir(ftpobj);
% RTEFIT_path=cd(ftpobj,'RTEFIT_ENN/RTEFIT'); %path to the RTEFIT directory
RTEFIT_path='/home/public/RTEFIT_ENN/RTEFIT';
% mkdir(ftpobj,'RT_input/input_file');
input_file_path=cd(ftpobj,'RTEFIT_ENN/RTEFIT/RT_input/input_file'); %path to the input_file directory
listing=dir(ftpobj);

for i=1:length(listing)
    if strcmp(listing(i).name,['shot',num2str(ISHOT,'%06d')])
        delete(ftpobj,['shot',num2str(ISHOT,'%06d'),'/*'])
        rmdir(ftpobj,['shot',num2str(ISHOT,'%06d')]);
    end
end

mkdir(ftpobj,['shot',num2str(ISHOT,'%06d')]);
input_shot_path=cd(ftpobj,['shot',num2str(ISHOT,'%06d')]); %path to the input_file directory
mput(ftpobj,[path_inputshot,'\*']); %copy m files to the server
close(ftpobj);

%% SIMPLE CONNECTION: To run the script "runefit" on the Debye cluster
%Copyright (c) 2018, David S. Freedman
addpath('..\ssh2_v2_m1_r7');
command_output = ssh2_simple_command('10.1.141.210','public','public',['cd ',RTEFIT_path,'&& ./runefit ',num2str(ISHOT),' ',num2str(ntime)]);
% this command makes a remote ssh connection to the host HOSTNAME, 
% with username USERNAME, and password PASSWORD. 
% Once connected the command 'ls -la' looked
% for files on the remote host. The output of this command will be
% printed in the matlab command window and returned to the variable command_output.

%% copy g[ISHOT].[ITIME] and a[ISHOT].[ITIME] files from the Debye cluster to the local machine
ftpobj = ftp('10.1.141.210','public','public','System','UNIX');
path_outputshot = ['..\RT_output\Raw\shot',num2str(ISHOT,'%06d')]; %path on local machine

if ~exist(path_outputshot,'dir')
   mkdir(path_outputshot);
else
   fileinfo = dir(path_outputshot);
    if length(fileinfo) ~= 2 %判断是否为空，因为matlab有.和..，所以空文件夹的信息长度为2
        rmdir([path_outputshot,'\*'],'s'); %delete all the directories in the shot[ISHOT] directory on local machine.
    end
end

output_shot_path=cd(ftpobj,['RTEFIT_ENN/RTEFIT/RT_output/Raw/shot',num2str(ISHOT,'%06d')]); %path to the output_file directory on the Debye cluster
mget(ftpobj,'*',path_outputshot);

%% read gfiles from the equilibrium reconstructions and make a movie of the equilibrium evolution
path_movie=['..\RT_output\Plot\shot',num2str(ISHOT,'%06d')];

if ~exist(path_movie,'dir')
   mkdir(path_movie);
end

FrameRate=10; %number of frames per second
movie(path_outputshot,path_movie,ISHOT,t_opt,FrameRate);
