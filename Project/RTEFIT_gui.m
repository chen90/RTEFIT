function varargout = RTEFIT_gui(varargin) 
%The REAL TIME EFIT　Equilibrium reconstruction code for　EXL50
% RTEFIT_gui MATLAB code for RTEFIT_gui.fig
%      RTEFIT_gui, by itself, creates a new RTEFIT_gui or raises the existing
%      singleton*.
%
%      H = RTEFIT_gui returns the handle to a new RTEFIT_gui or the handle to
%      the existing singleton*.
%
%      RTEFIT_gui('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RTEFIT_gui.M with the given input arguments.
%
%      RTEFIT_gui('Property','Value',...) creates a new RTEFIT_gui or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RTEFIT_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RTEFIT_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RTEFIT_gui

% Last Modified by GUIDE v2.5 12-Dec-2019 16:50:58
% Edited by Bin Chen in 2019/06/16: Programed the base of the GUI code
% Edited by Bin Chen in 2019/11/05: Changed the code with the FLUX Loops data used
% Contact: chenbino@enn.cn (Bin Chen), sunshuyingc@enn.cn (Shuying Sun)
% ENN Group 1989-2019, all rights reserved.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RTEFIT_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @RTEFIT_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before RTEFIT_gui is made visible.
function RTEFIT_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RTEFIT_gui (see VARARGIN)

%% initialization of the basic setting
clc;
handles.pwd_path=pwd;
path(path,strcat(pwd,'\..\DataProc'));  %   add mdsplus path
pmds; %set the path for the mdsplus
sct('exl50',-1) %update the data base
% DP; % open the Data view GUI by Dr. SONG Xianming
cd(handles.pwd_path);
% path(path,strcat('D:\Program Files\MATLAB'));  %add mdsplus path, whoever want to use the mdsplus should add the path of mdsplus to Matlab path
handles.path_mfile='..\RT_input\FreeBequ\m_EXL.ini'; %path on local machine

% read data and basic parameters setting
handles.tstart=0.0;
handles.tend=10000.0; %unit:ms
handles.CUTIP=1000.0; %minimum Ip (unit:Amp)
handles.FrameRate=5;  %number of frames per second
handles.ntime=50.0;  %number of time slice choosed from the raw data
handles.ISHOT=currentshot()-1;
set(handles.shot_number, 'string',handles.ISHOT);
handles.PF_num=6; %number of PF coils
handles.Systems=1; %Operation systems for the local computer:1 or 2，which implies Windows or Unix, respectively
handles.Play = 0; %0 or 1，which implies not display or display the movie after the movie was created, respectively
handles.Smooth = 0; %0 or 1，which implies not smooth or smooth the experimental data, respectively
handles.Method=1; %1 or 2, run the code with equilibrium mode 1 or with fitting mode 2

% read data from Server-----------------------------
handles.Data_source=1; %Data_source:1 or 2，which implies read data from experiment or create fake data, respectively
handles.server='192.168.20.22';
handles.iserver=1;

% create fake data-----------------------------
handles.fs=1.0e5; % sample rate, unit:Hz
handles.t_length=10.0; %length of the signal, unit: s
handles.Error_Relative=0.03; %relative error for Ip, MP, flux loops and PF coils

% Choose default command line output for RTEFIT_gui
handles.output = hObject;

% remove the drift of the flux loop data
handles.Rdrift_on=0; % 0/1, if 0 the raw data is used, if 1 remove the drift of the flux loop data
handles.t_start_Rdrift = 2; %unit: s, start time to do the linear fiting of the data to remove the drift from the raw data
handles.t_end_Rdrift = 5; %unit: s, end time to do the linear fiting of the data to remove the drift from the raw data

% plot the outline in the figure of figure_GUI
axes(handles.img);
ax{1}=gca;
set(ax{1},'units','normalized','position',[0.5685,0.1383,0.2128,0.6816]);
outline; % plot the outline for the device
axes('units','normalized','position',[0.5685,0.1383,0.2128,0.6816]);
outline; % plot the outline for the device
ax{2}=gca;
info_equ('',0,ax{1});
handles.ax=ax; %the two axeses are saved into the GUI handles for the usage of movie.m

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes RTEFIT_gui wait for user response (see UIRESUME)
% uiwait(handles.figure_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = RTEFIT_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in server_popupmenu.
function server_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to server_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
switch val
    case 1
        handles.iserver=1;
        handles.server='192.168.20.151'; %server 1
    case 2  
        handles.iserver=2;
        handles.server='192.168.20.22'; %server 2
end
guidata(hObject, handles); % Update handles structure
% Hints: contents = get(hObject,'String') returns server_popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from server_popupmenu


% --- Executes during object creation, after setting all properties.
function server_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to server_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Datasource.
function Datasource_Callback(hObject, eventdata, handles)
% hObject    handle to Datasource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
switch val
    case 1
        handles.Data_source=1; %Read data from experiment
    case 2  
        handles.Data_source=2; %Create fake data
end
guidata(hObject, handles); % Update handles structure
% Hints: contents = cellstr(get(hObject,'String')) returns Datasource contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Datasource


% --- Executes during object creation, after setting all properties.
function Datasource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Datasource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% ----------------------------------------------------------------------------------
%  **********************************Main Program********************************************
%------------------------------------------------------------------------------------------------
%
% --- Executes on button press in Movie_pushbutton.
function Movie_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Movie_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
    % Edited by Bin Chen in 2019/11/05: Changed the code with the FLUX Loops data used
    % Contact: chenbino@enn.cn (Bin Chen), sunshuyingc@enn.cn (Shuying Sun)
    % ENN Group 1989-2019, all rights reserved.
    h = msgbox('Start running the main program of the real time EFIT');
    %% Prepare the inputs and initialisation
    % create useful variable according to the diagnostic data
    % create a structure varialbe mfile with all input parameters from the m_EXL.ini namelist file
    path_mfile=handles.path_mfile; %path on local machine
    mfile=read_mfile(path_mfile);
    listname_mfile=fieldnames(mfile); %fields name of the structure of mfile
    varname_IN1=mfile.IN1.varname; %variable names of the IN1 namelist of mfile
    varvalue_IN1=mfile.IN1.varvalue; %variable values of the IN1 namelist of mfile
    format long;

    for i=1:length(varname_IN1)
         eval([varname_IN1{i}, '= varvalue_IN1{i};']); % create varialbes from the IN1 namelist
    end
    ISHOT=handles.ISHOT;
    delete(h);

    if handles.Systems==1 %Windows operation system
        path_movie=['..\RT_output\Plot\shot',num2str(ISHOT,'%06d')];
    elseif handles.Systems==2 %Unix operation system
        path_movie=['../RT_output/Plot/shot',num2str(ISHOT,'%06d')];
    end

    if ~exist(path_movie,'dir')
       mkdir(path_movie);
    end

    %% read data from the MDSplus or create fake data for the reconstruction
    if handles.Data_source==1  % read data from the MDSplus
        % read data from the MDSplus, ISHOT_MDSplus, TIME_t_MDSplus(ntime), PLASMA_t_MDSplus(ntime), EXPMP2_t_MDSplus(ntime,:), COILS_t_MDSplus(ntime,:), BRSP_t_MDSplus(ntime,:)     
        h = msgbox('Read data from the MDSplus');
        %errordlg('Error: to do...','File Error');
        %return;
        read_data;
        delete(h);
    elseif handles.Data_source==2 %create fake data
        % create fake data by adding random error to the simulated data from a EFIT equilibrim simulation
        % only use, if no data from the MDSplus available
        h = msgbox('Create fake data');
        fake_data;
        delete(h);
    end

    %% Create the series of m[ISHOT].[1:1:ntime] files for the EFIT equlibrium reconstruction, update the mfile according to the data
    % update the mfile according to the fake data
    h = msgbox('Create input for the EFIT equlibrium reconstruction');
    if handles.Systems==1 %Windows operation system
        path_inputshot = ['..\RT_input\input_file\shot',num2str(ISHOT,'%06d')]; %path on local machine
    elseif handles.Systems==2 %Unix operation system
        path_inputshot = ['../RT_input/input_file/shot',num2str(ISHOT,'%06d')]; %path on local machine
    end     

    if ~exist(path_inputshot,'dir')
       mkdir(path_inputshot);
    end

    fileinfo = dir(path_inputshot);
    if length(fileinfo) ~= 2 %判断是否为空，因为matlab有.和..，所以空文件夹的信息长度为2
        n=length(fileinfo);
        for i=1:n
            if (~strcmp(fileinfo(i).name,'.') && ~strcmp(fileinfo(i).name,'..') )
                %disp(fileinfo(i).name);
                    if handles.Systems==1 %Windows operation system
                        %exist([path_inputshot,'\',fileinfo(i).name])
                        delete([path_inputshot,'\',fileinfo(i).name]); %delete all the files in the shot[ISHOT] dir.
                    elseif handles.Systems==2 %Unix operation system
                        %exist([path_inputshot,'/',fileinfo(i).name])
                        delete([path_inputshot,'/',fileinfo(i).name]); %delete all the files in the shot[ISHOT] dir.
                    end     
                
            end
        end
    end

    for j=1:ntime
        for i=1:length(mfile.IN1.varname)
            if j == 1
                if strcmpi(mfile.IN1.varname{i},'PLASMA')
                    mfile.IN1.varvalue{i}=PLASMA_t(j);
                elseif strcmpi(mfile.IN1.varname{i},'EXPMP2')
                    mfile.IN1.varvalue{i}=EXPMP2_t(j,:);
                elseif strcmpi(mfile.IN1.varname{i},'COILS')
                    mfile.IN1.varvalue{i}=COILS_t(j,:);
                elseif strcmpi(mfile.IN1.varname{i},'BRSP')
                    mfile.IN1.varvalue{i}=BRSP_t(j,:);
                elseif strcmpi(mfile.IN1.varname{i},'ITIME')
                     mfile.IN1.varvalue{i}=t_opt(j);
                elseif strcmpi(mfile.IN1.varname{i},'ISHOT')
                    mfile.IN1.varvalue{i}=ISHOT;
                elseif strcmpi(mfile.IN1.varname{i},'CUTIP')
                    mfile.IN1.varvalue{i}=handles.CUTIP;                    
                end
            else
                if strcmpi(mfile.IN1.varname{i},'PLASMA')
                    mfile.IN1.varvalue{i}=PLASMA_t(j);
                elseif strcmpi(mfile.IN1.varname{i},'EXPMP2')
                    mfile.IN1.varvalue{i}=EXPMP2_t(j,:);
                elseif strcmpi(mfile.IN1.varname{i},'COILS')
                    mfile.IN1.varvalue{i}=COILS_t(j,:);
                elseif strcmpi(mfile.IN1.varname{i},'BRSP')
                    mfile.IN1.varvalue{i}=BRSP_t(j,:);
                elseif strcmpi(mfile.IN1.varname{i},'ITIME')
                     mfile.IN1.varvalue{i}=t_opt(j);
                end
            end
        end
            if handles.Systems==1 %Windows operation system
                outputstr=[path_inputshot,'\m',num2str(ISHOT),'.',num2str(j)];
            elseif handles.Systems==2 %Unix operation system
                outputstr=[path_inputshot,'/m',num2str(ISHOT),'.',num2str(j)];
            end     
        write_mfile(mfile,outputstr);
    end
    delete(h);
    
    %% connecting to the server (default: ENN/Debye)
    h = msgbox('Connecting to the server');
    ftpobj = ftp('192.168.20.22','public','public','System','UNIX');
    % dir(ftpobj);
    % RTEFIT_path=cd(ftpobj,'RTEFIT_ENN_pf78/RTEFIT'); %path to the RTEFIT directory
    RTEFIT_path='/home/public/RTEFIT_ENN_pf78/RTEFIT';
    % mkdir(ftpobj,'RT_input/input_file');
    input_file_path=cd(ftpobj,'RTEFIT_ENN_pf78/RTEFIT/RT_input/input_file'); %path to the input_file directory
    listing=dir(ftpobj);

    for i=1:length(listing)
        if strcmp(listing(i).name,['shot',num2str(ISHOT,'%06d')])
            delete(ftpobj,['shot',num2str(ISHOT,'%06d'),'/*']);
            rmdir(ftpobj,['shot',num2str(ISHOT,'%06d')]);
        end
    end

    mkdir(ftpobj,['shot',num2str(ISHOT,'%06d')]);
    input_shot_path=cd(ftpobj,['shot',num2str(ISHOT,'%06d')]); %path to the input_file directory
    
    if handles.Systems==1 %Windows operation system
        mput(ftpobj,[path_inputshot,'\*']); %copy m files to the server
    elseif handles.Systems==2 %Unix operation system
        mput(ftpobj,[path_inputshot,'/*']); %copy m files to the server
    end

    close(ftpobj);
    clear listing mfile;
    delete(h);
    
    %% SIMPLE CONNECTION: To run the script "runefit" on the Debye cluster
    %Copyright (c) 2018, David S. Freedman
    h = msgbox('Run EFIT on the Debye cluster');
    
    if handles.Systems==1 %Windows operation system
         addpath('..\ssh2_v2_m1_r7');
    elseif handles.Systems==2 %Unix operation system
         addpath('../ssh2_v2_m1_r7');
    end
   
    Run_output = ssh2_simple_command('192.168.20.22','public','public',['cd ',RTEFIT_path,'&& ./runefit ',num2str(ISHOT),' ',num2str(handles.ntime)]);
    delete(h);
    % this command makes a remote ssh connection to the host HOSTNAME, 
    % with username USERNAME, and password PASSWORD. 
    % Once connected the command 'ls -la' looked
    % for files on the remote host. The output of this command will be
    % printed in the matlab command window and returned to the variable command_output.

    %% copy g[ISHOT].[ITIME] and a[ISHOT].[ITIME] files from the Debye cluster to the local machine
    h = msgbox('Copy data from the Debye cluster to the local machine');
    ftpobj = ftp('192.168.20.22','public','public','System','UNIX');
    
    if handles.Systems==1 %Windows operation system
         path_outputshot = ['..\RT_output\Raw\shot',num2str(ISHOT,'%06d')]; %path on local machine
    elseif handles.Systems==2 %Unix operation system
         path_outputshot = ['../RT_output/Raw/shot',num2str(ISHOT,'%06d')]; %path on local machine
    end

    if ~exist(path_outputshot,'dir')
       mkdir(path_outputshot);
    else
       fileinfo = dir(path_outputshot);
        if length(fileinfo) ~= 2 %判断是否为空，因为matlab有.和..，所以空文件夹的信息长度为2
                if handles.Systems==1 %Windows operation system
                     rmdir([path_outputshot,'\*'],'s'); %delete all the directories in the shot[ISHOT] directory on local machine.
                elseif handles.Systems==2 %Unix operation system
                     rmdir([path_outputshot,'/*'],'s'); %delete all the directories in the shot[ISHOT] directory on local machine.
                end
        end
    end

    output_shot_path=cd(ftpobj,['RTEFIT_ENN_pf78/RTEFIT/RT_output/Raw/shot',num2str(ISHOT,'%06d')]); %path to the output_file directory on the Debye cluster
    mget(ftpobj,'*',path_outputshot);
    delete(h);

    %% read gfiles from the equilibrium reconstructions and make a movie of the equilibrium evolution
    h = msgbox('Make a movie of the equilibrium evolution'); 
% % path_outputshot = ['..\RT_output\Raw\shot',num2str(ISHOT,'%06d')]; %path on local machine
% % t_opt=0:102:5000;
    FrameRate=handles.FrameRate; %number of frames per second
    movie(path_outputshot,path_movie,ISHOT,t_opt,FrameRate,handles.Systems,handles.Play,handles.ax);
    delete(h);
%
%------------------------------------------------------------------------------------------------
%************************************End of Main Program*************************************
%------------------------------------------------------------------------------------------------



% --- Executes on button press in time_selection.
function time_selection_Callback(hObject, eventdata, handles)
% hObject    handle to time_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % create useful variable according to the diagnostic data
    % create a structure varialbe mfile with all input parameters from the m_EXL.ini namelist file
    path_mfile=handles.path_mfile; %path on local machine
    mfile=read_mfile(path_mfile);
    listname_mfile=fieldnames(mfile); %fields name of the structure of mfile
    varname_IN1=mfile.IN1.varname; %variable names of the IN1 namelist of mfile
    varvalue_IN1=mfile.IN1.varvalue; %variable values of the IN1 namelist of mfile
    format long;

    for i=1:length(varname_IN1)
         eval([varname_IN1{i}, '= varvalue_IN1{i};']); % create varialbes from the IN1 namelist
    end
    
    % choose time range to do the simulations
    if handles.Data_source==2 %fake data
        h = msgbox('Select time interval from fake data');
        fs=handles.fs/100.0; % sample rate, unit:100 Hz (Only for demonstration, divede by 100 to short the running time)
        t_length=handles.t_length; %length of the signal, unit: s
        t_signal=0:(t_length*1000.0)/(fs*t_length-1):t_length*1000.0; %time seires of the signals, unit: ms, linspace(0,t_length*1000.0,fs*t_length);
        MP_num=length(EXPMP2); %number of magnetic probe  unit T
        COILS_num=length(COILS); %number of flux loops  unit v.s/rad(wb/rad)
        PF_num=handles.PF_num; %number of PF coils
        ntime = handles.ntime; %number of time slice choosed from the raw data
        Error_Relative=handles.Error_Relative;
        Error_Relative_MP = Error_Relative; %relative error for MP
        Error_Relative_COIL = Error_Relative; %relative error for flux loops
        Error_Relative_PF = Error_Relative; %relative error for PF coils

        %create the raw data: add the error to data with standard normal distribution
        PLASMA_t_fake=ones(length(t_signal),1)*PLASMA;
        PLASMA_t_fake=PLASMA_t_fake+PLASMA_t_fake.*randn(length(t_signal),length(PLASMA)).*Error_Relative_PF;
        EXPMP2_t_fake=ones(length(t_signal),1)*EXPMP2;
        EXPMP2_t_fake=EXPMP2_t_fake+EXPMP2_t_fake.*randn(length(t_signal),MP_num).*Error_Relative_MP;
        COILS_t_fake=ones(length(t_signal),1)*COILS;
        COILS_t_fake=COILS_t_fake+COILS_t_fake.*randn(length(t_signal),COILS_num).*Error_Relative_COIL;
        BRSP_length=20;
        BRSP_t_fake_tmp=ones(length(t_signal),1)*BRSP(1:BRSP_length);
        BRSP_t_fake=BRSP_t_fake_tmp+BRSP_t_fake_tmp.*randn(length(t_signal),BRSP_length).*Error_Relative_PF;
        BRSP_t_fake(:,PF_num+1:end)=BRSP_t_fake_tmp(:,PF_num+1:end);

        figure(1);
        clf;
        set(gcf,'units','normalized','position',[0.1 0.3 0.8 0.4]);
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
        title('click twice on the plot to choose time range: first click is the starting time, second click is the ending time')

        [tstart,y1]=ginput(1);
        [tend,y2]=ginput(1);
        close(1);

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

        handles.tstart=tstart;
        set(handles.tstart_edit, 'string', roundn(tstart/1000.0,-3));
        handles.tend=tend;
        set(handles.tend_edit, 'string', roundn(tend/1000.0,-3));
        delete(h);
        
    elseif handles.Data_source==1 % select time interval from experiment 
        % select time interval from experiment 
        h = msgbox('Select time interval from experiment');
        %errordlg('Error: to do...','File Error');
        %return;
        
        ntime = handles.ntime; %number of time slice choosed from the raw data
        ISHOT=handles.ISHOT;

        try
            [ip,t_signal] = exl50db(ISHOT,'IP','0:10:0.05'); %plasma current in unit PLASMA_unit(kA), t_signal in unit s, channels from ip01 to ip04 in different toroidal locations
        catch
            errordlg('Error: Can not connect to the EXL50 data base','File Error');
            delete(h);
            error('Error: Can not connect to the EXL50 data base');
        end

        if ~isempty(find(isnan(ip)==1,1)) || ~isempty(find(isempty(ip)==1,1))
            error('Error: The plasma current is NaN or empty');
        end

        [ipf1,t2] = exl50db(ISHOT,'IPF1','0:10:0.05');
        [ipf6,t3] = exl50db(ISHOT,'IPF6','0:10:0.05');
        [ipf3,t4] = exl50db(ISHOT,'IPF3','0:10:0.05');
        [ipf4,t5] = exl50db(ISHOT,'IPF4','0:10:0.05');
        [mwi,t6] = exl50db(ISHOT,'MWI_NE001','0:10:0.05');
%         [blockgas,t7] = exl50db(ISHOT,'BLOCKGAS','0:10:0.05');
        [m1pin,t8] = exl50db(ISHOT,'M1_Pin','0:10:0.05');
        [m1pout,t9] = exl50db(ISHOT,'M1_Pout','0:10:0.05');
        [m2pin,t10] = exl50db(ISHOT,'M2_Pin','0:10:0.05');
        [m2pout,t11] = exl50db(ISHOT,'M2_Pout','0:10:0.05');
        [m3pin,t12] = exl50db(ISHOT,'M3_Pin','0:10:0.05');
        [m3pout,t13] = exl50db(ISHOT,'M3_Pout','0:10:0.05');
        cd(handles.pwd_path);
        
        fig=figure(1);
        set(fig,'units','normalized','position',[0.1 0.05 0.8 0.85]);
        ax1=axes('position',[0.1,0.77,0.8,0.18]);
        plot(t_signal,ip,'k','linewidth',3);
        grid on;axis tight;
        xlim([0,10]);
        % leg=legend('IP(kA)','location','northwest');
        % set(leg,'color','none','box','off');
        title(['Shot #',num2str(ISHOT)]);
        ylabel('Ip (kA)','fontsize',18,'fontname','times');
        % text(11,2.8,['shot #',num2str(shot)]);
        set(gca,'fontweight','bold','fontsize',18,'fontname','times','xTicklabel',{});

        axes('position',[0.1,0.59,0.8,0.18]);
        plot(t2,ipf1,'k',t3,ipf6,'r',t4,ipf3,'b',t5,ipf4,'m','linewidth',3);
        grid on;axis tight;
        xlim([0,10]);
        leg=legend('IPF1','IPF6','IPF3','IPF4','location','northwest','orientation','vertical');
        set(leg,'color','none','box','off');
        ylabel('Current (A)');
        set(gca,'fontweight','bold','fontsize',18,'fontname','times','xTicklabel',{});

        axes('position',[0.1,0.41,0.8,0.18]);
        plot(t6,mwi,'b','linewidth',3);
        grid on;axis tight;
        xlim([0,10]);
        leg=legend('MWI','location','northwest','orientation','horizontal');
        set(leg,'color','none','box','off');
        ylabel('n_e\times 10^1^7m^-^2');
        set(gca,'fontweight','bold','fontsize',18,'fontname','times','xTicklabel',{});

        %axes('position',[0.1,0.23,0.8,0.18]);
        %plot(t7,blockgas,'m','linewidth',3);
        %grid on;axis tight;
        %xlim([-2,12]);
        %leg=legend('Blockgas','location','north','orientation','horizontal');
        %set(leg,'color','none','box','off');
        % % ylabel('Blockgas');
        %set(gca,'fontweight','bold','fontsize',18,'fontname','times','xTicklabel',{});

        %axes('position',[0.1,0.08,0.8,0.15]);
        axes('position',[0.1,0.23,0.8,0.18]);
        plot(t8,m1pin,'k',t9,m1pout,'k-.',t10,m2pin,'r',t11,m2pout,'r-.',t12,m3pin,'b',t13,m3pout,'b-.','linewidth',3);
        grid on;axis tight;
        xlim([0,10]);
        %axis([0 12 0 20]);
        leg=legend('M1_-PIN','M1_-POUT','M2_-PIN','M2_-POUT','M3_-PIN','M3_-POUT','location','north','orientation','horizontal');
        set(leg,'color','none','box','off');
        ylabel('Power (kW)');
        % xlabel('Time(s)','position',[6,-3]);
        xlabel('Time(s)');
        set(gca,'fontweight','bold','fontsize',18,'fontname','times');
        
        if handles.Systems==1 %Windows operation system
           path_movie=['..\RT_output\Plot\shot',num2str(ISHOT,'%06d')];
            if ~exist(path_movie,'dir')
               mkdir(path_movie);
            end
           saveas(fig,[path_movie,'\shot',num2str(ISHOT,'%06d'),'.png']);
        elseif handles.Systems==2 %Unix operation system
           path_movie=['../RT_output/Plot/shot',num2str(ISHOT,'%06d')];
            if ~exist(path_movie,'dir')
               mkdir(path_movie);
            end
           saveas(fig,[path_movie,'/shot',num2str(ISHOT,'%06d'),'.png']);
        end
        
        axes(ax1);
        title(['Shot #',num2str(ISHOT), ...
            ': click twice on the plot to choose time range']);
        
        [tstart,y1]=ginput(1);
        [tend,y2]=ginput(1);
        close(1);
        t_signal=t_signal*1000.0;
        tstart=tstart*1000.0; %units from s to ms
        tend=tend*1000.0;

        if tend<tstart
            error('error: tend < tstart when selecting the time range for simulation!');
        end
        
        if tend<t_signal(1)
            warning('Warning: tend < minimum time ! Set tend = minimum time');
            tend=t_signal(1);
        end
        
        if tstart>t_signal(end)
            warning('Warning: tstart > maximum time ! Set tend = maximum time');
            tstart=t_signal(end);
        end

        if tstart<=t_signal(1)
            tstart=t_signal(1);
        end

        if tend>=t_signal(end)
            tend=t_signal(end);
        end

        if abs(tend-tstart)/(ntime+1)<=0.01 %the sample rate is 100 kHz for flux loop data, so the delta(t_opt) must larger than 10 \mu s
            warning(['Warning: The ntime is too big, delta(t_opt)<10' char(956) 's ! Set Sample=1.0']);
            handles.ntime=1.0;
            set(handles.ntime_edit, 'string', 1.0);
        end

        handles.tstart=tstart;
        set(handles.tstart_edit, 'string', roundn(tstart/1000.0,-3));
        handles.tend=tend;
        set(handles.tend_edit, 'string', roundn(tend/1000.0,-3));
        delete(h);
    end
guidata(hObject, handles);


function fs_Callback(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.fs=eval(a)*1.0e+3;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of fs as text
%        str2double(get(hObject,'String')) returns contents of fs as a double


% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function t_length_Callback(hObject, eventdata, handles)
% hObject    handle to t_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.t_length=eval(a);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of t_length as text
%        str2double(get(hObject,'String')) returns contents of t_length as a double


% --- Executes during object creation, after setting all properties.
function t_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function error_Callback(hObject, eventdata, handles)
% hObject    handle to error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.Error_Relative=eval(a)*0.01;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of error as text
%        str2double(get(hObject,'String')) returns contents of error as a double


% --- Executes during object creation, after setting all properties.
function error_CreateFcn(hObject, eventdata, handles)
% hObject    handle to error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function shot_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shot_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function shot_number_Callback(hObject, eventdata, handles)
% hObject    handle to shot_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.ISHOT=eval(a);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of shot_number as text
%        str2double(get(hObject,'String')) returns contents of shot_number as a double


% --- Executes during object creation, after setting all properties.
function tstart_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tstart_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function tstart_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tstart_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.tstart=eval(a)*1000.0;

if handles.tstart==handles.tend || abs(handles.tend-handles.tstart)/(handles.ntime+1)<=0.01 %the sample rate is 100 kHz for flux loop data, so the delta(t_opt) must larger than 10 \mu s
    warning(['Warning: The ntime is too big, delta(t_opt)<10' char(956) 's ! Set Sample=1.0']);
    handles.ntime=1.0;
    set(handles.ntime_edit, 'string', 1.0);
end 

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of tstart_edit as text
%        str2double(get(hObject,'String')) returns contents of tstart_edit as a double


% --- Executes during object creation, after setting all properties.
function tend_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tend_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function tend_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tend_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.tend=eval(a)*1000.0;

if handles.tstart==handles.tend || abs(handles.tend-handles.tstart)/(handles.ntime+1)<=0.01 %the sample rate is 100 kHz for flux loop data, so the delta(t_opt) must larger than 10 \mu s
    warning(['Warning: The ntime is too big, delta(t_opt)<10' char(956) 's ! Set Sample=1.0']);
    handles.ntime=1.0;
    set(handles.ntime_edit, 'string', 1.0);
end 

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of tend_edit as text
%        str2double(get(hObject,'String')) returns contents of tend_edit
%        as a double


% --- Executes during object creation, after setting all properties.
function minip_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minip_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function minip_edit_Callback(hObject, eventdata, handles)
% hObject    handle to minip_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.CUTIP=eval(a)*1.0e+3;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of minip_edit as text
%        str2double(get(hObject,'String')) returns contents of minip_edit as a double


% --- Executes during object creation, after setting all properties.
function framerate_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function framerate_edit_Callback(hObject, eventdata, handles)
% hObject    handle to framerate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.FrameRate=eval(a);
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of framerate_edit as text
%        str2double(get(hObject,'String')) returns contents of framerate_edit as a double


% --- Executes during object creation, after setting all properties.
function ntime_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ntime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function ntime_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ntime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=get(hObject, 'String');
handles.ntime=eval(a);

if handles.tstart==handles.tend;
    handles.ntime=1.0;
    set(handles.ntime_edit, 'string', 1.0);
end

if handles.ntime==1.0;
    handles.tend=handles.tstart;
    set(handles.tend_edit, 'string', roundn(handles.tend/1000.0,-3));
end 
    
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of ntime_edit as text
%        str2double(get(hObject,'String')) returns contents of ntime_edit as a double


% --- Executes during object creation, after setting all properties.
function figure_GUI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in System.
function System_Callback(hObject, eventdata, handles)
% hObject    handle to System (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
switch val
    case 1
        handles.Systems=1; %The local machine is Windows operation system
    case 2  
        handles.Systems=2; %The local machine is Unix based operation system
        if handles.Method==1
            handles.path_mfile='../RT_input/FreeBequ/m_EXL.ini'; %path on local machine
        elseif handles.Method==2
            handles.path_mfile='../RT_input/Fitting/m_EXL.ini'; %path on local machine
        end
end
guidata(hObject, handles); % Update handles structure
% Hints: contents = cellstr(get(hObject,'String')) returns System contents as cell array
%        contents{get(hObject,'Value')} returns selected item from System


% --- Executes during object creation, after setting all properties.
function System_CreateFcn(hObject, eventdata, handles)
% hObject    handle to System (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Play = get(hObject,'Value'); %0 or 1，which implies not display or display the movie after the movie was created, respectively
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes during object creation, after setting all properties.
function img_CreateFcn(hObject, eventdata, handles)
% creating an axes named img for the movie in the GUI figure_GUI
% hObject    handle to img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% set(hObject,'xticklabel',[]);
% set(hObject,'yticklabel',[]);
%guidata(hObject, handles);
% set(hObject,'box','off');
% Hint: place code in OpeningFcn to populate img


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Smooth = get(hObject,'Value'); %0 or 1，which implies not smooth or smooth the experimental data, respectively
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in Rdrift.
function Rdrift_Callback(hObject, eventdata, handles)
% hObject    handle to Rdrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% select time interval of the current flat top  
% h = msgbox('Remove the drift from the raw data');
handles.Rdrift_on=1; % 0/1, if 0 the raw data is used, if 1 remove the drift of the flux loop data
set(handles.Rdrift_select, 'Value', 1);
% ISHOT=handles.ISHOT;
% 
% try
%     [ip,t_signal] = exl50db(ISHOT,'IP','0:10:0.05'); %plasma current in unit PLASMA_unit(kA), t_signal in unit s, channels from ip01 to ip04 in different toroidal locations
% catch
%     errordlg('Error: Can not connect to the EXL50 data base','File Error');
%     delete(h);
%     error('Error: Can not connect to the EXL50 data base');
% end
% 
% if ~isempty(find(isnan(ip)==1,1)) || ~isempty(find(isempty(ip)==1,1))
%     error('Error: The plasma current is NaN or empty');
% end
% 
% [FLUX003,t2] = exl50db(ISHOT,'FLUX003','0:10:0.05');
% [FLUX007,t3] = exl50db(ISHOT,'FLUX007','0:10:0.05');
% [FLUX017,t4] = exl50db(ISHOT,'FLUX017','0:10:0.05');
% [FLUX026,t5] = exl50db(ISHOT,'FLUX026','0:10:0.05');
% 
% fig=figure(1);
% set(fig,'units','normalized','position',[0.1 0.05 0.8 0.8]);
% axes('position',[0.1,0.77,0.8,0.18]);
% plot(t_signal,ip,'linewidth',2);
% grid on;xlim([-5,13]);
% leg=legend('IP(kA)','location','north');
% set(leg,'color','none','box','off');
% axis tight;
% title(['Shot #',num2str(ISHOT), ...
%     ': click twice on the plot to choose the flat top of Ip']);
% ylabel('Ip (kA)','fontsize',18,'fontname','times');
% % text(11,2.8,['shot #',num2str(shot)]);
% set(gca,'fontweight','bold','fontsize',18,'fontname','times');
% 
% axes('position',[0.1,0.59,0.8,0.18]);
% plot(t2,FLUX003,t3,FLUX007,t4,FLUX017,t5,FLUX026,'linewidth',2);
% grid on;axis tight; %xlim([-5,13]);
% leg=legend('FLUX003','FLUX007','FLUX017','FLUX026','location','south','orientation','horizontal');
% set(leg,'color','none','box','off');
% ylabel('a.u.');
% set(gca,'fontweight','bold','fontsize',18,'fontname','times');
% 
% [tstart,y1]=ginput(1);
% [tend,y2]=ginput(1);
% close(1);
% 
% if tend<tstart
%     error('error: tend < tstart!');
% end
% 
% if tstart<=t_signal(1) || tend>=t_signal(end)
%     error('error: out of the time range');
% end
% 
% handles.Rdrift_on=1; % 0/1, if 0 the raw data is used, if 1 remove the drift of the flux loop data
% handles.t_start_Rdrift = tstart; %unit: s, start time to do the linear fiting of the data to remove the drift from the raw data
% handles.t_end_Rdrift = tend; %unit: s, end time to do the linear fiting of the data to remove the drift from the raw data
% delete(h);
guidata(hObject, handles);


% --- Executes on selection change in Method.
function Method_Callback(hObject, eventdata, handles)
% hObject    handle to Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
switch val
    case 1
        handles.Method=1; %run the code with equilibrium mode
        if handles.Systems==1
            handles.path_mfile='..\RT_input\FreeBequ\m_EXL.ini'; %path on local machine
        elseif handles.Systems==2
            handles.path_mfile='../RT_input/FreeBequ/m_EXL.ini'; %path on local machine
        end
    case 2  
        handles.Method=2; %run the code with fitting mode
        if handles.Systems==1
            handles.path_mfile='..\RT_input\Fitting\m_EXL.ini'; %path on local machine
        elseif handles.Systems==2
            handles.path_mfile='../RT_input/Fitting/m_EXL.ini'; %path on local machine
        end
end
guidata(hObject, handles); % Update handles structure
% Hints: contents = cellstr(get(hObject,'String')) returns Method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Method


% --- Executes during object creation, after setting all properties.
function Method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text100_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text101_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text101 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in Rdrift_select.
function Rdrift_select_Callback(hObject, eventdata, handles)
% hObject    handle to Rdrift_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Rdrift_on = get(hObject,'Value'); % % 0/1, if 0 the raw data is used, if 1 remove the drift of the flux loop data
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Rdrift_select
