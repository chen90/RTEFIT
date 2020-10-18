clear;
clc;

format long;
filename = '.\m000086.00002';
MP_num=63; %number of magnetic probe  unit T
COILS_num=41; %number of flux loops  unit v.s/rad(wb/rad)
PF_num=6; %number of PF coils unit A
Error_Relative_MP = 0.03; %relative error for MP
Error_Relative_COIL = 0.03; %relative error for flux loops
Error_Relative_PF = 0.15; %relative error for PF coils

EXPMP2(1)=textread(filename,'%*s %*s %f %*s',1,'headerlines', 5);
for i=2:MP_num
    EXPMP2(i)=textread(filename,'%f %*s',1,'headerlines', 4+i);
end
EXPMP2_equ=EXPMP2;
% EXPMP2=EXPMP2*1e4; % unit from T to Gs

COILS(1)=textread(filename,'%*s %*s %f %*s',1,'headerlines', 5+MP_num);
for i=2:COILS_num
    COILS(i)=textread(filename,'%f %*s',1,'headerlines', 4+MP_num+i);
end
COILS_equ=COILS;

BRSP(1)=textread(filename,'%*s %*s %f %*s',1,'headerlines', 537);
for i=2:PF_num
    BRSP(i)=textread(filename,'%f %*s',1,'headerlines', 536+i);
end
BRSP_equ=BRSP;

% % %%if the absolute errors for all the magnetic probes are all the same, which implies the relative errors
% % %%will change with the value of the magnetic field
% % 
% % error_MP=max(abs(EXPMP2))*Error_Relative_MP; %absolute error T, here assume the same absolute errors for all the MPs
% % error_COIL=max(abs(COILS))*Error_Relative_COIL; %absolute error wb/rad
% % 
% % EXPMP2=EXPMP2+rand(MP_num,1)'.*2*error_MP-error_MP;
% % COILS=COILS+rand(COILS_num,1)'.*2*error_COIL-error_COIL;

% % %%if the relative errors for all the magnetic probes are all the same, which implies the relative errors
% % %%will not change with the value of the magnetic field

%%
%calculate the error with uniform distribution
% % error_MP=abs(EXPMP2).*Error_Relative_MP; %absolute error T, here calculate the absolute errors for different MPs
% % error_COIL=abs(COILS).*Error_Relative_COIL; %absolute error wb/rad, here calculate the absolute errors for different flux loops
% % error_PF=abs(BRSP).*Error_Relative_PF; %absolute error A, here calculate the absolute errors for different PF coils
% % 
% % EXPMP2=EXPMP2+rand(MP_num,1)'.*2.*error_MP-error_MP;
% % COILS=COILS+rand(COILS_num,1)'.*2.*error_COIL-error_COIL;
% % BRSP=BRSP+rand(PF_num,1)'.*2.*error_PF-error_PF;

%%
%calculate the error with standard normal distribution
EXPMP2=EXPMP2+EXPMP2.*randn(MP_num,1)'.*Error_Relative_MP;
COILS=COILS+COILS.*randn(COILS_num,1)'.*Error_Relative_COIL;
BRSP=BRSP+BRSP.*randn(PF_num,1)'.*Error_Relative_PF;
delta_EXPMP2=sqrt(sum(((EXPMP2-EXPMP2_equ)./EXPMP2_equ).^2));
delta_COILS=sqrt(sum(((COILS-COILS_equ)./COILS_equ).^2));
% delta_BRSP=sqrt(sum(((BRSP-BRSP_equ)./BRSP_equ).^2));
delta_BRSP=sqrt(sum(((BRSP-BRSP_equ)./max(abs(BRSP_equ))).^2));

%calculate the deviation

a=fprintf(',%16.10f',EXPMP2)
b=fprintf(',%16.10f',COILS)
c=fprintf(',%16.10f',BRSP)

% title(['\delta_{MP}=',num2str(Error_Relative_MP),',\delta_{COIL}=',num2str(Error_Relative_COIL),',\delta_{PF}=',num2str(Error_Relative_PF)]);
title(['\delta_{Ip}=0.15,','\delta_{MP}=',num2str(Error_Relative_MP),',\delta_{COIL}=',num2str(Error_Relative_COIL),',\delta_{PF}=',num2str(Error_Relative_PF),10,...
    '\Delta_{Ip}=0.15,','\Delta_{MP}=',num2str(delta_EXPMP2),',\Delta_{COIL}=',num2str(delta_COILS),',\Delta_{PF}=',num2str(delta_BRSP)]);

% st=suptitle(['(a) The configuration size (cm), and (b) the confinement time (s) ',10,...
%     'for FRC geometry($\beta=1$) operating at $\phi\le$',num2str(delta_EXPMP2)]);
% set(st,'Interpreter','latex');
