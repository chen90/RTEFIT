function info_equ(path_afiles,Wdata,ax)
%% get basic information of the equilibrium plasmas from EFIT afile
% -inputs:
% path_afiles                      string,  path where the afiles are saved (set it the same as the path of the gfiles in default)
% Wdata                            logical, 0/1, 0 for the parameters list, 1 for the parameters list with data
% ax                               axes object, ax is the axes of the outline

% -outputs:
% it provides various global plasma parameters of interest and displaying on the figures.

% Edited by Bin Chen in 2019/06/05
% Contact: chenbino@enn.cn
% ENN Group 1989-2019, all rights reserved.

%% make 2 movies: one is in format avi, the other is in format GIF
if Wdata
    afile=read_afile(path_afiles);

    xdim=sum(get(ax,'xlim'));
    ylim=get(ax,'ylim');
    ydim=ylim(2);
    N=7; %number of parameters to be shown on the figure
    delta_x=xdim/10;
    delta_y=ydim/N;
    text(xdim+delta_x,ydim-0.5*delta_y,['I_p =',num2str(afile.pasmat/1000,'%6.2f'), ' kA'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-1.5*delta_y,['a =',num2str(afile.aout/100,'%6.2f'),' m'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-2.5*delta_y,['R_{MAGX} =',num2str(afile.rmagx/100,'%6.2f'),' m'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-3.5*delta_y,['R_{OUT} =',num2str(afile.rout/100,'%6.2f'),' m'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-4.5*delta_y,['A =',num2str(afile.rout/afile.aout,'%6.2f')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-5.5*delta_y,['V_{OUT} =',num2str(afile.vout/1e6,'%6.2f'), ' m^3'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-6.5*delta_y,['\beta_p =',num2str(afile.betap,'%6.2f')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-7.5*delta_y,['l_i =',num2str(afile.ali,'%6.2f')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-8.5*delta_y,['\delta_U =',num2str(afile.doutl,'%6.2f')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-9.5*delta_y,['\delta_L =',num2str(afile.doutu,'%6.2f')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-10.5*delta_y,['k =',num2str(afile.eout,'%6.2f')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-11.5*delta_y,['\Psi_{BDRY} =',num2str(afile.sibdry,'%9.2e'),' Wb/rad'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-12.5*delta_y,['\Psi_{MAGX} =',num2str(afile.simagx,'%9.2e'),' Wb/rad'],'horiz','left','color','k','fontsize',16);
    % text(xdim+delta_x,ydim-12*delta_y,['VLOOPT =',num2str(afile.vloopt,'%6.2f'), ' Volt'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-13.5*delta_y,['error =',num2str(afile.terror,'%9.3e')],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-14.5*delta_y,['chi2 =',num2str(afile.tsaisq,'%9.3e')],'horiz','left','color','k','fontsize',16);
else
    xdim=sum(get(ax,'xlim'));
    ylim=get(ax,'ylim');
    ydim=ylim(2);
    N=7; %number of parameters to be shown on the figure
    delta_x=xdim/10;
    delta_y=ydim/N;
    text(xdim+delta_x,ydim-0.5*delta_y,['I_p =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-1.5*delta_y,['a =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-2.5*delta_y,['R_{MAGX} =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-3.5*delta_y,['R_{OUT} =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-4.5*delta_y,['A =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-5.5*delta_y,['V_{OUT} =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-6.5*delta_y,['\beta_p =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-7.5*delta_y,['l_i =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-8.5*delta_y,['\delta_U =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-9.5*delta_y,['\delta_L =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-10.5*delta_y,['k =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-11.5*delta_y,['\Psi_{BDRY} =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-12.5*delta_y,['\Psi_{MAGX} =      '],'horiz','left','color','k','fontsize',16);
    % text(xdim+delta_x,ydim-12*delta_y,['VLOOPT =',num2str(afile.vloopt,'%6.2f'), ' Volt'],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-13.5*delta_y,['error =      '],'horiz','left','color','k','fontsize',16);
    text(xdim+delta_x,ydim-14.5*delta_y,['chi2 =      '],'horiz','left','color','k','fontsize',16);
end
end