function gfile=fluxsurface(path_gfile,varargin)
%% function: plot the flux surface based on the gfile
% -inputs:
% path_gfile                       string,  path where the gfile is saved
% varargin                         key, if varargin='contour' plot the contour of the flux surface, else just plot the lines of the flux surface
% nw_refine/nh_refine              int, parameters which are used to refine the grids and varialbes
% position_coil                    matrix, position of the PF coils (the columns are:R, Z, Rwidth, Zwidth, flag for coil_grid_r_z.m, turns in R direction (used just for the number of grids in coil_grid_r_z.m,which can be arbitrary value), turns in Z direction, total turns.)
% vv                               vector, the dimensions of the VV, [defaut for EXL: 0.158 1.655 -1.405 1.405]
%rlim/zlimter

% -outputs:
% gfile              struct, the structure variable of gfile
% plot of the flux surface based on the gfile

% Edited by Bin Chen in 2019/06/13
% Contact: chenbino@enn.cn
% ENN Group 1989-2019, all rights reserved.

gfile=read_gfile_struct(path_gfile); % create a structure variable for the gfile

% refine the grids and varialbes
nw_refine=129; %refine the grids
nh_refine=129;
rgrid_refine=linspace(gfile.rleft,gfile.rleft+gfile.rdim,nw_refine);
zgrid_refine=linspace(gfile.zmid-gfile.zdim/2,gfile.zmid+gfile.zdim/2,nh_refine);
[rmesh_refine,zmesh_refine]=meshgrid(rgrid_refine,zgrid_refine);
flux_efit_refine=imresize(gfile.flux_efit0,[nw_refine nh_refine],'bicubic');

hold off;
plot(gfile.rmaxis,gfile.zmaxis,'r+','linewidth',14);

%fig=figure(1);
%set(fig,'units','normalized','position',[0.1 0.05 0.4 0.87]);
if strcmpi(varargin, 'contour')
    if gfile.sibry<0
        flux_surface=flux_efit_refine'-gfile.sibry;
        hold on;
        surf(rmesh_refine,zmesh_refine,min(min(flux_surface))*ones(size(rmesh_refine))*1.05,ones(size(rmesh_refine))*1.0);
        surf(rmesh_refine,zmesh_refine,flux_surface);
        colormap(flipud(hot));
        shading interp;
        grid off
        box off
        view([0,0,1]);
        zlim([min(min(flux_surface))*1.1,0]);
        caxis([min(min(flux_surface))*1.1,0]);
    else
        flux_surface=flux_efit_refine';
        hold on;
        pcolor(rmesh_refine,zmesh_refine,flux_surface);
        colormap(flipud(hot));
        shading interp;
        grid off
        box off
        caxis([min(min(flux_surface)),max(max(flux_surface))]);
    end
end
hold on;
contour(rmesh_refine,zmesh_refine,flux_efit_refine',80,'Color',[0.5 0.5 0.5]);   
%%%%%%%%%%%%%%%%%%%if configure divertor%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % plot(gfile.rbbbs,gfile.zbbbs,'r-','linewidth',2);
contour(rmesh_refine,zmesh_refine,flux_efit_refine',[gfile.sibry,gfile.sibry],'r','linewidth',2);
set(gca,'fontweight','bold','fontsize',16,'fontname','times','XTick',[],'YTick',[],'Color', 'None');
axis equal;
axis([0 2.2 -2.0 2.0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%if configure limiter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  [cc,hh]=contour(rmesh_refine,zmesh_refine,flux_efit_refine',[gfile.sibry,gfile.sibry],'w');
% %     
% %     lentmp=size(cc);
% %     count=1:3;
% %     sumtmp=1;
% %     for j=1:count(end)
% %         if sumtmp>lentmp(2)
% %               continue
% %         end 
% %          count(j)=cc(2,sumtmp);
% %          tmpx=cc(1,sumtmp+1:sumtmp+count(j));
% %          tmpy=cc(2,sumtmp+1:sumtmp+count(j));
% %          tmpz=tmpx*0.0;
% %          if max(tmpy)<max(zlimter)+0.1 && min(tmpy)>min(zlimter)-0.1;
% %            plot3(tmpx,tmpy,tmpz,'r','linewidth',3.0);
% %          end
% %          sumtmp=sumtmp+count(j)+1
% %     end

fclose all;
end