function A = read_gfile_struct(path_gfile)
%% read the EFIT output data from gfile
% -inputs:
% path_gfile    string, the path of the gfile
% -outputs:
% A             struct, the data passed from gfile

% Edited by Chenbin in 2018/10/25
% Modified by Shuying SUN in 2019/05/10

fid=fopen(path_gfile);
% read the file header, 40 characters
fscanf(fid,'%c',40);
% read the simulation grid info, 3 float numbers
temp=fscanf(fid,'%g',3);
A.nw=temp(2);   % number of horizontal R grid points
A.nh=temp(3);   % number of vertical Z grid points

temp=fscanf(fid,'%g',4*5);
A.rdim  =temp(1);   % Horizontal dimension in meter of computational box
A.zdim  =temp(2);   % Vertical dimension in meter of computational box
A.rcentr=temp(3);   % R in meter of vacuum toroidal magnetic field BCENTR
A.rleft =temp(4);   % Minimum R in meter of rectangular computational box
A.zmid  =temp(5);   % Z of center of computational box in meter
A.rmaxis=temp(6);   % R of magnetic axis in meter
A.zmaxis=temp(7);   % Z of magnetic axis in meter
A.simag =temp(8);   % poloidal flux at magnetic axis in Weber /rad
A.sibry =temp(9);   % poloidal flux at the plasma boundary in Weber /rad
A.bcentr=temp(10);  % Vacuum toroidal magnetic field in Tesla at RCENTR
A.current=temp(11); % Plasma current in Ampere
A.simag =temp(12);
A.xdum  =temp(13);
A.rmaxis=temp(14);
A.xdum  =temp(15);
A.zmaxis=temp(16);
A.xdum  =temp(17);
A.sibry =temp(18);
A.xdum  =temp(19);
A.xdum  =temp(20);
A.fpol=fscanf(fid,'%g',A.nw);   % Poloidal current function in m-T, F = RB_T on flux grid
A.pres=fscanf(fid,'%g',A.nw);   % Plasma pressure in nt / m on uniform flux grid
A.ffprim=fscanf(fid,'%g',A.nw); % FF' in (mT)^2 / (Weber /rad) on uniform flux grid
A.pprime=fscanf(fid,'%g',A.nw); % P' in (nt /m2) / (Weber /rad) on uniform flux grid
% read 2D flux mesh
% Poloidal flux in Weber / rad on the rectangular grid points
temp=fscanf(fid,'%g',A.nw*A.nh);
% A.flux=reshape(temp,A.nw,A.nh);
for i=1:A.nw
    flux_efit0(:,i)=temp((i-1)*A.nh+1:i*A.nh);
end
A.flux_efit0=flux_efit0;
A.psirz=flux_efit0;
% read the q profile against psi
% q values on uniform flux grid from axis to boundary
A.qpsi=fscanf(fid,'%g',A.nw);
% read the boundary and limiter parameters
temp=fscanf(fid,'%g',2);
A.nbbbs =temp(1);   % Number of boundary points
A.nlimtr=temp(2);   % Number of limiter points
temp=fscanf(fid,'%g %g',[2,A.nbbbs]);
A.rbbbs=temp(1,:);  % R of boundary points in meter
A.zbbbs=temp(2,:);  % Z of boundary points in meter
temp=fscanf(fid,'%g %g',[2,A.nlimtr]);
A.rlim=temp(1,:);   % R of surrounding limiter contour in meter
A.zlim=temp(2,:);   % Z of surrounding limiter contour in meter
%%
A.rgrid=linspace(A.rleft,A.rleft+A.rdim,A.nw);
A.zgrid=linspace(A.zmid-A.zdim/2,A.zmid+A.zdim/2,A.nh);
[A.rmesh,A.zmesh]=meshgrid(A.rgrid,A.zgrid);
fclose(fid);
