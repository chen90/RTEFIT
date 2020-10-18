function [ze_max,times]=edges(shotnum,t1,t2,dframe)
% This Function is used to compute the plasma surface LCSF, In this
% function, the camera lens distortion and position is considered, And the
% coordinate transforma from  Cartesian to cylindrical is computed.
% CAUTION£¡£¡£¡ The calibration file should be put in right path
if nargin<4
    dframe=1;
end
[info,frame,times]=downloadcine(shotnum,t1,t2,dframe);
M120_cal_path='G:\matlab code\ENN\M120_Calibration.mat'; %You should correct it depending on your own path
load(M120_cal_path);
for i=1:size(frame,1)
waitbar(i/size(frame,1));
img=reshape(frame(i,:,:),info.Height,info.Width); 
%-----Use the Calibration parameters to correct the image distortion------
img=undistortImage(img,cameraParams);
threshold=find_threshold(img,0); %Use double Gaussian fitting to get the proper luminess threshold
%%
ref_color=0;                                  % the backbround lightness is 67
bin_mask = zeros(info.Height,info.Width);
bin_mask = bin_mask | (img - ref_color).^2 <= threshold^2; % MagicWand method

Blobs=bwconncomp(~bin_mask);                   % Find the connected blobs in 2-D image
numPixels=cellfun(@numel,Blobs.PixelIdxList);  % Get the pixel numbers of each blob
[~,idx]=max(numPixels);                        % The index of the biggest Blob
biggest_Blob=zeros(info.Height,info.Width);
biggest_Blob(Blobs.PixelIdxList{idx})=1;       % biggest_Blob is the selection of biggest blob
boundaries = bwboundaries(biggest_Blob);       % Get the boundary of biggest blob

firstBoundary = boundaries{1};
x = firstBoundary(:, 2); 
y = firstBoundary(:, 1); 
%%
windowWidth = 199;
polynomialOrder = 2;
overlap=100;
x(end+1:end+overlap)=x(1:overlap);
y(end+1:end+overlap)=y(1:overlap);
x1 = sgolayfilt(x, polynomialOrder, windowWidth); % Curvefit of the discrete points in Img
y1 = sgolayfilt(y, polynomialOrder, windowWidth);
boundx=x1(1+overlap/2:end-overlap/2);
boundy=y1(1+overlap/2:end-overlap/2);

%%
cx=664;                  % Horizontal pixel center
cy=391;                  % Vertical pixel center
ratio=4.1/1e3;           % unit mm/pixel;
focalLength=14.7;        % unit mm
distance=2.020;          % unit mm
px=(boundx-cx)*ratio;    % Image x axis in meter unit
py=(boundy-cy)*ratio;    % Image x axis in meter unit
%%
center_width=0.5;        % Select part of the Curve to avoid the Center Rod
ue=px(px>center_width);
ve=py(px>center_width);
uc=0;                    % World coordinate of camera
vc=0;
wc=2.02;
%---------Coordinate transform From Cartesian to cylindrical------------
duedve=diff(ue)./diff(ve);
duedve=[duedve;duedve(end)];
phie=atan(((ue-uc)-(ve-vc).*duedve)./wc);
Re=ue.*wc./(wc*cos(phie)+(ue-uc).*sin(phie));
Ze=ve-(ve-vc).*(Re.*sin(phie))./wc;
Re2=Re(round(0.3*length(Re)):round(0.7*length(Re)));
sort_re=sort(Re2,'descend');
ze_max(i)=mean(sort_re(1:5));
end
end
