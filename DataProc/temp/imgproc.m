function y = imgproc
img1=imread('songxm1.jpg');
img1=rgb2gray(img1);
img2=imread('songxm1.jpg');
img2=rgb2gray(img2);
g=normxcorr2(img1,img2);







return
isCCD=1;
if isCCD
    load('img_data');
%     load('img2')
    
%     myMaxInteger=max(img(:));
%     k=255/myMaxInteger;
%     img=(img*k);
else
    img=imread('songxm.jpg');
    img=rgb2gray(img);
end

%% img is grayscale
imshow(img,[])
pause(1)
hold on

imgSum=sum(img,1);
[s1,idx1]=min(diff(imgSum));
[s2,idx2]=max(diff(imgSum));


myY=get(gca,'ylim');
line([idx1,idx1],myY,'Color','r')
line([idx2,idx2],myY)

img(:,idx1:idx2)=0;
imshow(img,[])
pause(1)

isOTSU=1;
if isOTSU
    thresh=multithresh(img,2);
%     img=imquantize(img,thresh);
end






% top threshold

idx=img>thresh(2);
img(:,:)=0;
img(idx)=1;

contour(img)
imshow(img,[])
pause(1)
% median filter
img = medfilt2(img);




imshow(img,[])
pause(1)
img=bwareaopen(img,50);
imshow(img,[])
pause(1)





img=bwmorph(img,'clean');
imshow(img,[])
pause(1)

img=bwmorph(img,'fill');
imshow(img,[])
pause(1)



% 
% img=imclearborder(img);
% imshow(img,[])
% pause(1)

%% morphology
isMorphology=1;

% se=offsetstrel('ball',5,6);
se=strel('disk',5);




for i=1:1
    
if isMorphology
    imgOPEN=imopen(img,se);
    imgCLOSE=imclose(img,se);
    
    img=(imgOPEN+imgCLOSE)/2;
end
imshow(img,[])
drawnow
pause(0.1)

end



% for i=1:5
% img=bwmorph(img,'skel');
% imshow(img,[])
% end

% img=imerode(img,se);    
% imshow(img,[])

% CC=bwconncomp(img,8);
% stats=regionprops(CC,'Perimeter','PixelIdxList','ConvexHull');

C = contourc(double(img),[1 1]); %limiter contour
m=C(2,1);

C(:,1)=[];

x = C(1,1:m);
y = C(2,1:m);

% x0=idx2;
% y0=size(img,1)/2;
% 
% L=sqrt((x-x0).^2+(y-y0).^2);
% 
% [minY,idx]=min(L);
[minY,idx]=min(y);

x=circshift(x,idx,2);

y=circshift(y,idx,2);


%%
windowWidth = 31;
polynomialOrder = 2;
overlap=20;
x(end+1:end+overlap)=x(1:overlap);
y(end+1:end+overlap)=y(1:overlap);
x1 = sgolayfilt(x, polynomialOrder, windowWidth); % Curvefit of the discrete points in Img
y1 = sgolayfilt(y, polynomialOrder, windowWidth);

% x2=smooth(x,windowWidth,'sgolay');
% y2=smooth(y,windowWidth,'sgolay');
% x3=smooth(x,windowWidth,'lowess');
% y3=smooth(y,windowWidth,'lowess');
% x4=smooth(x,windowWidth,'rlowess');
% y4=smooth(y,windowWidth,'rlowess');
% x5=smooth(x,windowWidth,'rloess');
% y5=smooth(y,windowWidth,'rloess');
% 
% 
% 
% figure;hold on;plot(x1,'.r');plot(x2,'.b');plot(x3,'.m');plot(x4,'.g');plot(x5,'.c')
% figure;hold on;plot(y1,'.r');plot(y2,'.b');plot(y3,'.m');plot(y4,'.g');plot(y5,'.c')
% 
% 
% figure;hold on;plot(y,'.r');plot(y1,'.b');plot(y2,'.b')

boundx=x1(1+overlap/2:end-overlap/2);
boundy=y1(1+overlap/2:end-overlap/2);
hold on
plot(boundx,boundy,'.r')


return

for i=1:100
img=imdilate(img,se);    
imshow(img,[])
drawnow
pause(0.1)
end


%% ultimate erosion
BW2=bwulterode(img);
imshow(BW2,[])
pause(1)





% CC = bwconncomp(img,4);
% 


kernel=[1 1 1;1 1 1;1 1 1]/9;
for i=1:1
img=conv2(img,kernel,'same');
end

%% morphology
isMorphology=1;
if isMorphology
    se=offsetstrel('ball',5,6);
    imgOPEN=imopen(img,se);
    imgCLOSE=imclose(img,se);
    
    img=(imgOPEN+imgCLOSE)/2;
end
imshow(img,[])
pause(1)

%% sobel method
kernel = [1 2 1; 0 0 0; -1 -2 -1];


for i=1:1
    imgH=conv2(img,kernel,'same');
    imgV=conv2(img,kernel','same');
    img=sqrt(imgH.*imgH + imgV.*imgV)/sqrt(2);
end

imshow(img,[])
pause(1)



ind=img>20;
img(ind)=1;
img(~ind)=0;
imshow(img,[])
CC=bwconncomp(img,8);
stats=regionprops(CC,'Perimeter','PixelIdxList');
hold on
img(:)=0;

myPerimeter=cat(1,stats.Perimeter);
[maxL,idx]=sort(myPerimeter,'descend');
img(stats(idx(1)).PixelIdxList)=1;
img(stats(idx(2)).PixelIdxList)=1;
img(stats(idx(3)).PixelIdxList)=1;

imshow(img,[])
pause(1)

isOTHER=0;
if isOTHER
    [imgLOG,thresh]=edge(img,'log',[],2.3);
    imshow(imgLOG)
    
    imgWater=watershed(imgLOG,8);
    imshow(imgWater)
    
    k=graythresh(img/256);
    
    imgOTSU=im2bw(img,k);
    
    imshow(imgOTSU,[])
end






return

kernal=[1 1 1;1 -8 1;1 1 1];
img=conv2(img,kernel,'same');
image(img)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
return

% figure;imshow(img,[])

% image(img)
if 0
    kernal=[1 1 1;1 -8 1;1 1 1];
        img=conv2(img,kernel,'same');
    image(img)
    return

end

% figure;imshow(img,[])
% 
% 
% img=imread('songxm.jpg');
% image(img)
if 1
    kernal=[1 1 1;1 -8 1;1 1 1];
    for i=1:3
%         img(:,:,i)=conv2(img(:,:,i),kernal,'same');
        img(:,:,i)=conv2(img(:,:,i),kernal,'same');
    end
    image(img)
end
%% sobel method
kernal = [1 2 1; 0 0 0; -1 -2 -1];
for i=1:3
    imgH=conv2(img(:,:,i),kernal,'same');
    imgV=conv2(img(:,:,i),kernal','same');
    img(:,:,i)=sqrt(imgH.*imgH + imgV.*imgV)/sqrt(2);
end
image(img)
end

