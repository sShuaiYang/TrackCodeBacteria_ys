function imageStack= import_tiff_stack( fname )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
warning off all
infoImage=imfinfo(fname);
frameNum=size(infoImage,1);
imageWith=infoImage(1).Width;
imageHeight=infoImage(1).Height;
imageBit=infoImage(1).BitsPerSample;
imageBit=strcat('uint',num2str(imageBit));
imageStack=zeros(imageHeight,imageWith,imageBit);
imageCurrent=Tiff(fname,'r');
for iframe=1:frameNum
    imageCurrent.setDirectory(iframe);
    imageStack(:,:,iframe)= imageCurrent.read();
end
imageCurrent.close();
end
function [afterProcessingImages,imageProcessingInfo]=myImageProcessing(beforeProcessingImages,picType) %this function can segreate orignal images as you want and returen the mask images
imageType='uint8'; %here you can chenge your image type
cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
% edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
afterProcessingImages1=zeros(size(beforeProcessingImages),imageType);
% here you can change the image parameters
grayThresh=30;
areaThreshold=40;  % proper 60, overlap 60+40
maxIntensity=235;      %%%%%%% long time 100X neo

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
cropArea=~cropInfo;

for iframe=1:size(beforeProcessingImages,3)
    % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter);
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter);
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove');
    afterProcessingImages(:,:,iframe)=imdilate(afterProcessingImages(:,:,iframe),ones(2));
%     afterProcessingImages(:,:,iframe)=afterProcessingImages(:,:,iframe)-afterProcessingImages1(:,:,iframe);
end
afterProcessingImages=logical(afterProcessingImages);
end
% function imageStack= import_tiff_stack( fname )
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% warning off all
% infoImage=imfinfo(fname);
% frameNum=size(infoImage,1);
% imageWith=infoImage(1).Width;
% imageHeight=infoImage(1).Height;
% imageBit=infoImage(1).BitsPerSample;
% imageBit=strcat('uint',num2str(imageBit));
% imageStack=zeros(imageHeight,imageWith,3,'uint8');
% imageCurrent=Tiff(fname,'r');
% for iframe=1:frameNum
%     imageCurrent.setDirectory(iframe);
%     imageStack(:,:,:,iframe)= imageCurrent.read();
% end
% imageCurrent.close();
% end
