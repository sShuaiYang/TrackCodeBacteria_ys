function maskImage=gfpChannelImageProcessing(image)
% find those
backGround=getBackGround(image(:,:,1));
image=image-backGround;
[maskImage,removeMask,para]=imageProcessing(image);
dirFile='F:\2013-09-30 jzy F1-T-dimer-EGFP longTime\fluoImage\image1';
mkdir(dirFile);
cd(dirFile)
save('paraImageProcessing','para')
getOrderNum(maskImage,image,removeMask);
end
function backGround=getBackGround(image)
image=image(:);
image=sort(image);
imageSize=numel(image);
backGround=mean(image(1:imageSize/10));
end
function [afterProcessingImages,afterProcessingImages1,para]=imageProcessing(beforeProcessingImages)
gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),'uint16');% here intiatlize the maskImages stacks
grayThresh=80;
areaThreshold=100;  % proper 60, overlap 60+40
maxIntensity=200;      %%%%%%% long time 100X neo

para.grayThresh=grayThresh;
para.areaThreshold=areaThreshold;
para.maxIntensity=maxIntensity;

parfor iframe=1:size(beforeProcessingImages,3)
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/65535);%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
%     %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open');
%     %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'bridge'); % bridge process
    image=afterProcessingImages(:,:,iframe);
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    afterProcessingImages(:,:,iframe)=image;
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
% %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes');
    afterProcessingImages1(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % fill holes process
end
afterProcessingImages=logical(afterProcessingImages);
afterProcessingImages1=logical(afterProcessingImages1);
end
function getOrderNum(maskImage,image,removeMask)
for i=1:size(maskImage,3)
    cc=regionprops(maskImage(:,:,i),'Centroid');
    imshow(segrationImage(imadjust(image(:,:,i)),removeMask(:,:,i)));
    hold on
    for iCC=1:numel(cc)
        text(cc(iCC).Centroid(1),cc(iCC).Centroid(2),num2str(iCC),'Color','r','FontSize',10);
    end
    saveas(gcf,num2str(i),'tif');
    close all
end
end
function  afterProcessingImages=segrationImage(beforeProcessingImages,maskImages) % this function can overlay your mask and orignal images 
maskColor=[0,1,0];
imageType='uint8'; %here you can chenge your image type
afterProcessingImages=zeros(size(beforeProcessingImages,1),size(beforeProcessingImages,2),3,size(beforeProcessingImages,3),imageType);
parfor iframe=1:size(beforeProcessingImages,3)
%     afterProcessingImages(:,:,iframe)=immultiply((~bwperim(maskImages(:,:,iframe),4)),beforeProcessingImages(:,:,iframe));
%      afterProcessingImages(:,:,iframe)=immultiply(~maskImages(:,:,iframe),beforeProcessingImages(:,:,iframe));
afterProcessingImages(:,:,:,iframe)=imoverlay(beforeProcessingImages(:,:,iframe),maskImages(:,:,iframe),maskColor);
end
end