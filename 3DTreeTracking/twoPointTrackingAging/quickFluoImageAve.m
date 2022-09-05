function [result,brightFieldImage]=quickFluoImageAve(image,dirFile,brightFieldImage)
if nargin==2
    brightFieldImage=loadingBFImage(dirFile);
    backGround=import_tiff_stack(strcat(dirFile,'\backGround.tif'));
    brightFieldImage=uint16(double(brightFieldImage)/65535*2048);
    backGround=uint16(double(backGround)/65535*2048);
    brightFieldImage=backGroundCorrection(brightFieldImage,backGround,'16bit');
end
brightMask=myImageProcessing(brightFieldImage);
for i=1:size(image,3)
    backGround=getBackGround(image(:,:,i));
    image(:,:,i)=image(:,:,i)-backGround;
    operateImage=image(:,:,i);
    result(i)=mean(operateImage(brightMask(:,:,i)==1));
end
end
function backGround=getBackGround(image)
image=image(:);
image=sort(image);
imageSize=numel(image);
backGround=mean(image(1:imageSize/5));
end
function brightFieldImage=loadingBFImage(dirFile)
dirFile1=strcat(dirFile,'\original data');
nameList=dir(dirFile1);
for i=1:numel(nameList)-2
    fileName=strcat(dirFile1,'\',nameList(i+2).name);
    brightFieldImage(:,:,i)=import_one_tiff(fileName);
end
finalFrame=i;
dirFile2=strcat(dirFile,'\information\otherImage');
nameList=dir(dirFile2);
for i=1:numel(nameList)-2
    fileName=strcat(dirFile2,'\',nameList(i+2).name);
    brightFieldImage(:,:,i+finalFrame)=import_one_tiff(fileName);
end
end
function afterProcessingImages=myImageProcessing(beforeProcessingImages) %this function can segreate orignal images as you want and returen the mask images
imageType='uint8'; %here you can chenge your image type
gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=90;
areaThreshold=60;  % proper 60, overlap 60+40
maxIntensity=220;      %%%%%%% long time 100X neo
parfor iframe=1:size(beforeProcessingImages,3)
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
        image=afterProcessingImages(:,:,iframe);
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList','MeanIntensity');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity 
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    cc=regionprops(logical(1-afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'PixelIdxList','MeanIntensity','FilledArea');
    for iCC=1:size(cc,1)
        if cc(iCC).MeanIntensity==255 && cc(iCC).FilledArea<=areaThreshold
            image(cc(iCC).PixelIdxList)=1;
        end
    end
    afterProcessingImages(:,:,iframe)=image;
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
end
afterProcessingImages=logical(afterProcessingImages);
parfor i=1:size(afterProcessingImages,3)
    afterProcessingImages(:,:,i)=imclearborder(afterProcessingImages(:,:,i));
    cc=regionprops(logical(afterProcessingImages(:,:,i)),beforeProcessingImages(:,:,i),'MaxIntensity','FilledArea','PixelIdxList');
    image=afterProcessingImages(:,:,i);
    for iCC=1:size(cc,1)
        if (cc(iCC).FilledArea<=areaThreshold+40 && cc(iCC).MaxIntensity<=maxIntensity) || cc(iCC).FilledArea<=areaThreshold || cc(iCC).MaxIntensity<=120
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    afterProcessingImages(:,:,i)=image;
end
end
