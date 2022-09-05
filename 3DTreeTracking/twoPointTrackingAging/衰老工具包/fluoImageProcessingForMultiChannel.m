function fluoImageProcessingForMultiChannel()
% protocal contain xyScan(field>=1) for GFP,RFP and brightFiled
% xyScan field has been divided, but the file number depend on the imagenum
dirFile=uigetdir();
fileNum=2;    %%% ceil(allFrame*3/189)   each field has how many files
dirOriginalData=strcat(dirFile,'\original data');
nameList=dir(dirOriginalData);
dirTiff=strcat(dirFile,'\tiffData');
mkdir(dirTiff);
backGround=import_tiff_stack(strcat(dirFile,'\backGround.tif'));
for i=1:(numel(nameList)-2)/fileNum
    smallFile=strcat(dirTiff,'\',num2str(i));
    mkdir(smallFile);
    cd(smallFile);
    dirFileName=[];
    for sFile=1:fileNum
        dirFileName{sFile}=strcat(dirOriginalData,'\',nameList(2+(i-1)*fileNum+sFile).name);
    end
    [RFPImage,GFPImage,BFImage]=loadImage(dirFileName,fileNum);
%     [newRFPImage,newGFPImage]=generateAdjustImage(RFPImage,GFPImage);
    save('RFP.mat','RFPImage');
    save('GFP.mat','GFPImage');
%     [result,BFImage]=quickImageAve(GFPImage,RFPImage,BFImage,backGround);
%     save('BF.mat','BFImage')
%     save('result.mat','result')
end
end
function [RFPImage,GFPImage,BFImage]=loadImage(dirFileName,fileNumber)
finalImage=import_tiff_stack(dirFileName{1});
for i=2:fileNumber
    finalImage1=import_tiff_stack(dirFileName{i});
    finalImage=cat(3,finalImage,finalImage1);
end
GFPImage=finalImage(:,:,1:3:end);
RFPImage=finalImage(:,:,2:3:end);
BFImage=finalImage(:,:,3:3:end);
end
function [newRFPImage,newGFPImage]=generateAdjustImage(RFPImage,GFPImage)
gfpMax=max(max(max(GFPImage)));
gfpMin=min(min(min(GFPImage)));
rfpMax=max(max(max(RFPImage)));
rfpMin=min(min(min(RFPImage)));
for i=1:size(RFPImage,3)
    newRFPImage(:,:,i)=imadjust(RFPImage(:,:,i),[double(rfpMin)/65535*1.2;double(rfpMax)/65535/1.1],[0;1]);
    newGFPImage(:,:,i)=imadjust(GFPImage(:,:,i),[double(gfpMin)/65535*1.2;double(gfpMax)/65535/1.1],[0;1]);
end
end
function [result,BFImage]=quickImageAve(GFPImage,RFPImage,BFImage,backGround)
BFImage=uint16(double(BFImage)/double(max(max(max(BFImage))))*0.5*2048);
backGround=uint16(double(backGround)/double(max(max(max(backGround))))*0.5*2048);
BFImage=backGroundCorrection(BFImage,backGround,'16bit');
brightMask=myImageProcessing(BFImage);
result(:,1)=getResult(brightMask,RFPImage);
result(:,2)=getResult(brightMask,GFPImage);
end
function result=getResult(brightMask,image)
for i=1:size(image,3)
    backGround=getBackGround(image(:,:,i));
    image(:,:,i)=image(:,:,i)-backGround;
    operateImage=image(:,:,i);
    result(i,1)=mean(operateImage(brightMask(:,:,i)==1));
end
end
function backGround=getBackGround(image)
image=image(:);
image=sort(image);
imageSize=numel(image);
backGround=mean(image(1:imageSize/10));
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
beforeProcessingImages=255-beforeProcessingImages;
grayThresh=80;
areaThreshold=60;  % proper 60, overlap 60+40
maxIntensity=240;      %%%%%%% long time 100X neo
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
    afterProcessingImages(:,:,iframe)=image;
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes');
end
afterProcessingImages=logical(afterProcessingImages);
end