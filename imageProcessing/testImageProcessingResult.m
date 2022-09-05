function [mergeImage,imageGet,myBackGround]=testImageProcessingResult(varargin)
% this funcion has lots of benifits
% varargin could be (),imageGet,'backGround',which means

% #################
% varargin=()
% you will choose a document where the original data locates
% and the you will get the afterProcessingImage--mergeImage
% use implay function to view it
% 
% ###########
% varargin='backGround';
% then this function will use the backGrond information we had get to create a new backGround image
% 
% ##########
% varargin=imageGet
% which means you have run this function more than one time and you have got the original picture,this will help you get you answer more quickly

imageGet=[];
if size(varargin,2)~=0
    if size(varargin{1},3)~=0
        imageGet=varargin{1};
    end
end
if isempty(imageGet)
    dataDir=uigetdir();
    imageDir=strcat(dataDir,'\original data');
    cd(imageDir);
    nameList=dir(imageDir);
    imageGet=[];
    backGround=import_one_tiff(strcat(dataDir,'\','backGround.tif'));
    for i=1:1:size(nameList,1)-2
        disp(i)
        nameInfo=nameList(i+2,1).name;
        fname=strcat(imageDir,'\',nameInfo);
        imageGet(:,:,i)=import_one_tiff(fname);
    end
    imageGet=uint16(imageGet);
end
if size(varargin,2)~=0 && strcmp(varargin{1},'backGround')==1
    rudeImage=backGroundCorrection(imageGet,backGround,'16bit');
    rudeMask=myImageProcessingforOriData(rudeImage);
    myBackGround=uint16(zeros(size(rudeMask,1),size(rudeMask,2)));
    myMask=[];
    for i=1:size(rudeMask,3)
        if i==1
            currentImage=imageGet(:,:,i);
            myBackGround(rudeMask(:,:,1)==0)=currentImage(rudeMask(:,:,1)==0);
            myMask=~rudeMask(:,:,1);
        else
            currentImage=imageGet(:,:,i);
            overlapMask=~rudeMask(:,:,i) & ~myMask;
            myBackGround(overlapMask==1)=currentImage(overlapMask==1);
            myMask=~rudeMask(:,:,i) | myMask;
        end
        if min(min(myMask))==1
            break
        end
    end
    imageGet=backGroundCorrection(imageGet,myBackGround,'16bit');
else
    if size(varargin,2)==0
        backGround=uint16(double(backGround)/65535*2048);
        imageGet=uint16(double(imageGet)/65535*2048);
%         backGround=backGround(1200:1200+723,1605:1605+723);
%         imageGet=imageGet(1200:1200+723,1605:1605+723,:,:);
        imageGet=backGroundCorrection(imageGet,backGround,'16bit');
    end
    myBackGround=[];
end
processingImage=myImageProcessing(imageGet,'w');
mergeImage=segrationImage(imageGet,processingImage);
end
function afterProcessingImages=myImageProcessingforOriData(beforeProcessingImages) % this function can segreate orignal images as you want and returen the mask images
imageType='uint8'; %here you can chenge your image type
gaussianFilter=fspecial('gaussian',[5, 5],5); %here create Gaussian Blur Filter
% edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=75;
areaThreshold=100;
for iframe=1:size(beforeProcessingImages,3)
    %     afterProcessingImages(:,:,iframe)=imadjust(beforeProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    image=afterProcessingImages(:,:,iframe);
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=220
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    afterProcessingImages(:,:,iframe)=image;
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes'); % fill holes process
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate',8); % find outLine process
end
afterProcessingImages=logical(afterProcessingImages);
end
