function B3TransformInfo=imtransformFluoImage2FL(dirFile)

disp('Used to transform mScarletI Venus TDsmURFP sfGFP image to CyOFP.');
ZylaPointsCyOFP2mScarletI=[];ZylaPointsmScarletI=[];
ZylaPointsCyOFP2Venus=[];ZylaPointsVenus=[];
ZylaPointsCyOFP2TDsmURFP=[];ZylaPointsTDsmURFP=[];
ZylaPointsCyOFP2sfGFP=[];ZylaPointssfGFP=[];
for i=1:1
imageCyOFP=importTiff(strcat(dirFile,'\feild',num2str(i,'%.4d'),'\CyOFP\imageCyOFP00001.tif'));
maskImageCyOFP=myImageProcessingFL1(imageCyOFP,'8bit');
imageRFP=importTiff(strcat(dirFile,'\feild',num2str(i,'%.4d'),'\mScarletI\imagemScarletI00001.tif'));
maskImagemScarletI=myImageProcessingFL2(imageRFP,'16bit');
imageVenus=importTiff(strcat(dirFile,'\feild',num2str(i,'%.4d'),'\Venus\imageVenus00001.tif'));
maskImageVenus=myImageProcessingFL2(imageVenus,'16bit');
imageTDsmURFP=importTiff(strcat(dirFile,'\feild',num2str(i,'%.4d'),'\TDsmURFP\imageTDsmURFP00001.tif'));
maskImageTDsmURFP=myImageProcessingFL2(imageTDsmURFP,'16bit');
imagesfGFP=importTiff(strcat(dirFile,'\feild',num2str(i,'%.4d'),'\sfGFP\imagesfGFP00001.tif'));
maskImagemsfGFP=myImageProcessingFL2(imagesfGFP,'16bit');

[ZylaPointsCyOFP2mScarletI,ZylaPointsmScarletI,bestPositionCyOFP2mScarletI]=getCalibrationPoints(ZylaPointsCyOFP2mScarletI,ZylaPointsmScarletI,maskImageCyOFP,maskImagemScarletI);
[ZylaPointsCyOFP2Venus,ZylaPointsVenus,bestPositionCyOFP2Venus]=getCalibrationPoints(ZylaPointsCyOFP2Venus,ZylaPointsVenus,maskImageCyOFP,maskImageVenus);
[ZylaPointsCyOFP2TDsmURFP,ZylaPointsTDsmURFP,bestPositionCyOFP2TDsmURFP]=getCalibrationPoints(ZylaPointsCyOFP2TDsmURFP,ZylaPointsTDsmURFP,maskImageCyOFP,maskImageTDsmURFP);
[ZylaPointsCyOFP2sfGFP,ZylaPointssfGFP,bestPositionCyOFP2sfGFP]=getCalibrationPoints(ZylaPointsCyOFP2sfGFP,ZylaPointssfGFP,maskImageCyOFP,maskImagemsfGFP);

end

transformInfoZylaCyOFP2mScarletI=getTransformInfoandCheckResult(ZylaPointsCyOFP2mScarletI,ZylaPointsmScarletI,imageRFP,maskImageCyOFP,bestPositionCyOFP2mScarletI);
B3TransformInfo.transformInfoZylaCyOFP2mScarletI=transformInfoZylaCyOFP2mScarletI;
B3TransformInfo.bestPositionCyOFP2mScarletI=bestPositionCyOFP2mScarletI;


transformInfoZylaCyOFP2Venus=getTransformInfoandCheckResult(ZylaPointsCyOFP2Venus,ZylaPointsVenus,imageVenus,maskImageCyOFP,bestPositionCyOFP2Venus);
B3TransformInfo.transformInfoZylaCyOFP2Venus=transformInfoZylaCyOFP2Venus;
B3TransformInfo.bestPositionCyOFP2Venus=bestPositionCyOFP2Venus;

transformInfoZylaCyOFP2TDsmURFP=getTransformInfoandCheckResult(ZylaPointsCyOFP2TDsmURFP,ZylaPointsTDsmURFP,imageTDsmURFP,maskImageCyOFP,bestPositionCyOFP2TDsmURFP);
B3TransformInfo.transformInfoZylaCyOFP2TDsmURFP=transformInfoZylaCyOFP2TDsmURFP;
B3TransformInfo.bestPositionCyOFP2TDsmURFP=bestPositionCyOFP2TDsmURFP;

transformInfoZylaCyOFP2sfGFP=getTransformInfoandCheckResult(ZylaPointsCyOFP2sfGFP,ZylaPointssfGFP,imagesfGFP,maskImageCyOFP,bestPositionCyOFP2sfGFP);
B3TransformInfo.transformInfoZylaCyOFP2sfGFP=transformInfoZylaCyOFP2sfGFP;
B3TransformInfo.bestPositionCyOFP2sfGFP=bestPositionCyOFP2sfGFP;




% transformInfoZylaCyOFP2mScarletI = cp2tform(ZylaPointsmScarletI,ZylaPointsCyOFP2mScarletI,'polynomial',3);
% transformInfoB3.transformInfoZylaCyOFP2mScarletI=transformInfoZylaCyOFP2mScarletI;
% transformInfoB3.bestPositionCyOFP2mScarletI=bestPositionCyOFP2mScarletI;
% imageTran=imtransform(imageRFP,transformInfoZylaCyOFP2mScarletI,'XData',[1 size(imageRFP,1)],'YData',[1 size(imageRFP,2)]);
% imageStack=cat(3,maskImageCyOFP,imageTran);
% imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
% imageTran=imageStack(:,:,2);
% imageFLAdjust=im2uint8(imadjust(imageTran));
% maskImage=bwmorph(maskImageCyOFP,'remove');
% imageSeg=segrationImage(imageFLAdjust,maskImage);
% figure,imshow(imageSeg);



% nameList=dir(dirFile);
% for i=1:length(nameList)-2
%     if strcmp(nameList(i+2).name,'GFP')||strcmp(nameList(i+2).name,'RFP')
%         dirImage=strcat(dirFile,'\',nameList(i+2).name);
%         dirNewImage=strcat(dirFile,'\',nameList(i+2).name,'new');
%         mkdir(dirNewImage);
%         imageList=dir(dirImage);
%         for j=1:length(imageList)-2
%             imageFluo=importTiff(strcat(dirImage),'\',imageList(j+2).name);
%             [~,imageFluo]=imCameraTransform(imageBF,imageFluo,transformInfoZyla,bestPosition);
%             imwrite(imageFluo,strcat(dirNewImage,'\image',nameList(i+2).name,num2str(j,'%05d'),'.tif'),'tif');
%         end
%     end
% end
end

function transformInfoZyla=getTransformInfoandCheckResult(ZylaPoints1,ZylaPoints2,image2,maskImage1,bestPosition)

transformInfoZyla= cp2tform(ZylaPoints2,ZylaPoints1,'polynomial',3);

imageTran=imtransform(image2,transformInfoZyla,'XData',[1 size(image2,1)],'YData',[1 size(image2,2)]);
imageStack=cat(3,maskImage1,imageTran);
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
imageTran=imageStack(:,:,2);
imageFLAdjust=im2uint8(imadjust(imageTran));
maskImage=bwmorph(maskImage1,'remove');
imageSeg=segrationImage(imageFLAdjust,maskImage);
figure,imshow(imageSeg);


end

function [ZylaPoints1,ZylaPoints2,bestPosition]=getCalibrationPoints(ZylaPoints1,ZylaPoints2,maskImage1,maskImage2)
maskImageA=cat(3,maskImage1,maskImage2);
[maskImageA,bestPosition]= fluoImageCorrection(maskImageA);
% maskImage=maskImageA(:,:,2);


[ ~ , bioTree] = bacteriaTracking( maskImageA(:,:,1) , maskImageA(:,:,2) );
bioTree=bioTreeMeasure(bioTree,0,size(maskImage1,1),size(maskImage1,2));
%
count=size(ZylaPoints1,1)+1;
% count=1;
for iroot=1:size(bioTree{1}.root,2)
    is2Node=bioTree{1}.root{iroot}.is2Node;
    leafInfo=bioTree{1}.root{iroot}.leafInfo;
    if is2Node == 0 && leafInfo(1)==2
        pixelIdxList1=bioTree{1}.root{iroot}.traceInfo.pixelIdxList{1};
        pixelIdxList2=bioTree{1}.root{iroot}.traceInfo.pixelIdxList{2};
        if numel(pixelIdxList1)<=1500 && numel(pixelIdxList1)>=100
            if size(bioTree{1}.root{iroot}.traceInfo.measurment{1},1)==1&&size(bioTree{1}.root{iroot}.traceInfo.measurment{2},1)==1
                Centroid1=ceil(bioTree{1}.root{iroot}.traceInfo.measurment{1}.Centroid);
                Centroid2=ceil(bioTree{1}.root{iroot}.traceInfo.measurment{2}.Centroid);
                if pdist2(Centroid1,Centroid2)<=3
                    ZylaPoints1(count,1:2)=Centroid1;
                    ZylaPoints2(count,1:2)=Centroid2;
                    count=count+1;
                end
            end
        end
    end
    
end


end




function [afterProcessingImages,imageProcessingInfo]=myImageProcessingBF(beforeProcessingImages,picType) %this function can segreate orignal images as you want and returen the mask images
imageType='uint8'; %here you can chenge your image type
% yCrop=210:670;
% xSize=size(beforeProcessingImages,1);
% ySize=size(beforeProcessingImages,2);
cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
% edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
% broderWithe=5; % here create broder frame
% broderFrame=true(xSize,ySize);broderFrame(:,1:broderWithe)=0;broderFrame(:,ySize-broderWithe:ySize)=0;broderFrame(1:broderWithe,:)=0;broderFrame(xSize-broderWithe:xSize,:)=0;
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=120;
areaThreshold=80;
maxIntensity=254;
% minIntensity=50;

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
for iframe=1:size(beforeProcessingImages,3)
    % adjust the images contraast
    %         afterProcessingImages(:,:,iframe)=beforeProcessingImages(:,:,iframe);
    %     beforeProcessingImages(:,:,iframe)=imadjust(beforeProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    %     thresholdImage=graythresh(afterProcessingImages(:,:,iframe)); %Otu thresholding images
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
    %     Max Entropy thresholding images
    %         afterProcessingImages(:,:,iframe)=im2bw_ent(afterProcessingImages(:,:,iframe));
    %     afterProcessingImages(:,:,iframe)=afterProcessingImages(:,:,iframe)&broderFrame; % remove pixels in the broder Frame
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe)); %remove the objectives in the border
    %     L=bwlabel(afterProcessingImages(:,:,iframe));
    %     postionsObjective=regionprops(L,beforeProcessingImages(:,:,iframe),'area','centroid','MaxIntensity');
    %     if ~isempty(postionsObjective), % empty stats causes crash
    %         afterProcessingImages(:,:,iframe)=removeFakeparticle(postionsObjective',afterProcessingImages(:,:,iframe));
    %     end
    %
    
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open');
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'bridge'); % bridge process
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image=afterProcessingImages(:,:,iframe);
            image(cc(iCC).PixelIdxList)=0;
            afterProcessingImages(:,:,iframe)=image;
        end
    end
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes'); % fill holes process
    %     cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MinIntensity','PixelIdxList');
    %     for iCC=1:size(cc,1)
    %         if cc(iCC).MinIntensity<=minIntensity
    %             image=afterProcessingImages(:,:,iframe);
    %             image(cc(iCC).PixelIdxList)=0;
    %             afterProcessingImages(:,:,iframe)=image;
    %         end
    %     end
    % afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    % afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % find outLine process
end
afterProcessingImages=logical(afterProcessingImages);
end


function [afterProcessingImages,imageProcessingInfo]=myImageProcessingFL1(beforeProcessingImages,picType) %this function can segreate orignal images as you want and returen the mask images
imageType='uint16'; %here you can chenge your image type
% yCrop=210:670;
% xSize=size(beforeProcessingImages,1);
% ySize=size(beforeProcessingImages,2);
cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
gaussianFilter=fspecial('gaussian',[5, 5],10); %here create Gaussian Blur Filter
% edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
% broderWithe=5; % here create broder frame
% broderFrame=true(xSize,ySize);broderFrame(:,1:broderWithe)=0;broderFrame(:,ySize-broderWithe:ySize)=0;broderFrame(1:broderWithe,:)=0;broderFrame(xSize-broderWithe:xSize,:)=0;
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=40;
areaThreshold=50;
maxIntensity=220;
% minIntensity=50;

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
for iframe=1:size(beforeProcessingImages,3)
    % adjust the images contraast
    %         afterProcessingImages(:,:,iframe)=beforeProcessingImages(:,:,iframe);
    %     beforeProcessingImages(:,:,iframe)=imadjust(beforeProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    %     thresholdImage=graythresh(afterProcessingImages(:,:,iframe)); %Otu thresholding images
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/(2^16-1));%1.2 for the image with few particles
    %     Max Entropy thresholding images
    %         afterProcessingImages(:,:,iframe)=im2bw_ent(afterProcessingImages(:,:,iframe));
    %     afterProcessingImages(:,:,iframe)=afterProcessingImages(:,:,iframe)&broderFrame; % remove pixels in the broder Frame
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe)); %remove the objectives in the border
    %     L=bwlabel(afterProcessingImages(:,:,iframe));
    %     postionsObjective=regionprops(L,beforeProcessingImages(:,:,iframe),'area','centroid','MaxIntensity');
    %     if ~isempty(postionsObjective), % empty stats causes crash
    %         afterProcessingImages(:,:,iframe)=removeFakeparticle(postionsObjective',afterProcessingImages(:,:,iframe));
    %     end
    %
    
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open');
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'bridge'); % bridge process
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image=afterProcessingImages(:,:,iframe);
            image(cc(iCC).PixelIdxList)=0;
            afterProcessingImages(:,:,iframe)=image;
        end
    end
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes'); % fill holes process
    %     cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MinIntensity','PixelIdxList');
    %     for iCC=1:size(cc,1)
    %         if cc(iCC).MinIntensity<=minIntensity
    %             image=afterProcessingImages(:,:,iframe);
    %             image(cc(iCC).PixelIdxList)=0;
    %             afterProcessingImages(:,:,iframe)=image;
    %         end
    %     end
    % afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    % afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % find outLine process
end
afterProcessingImages=logical(afterProcessingImages);
end

function [afterProcessingImages,imageProcessingInfo]=myImageProcessingFL2(beforeProcessingImages,picType) %this function can segreate orignal images as you want and returen the mask images
imageType='uint16'; %here you can chenge your image type
% yCrop=210:670;
% xSize=size(beforeProcessingImages,1);
% ySize=size(beforeProcessingImages,2);
cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
gaussianFilter=fspecial('gaussian',[5, 5],10); %here create Gaussian Blur Filter
% edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
% broderWithe=5; % here create broder frame
% broderFrame=true(xSize,ySize);broderFrame(:,1:broderWithe)=0;broderFrame(:,ySize-broderWithe:ySize)=0;broderFrame(1:broderWithe,:)=0;broderFrame(xSize-broderWithe:xSize,:)=0;
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=10;
areaThreshold=50;
maxIntensity=100;
% minIntensity=50;

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
for iframe=1:size(beforeProcessingImages,3)
    % adjust the images contraast
    %         afterProcessingImages(:,:,iframe)=beforeProcessingImages(:,:,iframe);
    %     beforeProcessingImages(:,:,iframe)=imadjust(beforeProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    %     thresholdImage=graythresh(afterProcessingImages(:,:,iframe)); %Otu thresholding images
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/(2^16-1));%1.2 for the image with few particles
    %     Max Entropy thresholding images
    %         afterProcessingImages(:,:,iframe)=im2bw_ent(afterProcessingImages(:,:,iframe));
    %     afterProcessingImages(:,:,iframe)=afterProcessingImages(:,:,iframe)&broderFrame; % remove pixels in the broder Frame
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe)); %remove the objectives in the border
    %     L=bwlabel(afterProcessingImages(:,:,iframe));
    %     postionsObjective=regionprops(L,beforeProcessingImages(:,:,iframe),'area','centroid','MaxIntensity');
    %     if ~isempty(postionsObjective), % empty stats causes crash
    %         afterProcessingImages(:,:,iframe)=removeFakeparticle(postionsObjective',afterProcessingImages(:,:,iframe));
    %     end
    %
    
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open');
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'bridge'); % bridge process
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image=afterProcessingImages(:,:,iframe);
            image(cc(iCC).PixelIdxList)=0;
            afterProcessingImages(:,:,iframe)=image;
        end
    end
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes'); % fill holes process
    %     cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MinIntensity','PixelIdxList');
    %     for iCC=1:size(cc,1)
    %         if cc(iCC).MinIntensity<=minIntensity
    %             image=afterProcessingImages(:,:,iframe);
    %             image(cc(iCC).PixelIdxList)=0;
    %             afterProcessingImages(:,:,iframe)=image;
    %         end
    %     end
    % afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    % afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % find outLine process
end
afterProcessingImages=logical(afterProcessingImages);
end


function [imageStack,bestPosition]= fluoImageCorrection(imageStack)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
imageSize=size(imageStack);
if min(imageSize(1:2))>=1500;
    imageStackNew=imageStack(200:1500,200:1500,:);
else
    imageStackNew=imageStack;
end
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);
end
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
% gfpImage=imageCorrectionWithBestPosition(gfpImage,bestPosition);
% rfpImage=imageCorrectionWithBestPosition(rfpImage,bestPosition);
end
function bestPosition=caculateCrossCorrelationForImage(image1,image2,step)
% for calculater the cross correlation of two image
% step is the searching range
% backGround should be calculater or set by user
[x,y]=meshgrid((-step:step)',(-step:step)');
x=x(:);
y=y(:);
correlationMatrix=zeros(size(x,1),1);
parfor i=1:numel(x)
    se=translate(strel(1),[x(i),y(i)]);
    image2New=imdilate(image2,se);  % 巧用imdilate实现平移
    sumImage=image1 & image2New;    % 利用逻辑矩阵的乘法相当于&
    correlationMatrix(i)=sum(sum(sumImage)); % 直接对逻辑矩阵求和速度比较快
end
bestPosition=[x(correlationMatrix==max(correlationMatrix)),y(correlationMatrix==max(correlationMatrix))];
bestPosition=bestPosition(1,:);
end
function image=imageCorrectionWithBestPosition(image,bestPosition)
% 已知漂移量后进行的较正
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
for i=2:size(bestPosition,1)
    se=translate(strel(1),bestPositionAccumulation(i,:));
    image(:,:,i)=imdilate(image(:,:,i),se);
end
end
