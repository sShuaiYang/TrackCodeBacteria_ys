function afterAddNew=multiFieldXyCorrection(beforeAdd,afterAdd,backGround)
% 专门为H2O2杀菌HJD1，xyScan多视野扫描而写的程序
% 默认channel为四个，第一个是BF,GFP,RFP,PI
% 由于xyScan的范围比较大，台子会发生一定的偏移，故使用该程序进行校正
% 需要backGRound是为了校正明场，再进行一些简易的识别轮廓
channelNum=4;
maskImage_before=getMaskImageFromBrightField(beforeAdd(:,:,1:channelNum:end),backGround);
maskImage_after=getMaskImageFromBrightField(afterAdd(:,:,1:channelNum:end),backGround);
for i=1:size(maskImage_before,3)
    correlationMatrix=caculateCrossCorrelationForImage(maskImage_before(:,:,i),maskImage_after(:,:,i),20);
    [a,b]=find(correlationMatrix==max(max(correlationMatrix)));
    afterAdd(:,:,4*i-3)=imageCorrectionWithBestPosition(afterAdd(:,:,4*i-3),[a,b]);
    afterAdd(:,:,4*i-2)=imageCorrectionWithBestPosition(afterAdd(:,:,4*i-2),[a,b]);
    afterAdd(:,:,4*i-1)=imageCorrectionWithBestPosition(afterAdd(:,:,4*i-1),[a,b]);
    afterAdd(:,:,4*i)=imageCorrectionWithBestPosition(afterAdd(:,:,4*i),[a,b]);
end
afterAddNew=afterAdd;
end
function maskImage=getMaskImageFromBrightField(bfImage,backGround)
bfImage=uint16(double(bfImage)/max(max(max(double(bfImage))))*1000);
backGround=uint16(double(backGround)/max(max(double(backGround)))*1000);
bfImage=backGroundCorrection(bfImage,backGround);
bfImage=255-bfImage;
gaussianFilter=fspecial('gaussian',[5, 5],20);
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
for iframe=1:size(bfImage,3)
bfImage(:,:,iframe)=imfilter(bfImage(:,:,iframe),gaussianFilter);
bfImage(:,:,iframe)=imfilter(bfImage(:,:,iframe),edgeFilter);
bfImage(:,:,iframe)=im2bw(bfImage(:,:,iframe),210/255);
bfImage(:,:,iframe)=imclearborder(bfImage(:,:,iframe));
end
maskImage=logical(bfImage);
end
function [imageStack,backGroundPara]=backGroundCorrection(imageStack,backGround)
costantHighIn=1.6; %16 bit image use 1.25; 8bit image use 1.5
costanLowIn=0.85; %16 bit image use 0.85; 8bit image use 0.6
ampIntensity=3.0; %16 bit image use 2.8; 8bit image use 1
c=0.6; %16 bit image use 0.6; 8bit image use 0.3
fprintf('\n');
backImage=double(backGround);
averageBack=mean(mean(backGround));
for iframe=1:size(imageStack,3)
    imageTemp=double(imageStack(:,:,iframe));
    averageI=mean(mean(imageTemp));
    ratio= averageI/averageBack;
    highIn=averageBack*costantHighIn;
    lowIn=averageBack*costanLowIn;
    slope=1/(highIn-lowIn);
    b=-lowIn/(highIn-lowIn);
    imageStack(:,:,iframe)=((((imageTemp./ratio)./backImage).*averageBack).*slope+b).*255; % backGround Normilzed
    imageStack(:,:,iframe)=(imageStack(:,:,iframe)-c*mean(mean( imageStack(:,:,iframe)))).*ampIntensity; % Intensity AMP
end
fprintf('\n');
imageStack=uint8(imageStack);
backGroundPara.costantHighIn=costantHighIn;
backGroundPara.costanLowIn=costanLowIn;
backGroundPara.ampIntensity=ampIntensity;
backGroundPara.c=c;
end
function image=imageCorrectionWithBestPosition(image,bestPosition)
imageSize=size(image);
deltaPosition=bestPosition-imageSize;
for i=1:size(deltaPosition,1)-1
    deltaPosition(i+1,:)=deltaPosition(i+1,:)+deltaPosition(i,:);
end
bestPosition1=deltaPosition+imageSize;
fullSize=[size(image,1)+size(image,1)-1,size(image,2)+size(image,2)-1];
centroid=floor((fullSize+1)/2);
delta1=floor((size(image,1)+1)/2)-1;
delta2=size(image,1)-floor((size(image,1)+1)/2);
fullImage=uint16(ones(fullSize)*100);
fullImage(bestPosition1(1,1)-delta1:bestPosition1(1,1)+delta2,bestPosition1(1,2)-delta1:bestPosition1(1,2)+delta2)=image;
image=fullImage(centroid(1)-delta1:centroid(1)+delta2,centroid(2)-delta1:centroid(2)+delta2);
end