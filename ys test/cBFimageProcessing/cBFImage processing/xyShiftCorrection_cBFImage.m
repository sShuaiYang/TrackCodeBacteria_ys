function [imageStack,bestPosition] = xyShiftCorrection_cBFImage(imageStack)
%bestPosition=stepCorrectionForImageStack(imageStack)
%   根据jzy上述函数修改而来 2020.01.02ys
%通过原始的相关明场图像进行相关计算xy漂移并进行矫正

% 计算出imageStack的偏移量，并对imageStack进行校正
tic;
imageSize=size(imageStack);
if min(imageSize(1:2))>=1000
    imageStackNew=imageStack(800:1100,800:1100,:);
else
    imageStackNew=imageStack;
end

bestPosition(1,:)=[size(imageStackNew,1),size(imageStackNew,2)];
parfor i=1:size(imageStackNew,3)-1
        correlationMatrix=caculateCrossCorrelationForBFImage(imageStackNew(:,:,i),imageStackNew(:,:,i+1),10);
    %     correlationMatrix=xcorr2(imageStack(:,:,i)-100,imageStack(:,:,i+1)-100);
%     correlationMatrix=normxcorr2(imageStackNew(:,:,i),imageStackNew(:,:,i+1));%normxcorr2
%     normxcorr2比xcorr2运算速快2-3倍，但是效果不好。normxcorr2对象对其，参考matlab help
    [a,b]=find(correlationMatrix==max(max(correlationMatrix)));
    bestPosition(i+1,:)=[a(1),b(1)];
end
bestPosition=bestPosition-bestPosition(1,:);
% imageStack1=imageStackCorrectionWithBestPosition(imageStack,bestPosition);
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
toc;
end
%%  图像相关计算
function correlationMatrix=caculateCrossCorrelationForBFImage(image1,image2,step)
% for calculater the cross correlation of two image
% step is the searching range
%来自jzy函数correlationMatrix=caculateCrossCorrelationForImage(image1,image2,step,backGround)
backGround=0;
fullSize=[size(image1,1)+size(image2,1)-1,size(image1,2)+size(image2,2)-1];
correlationMatrix=zeros(fullSize);
centroid=floor((fullSize+1)/2);
matrixSequence=(centroid(1)-step:centroid(2)+step)';
[x,y]=meshgrid(matrixSequence,matrixSequence);
x=x(:);
y=y(:);
delta1=floor((size(image2,1)+1)/2)-1;
delta2=size(image2,1)-floor((size(image2,1)+1)/2);
for i=1:numel(x)
    fullImage=double(ones(fullSize)*backGround);
    fullImage(x(i)-delta1:x(i)+delta2,y(i)-delta1:y(i)+delta2)=image2;
    image2New=fullImage(centroid(1)-delta1:centroid(1)+delta2,centroid(2)-delta1:centroid(2)+delta2);
    correlationMatrix(x(i),y(i))=sum(sum(double(image1).*image2New));
%     correlationMatrix(i)=sum(sum(double(image1).*image2New));
end

end
%%
function  imageStack1=imageStackCorrectionWithBestPosition(imageStack,bestPosition)
% old code 
backGround=0;
deltaPosition=bestPosition;
imageSize=[size(imageStack,1),size(imageStack,2)];
for i=1:size(deltaPosition,1)-1
    deltaPosition(i+1,:)=deltaPosition(i+1,:)+deltaPosition(i,:);
end
bestPosition1=deltaPosition+imageSize;
fullSize=[size(imageStack,1)+size(imageStack,1)-1,size(imageStack,2)+size(imageStack,2)-1];
centroid=floor((fullSize+1)/2);
delta1=floor((size(imageStack,1)+1)/2)-1;
delta2=size(imageStack,1)-floor((size(imageStack,1)+1)/2);
parfor i=2:size(imageStack,3)
    fullImage=uint16(ones(fullSize)*backGround);
    fullImage(bestPosition1(i,1)-delta1:bestPosition1(i,1)+delta2,bestPosition1(i,2)-delta1:bestPosition1(i,2)+delta2)=imageStack(:,:,i);
    imageStack1(:,:,i)=fullImage(centroid(1)-delta1:centroid(1)+delta2,centroid(2)-delta1:centroid(2)+delta2);
end
end
%% 
function image=imageCorrectionWithBestPosition(image,bestPosition)
% 已知漂移量后进行的较正
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
parfor i=2:size(bestPosition,1)
    se=translate(strel(1),bestPositionAccumulation(i,:));
    image(:,:,i)=imdilate(image(:,:,i),se);
end
end

