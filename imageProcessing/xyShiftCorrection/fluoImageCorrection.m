function [imageStack,gfpImage] = fluoImageCorrection(imageStack,gfpImage)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
tic
% imageSize=size(imageStack);
% if min(imageSize(1:2))>=1500
%     imageStackNew=imageStack(700:1400,700:1400,:);
% else
%     imageStackNew=imageStack;
% end
[imageStackNew]=cropImageForxyShiftCorrection(imageStack,250);
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),10);
end
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
% gfpImage=imageCorrectionWithBestPosition(gfpImage,bestPosition);
% rfpImage=imageCorrectionWithBestPosition(rfpImage,bestPosition);
toc
end
%% 
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
%% 
function [imageStackNew]=cropImageForxyShiftCorrection(imageStack,pixelRange)
%imageStack为二值图像 2020.01.03ys 寻找图像切割的范围
%   pixelRange=250;
imageSize=size(imageStack);
if min(imageSize(1:2))>=1000
    xCrop=[fix(imageSize(1)/2)-pixelRange,fix(imageSize(1)/2)+pixelRange];
    yCrop=[fix(imageSize(2)/2)-pixelRange,fix(imageSize(2)/2)+pixelRange];
    tempStack=imageStack(xCrop(1):xCrop(2),yCrop(1):yCrop(2),:);
    stats = regionprops(tempStack(:,:,1),'PixelList');
    if numel(stats)<7
        pixelRange=pixelRange+200;
        if pixelRange>=min(imageSize(1:2))/2
            imageStackNew=imageStack;
        else
            [imageStackNew]=cropImageForxyShiftCorrection(imageStack,pixelRange);
        end
    else
        imageStackNew=tempStack;        
    end
    
else
    imageStackNew=imageStack;
end

end

