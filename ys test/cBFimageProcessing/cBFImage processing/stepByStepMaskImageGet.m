function  stepByStepMaskImageGet()
%一边拍摄 一边对获得的数据进行图像识别imageProcessing和xyshifit矫正
%目前主要针对拍摄的多视野相关明场图像
% Shuai Yang 2020.03.27
dirAll='D:\2020-03-18 PAO1_IP32_100x ys_steptest';
% dirAll = uigetdir();
fieldList = dir(dirAll);
disp(['filed number = ',num2str(length(fieldList)-3)]);

for ifield = 1:length(fieldList)-3
    
    dirFile = [dirAll,'\',fieldList(ifield+2).name];
    dirTrack = [dirFile,'\Tracking'];
    mkdir(dirTrack);
%     dirMaskCheck = [dirFile,'\maskCheck'];
%     mkdir(dirMaskCheck);
%     dirMaskSave = strcat(dirFile,'\tiff2matlab');
%     mkdir(dirMaskSave);
    
end

disp('Ready for image processing and xy shift correction ')
imageAmount = 10;  % 每一次处理的image数目 stack方式处理
loopIndex = 0;

labelImage=cell(1,length(fieldList)-3);

while loopIndex <= 2 % loopIndex*imageAmount 决定总共需要处理image的数目
    
    prcessedImIdx = loopIndex*imageAmount;% 已经处理的image的数目 processed image index
    for ifield = 1:length(fieldList)-3
        disp(['Round ',num2str(loopIndex),' of ', fieldList(ifield+2).name]);
        
        dirFile = [dirAll,'\',fieldList(ifield+2).name];
        % 读取数量等于imageAmount的image 如果不够 等待
        waitingTime = 0;
%         disp('waiting')
        dirOriginalImage=[dirFile,'\','BF'];
        

        % 让系统一直等待，最多等30分钟，
        while waitingTime <= 30*60
            tic;imageList = dir(dirOriginalImage);
            newImageNum = numel(imageList)-3-prcessedImIdx;
            wait1 = toc;
            if newImageNum < imageAmount
                pause(60)
                waitingTime = wait1+60+waitingTime;    % 如果没有新的Image,等待60s，并进行下一次循环
%                 clc
                disp(['waitingTime=',num2str(waitingTime)]);
            else
                if newImageNum >= imageAmount      % 若新的照片数目大于最大规定数目，则强制转化为等于最大规定数目，剩下的下次再计算
                    newImageNum = imageAmount;
                end
                
                imageStack = zeros(2048,2048,newImageNum);
                parfor iImage = 1:newImageNum
                    imageStack(:,:,iImage) = import_tiff_stack([dirOriginalImage,'\',imageList(prcessedImIdx+iImage+3).name]);
                end
                
                [imageTrackings,~] = myImageProcessing_cBFys(imageStack,'w');
                %xy shift correction
                if loopIndex==0
                    [imageTrackings,~] = fluoImageCorrection(imageTrackings);
                else
                    imageTrackings(:,:,2:end+1) = imageTrackings;
                    imageTrackings(:,:,1) = labelImage{ifield};
                    [imageTrackings,~] = fluoImageCorrection(imageTrackings);
                    imageTrackings = imageTrackings(:,:,2:end);
                end
                labelImage{ifield} = imageTrackings(:,:,end);
                
                % 保存照片
                for iframe=1:newImageNum
                    filename=[dirFile,'\Tracking','\','imageTracking',imageList(loopIndex*imageAmount+iframe+3).name(end-8:end-4),'.mat'];
                    imageTracking=imageTrackings(:,:,iframe);
                    save(filename,'imageTracking');
                end
                
                break;%退出循环执行下一个field
                
            end
        end
        
    end
    
    loopIndex=loopIndex+1;
    
end

end
%% xyshift correction
function [imageStack,bestPositionAccumulation] = fluoImageCorrection(imageStack)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
% imageSize=size(imageStack);
% if min(imageSize(1:2))>=1500
%     imageStackNew=imageStack(1024-500:1024+500,1024-500:1024+500,:);
% %     imageStackNew=imageStack(600:1100,600:1100,:);
% else
%     imageStackNew=imageStack;
% end
[imageStackNew]=cropImageForxyShiftCorrection(imageStack,400);% 2020,01.03ys
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);%偏移最大[-15,15];ys2020.01.03
end
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
bestPositionAccumulation=zeros(size(bestPosition));
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
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
%% 
