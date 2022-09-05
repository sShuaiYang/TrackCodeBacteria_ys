function  stepByStepMaskImageGet()
%һ������ һ�߶Ի�õ����ݽ���ͼ��ʶ��imageProcessing��xyshifit����
%Ŀǰ��Ҫ�������Ķ���Ұ�������ͼ��
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
imageAmount = 10;  % ÿһ�δ����image��Ŀ stack��ʽ����
loopIndex = 0;

labelImage=cell(1,length(fieldList)-3);

while loopIndex <= 2 % loopIndex*imageAmount �����ܹ���Ҫ����image����Ŀ
    
    prcessedImIdx = loopIndex*imageAmount;% �Ѿ������image����Ŀ processed image index
    for ifield = 1:length(fieldList)-3
        disp(['Round ',num2str(loopIndex),' of ', fieldList(ifield+2).name]);
        
        dirFile = [dirAll,'\',fieldList(ifield+2).name];
        % ��ȡ��������imageAmount��image ������� �ȴ�
        waitingTime = 0;
%         disp('waiting')
        dirOriginalImage=[dirFile,'\','BF'];
        

        % ��ϵͳһֱ�ȴ�������30���ӣ�
        while waitingTime <= 30*60
            tic;imageList = dir(dirOriginalImage);
            newImageNum = numel(imageList)-3-prcessedImIdx;
            wait1 = toc;
            if newImageNum < imageAmount
                pause(60)
                waitingTime = wait1+60+waitingTime;    % ���û���µ�Image,�ȴ�60s����������һ��ѭ��
%                 clc
                disp(['waitingTime=',num2str(waitingTime)]);
            else
                if newImageNum >= imageAmount      % ���µ���Ƭ��Ŀ�������涨��Ŀ����ǿ��ת��Ϊ�������涨��Ŀ��ʣ�µ��´��ټ���
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
                
                % ������Ƭ
                for iframe=1:newImageNum
                    filename=[dirFile,'\Tracking','\','imageTracking',imageList(loopIndex*imageAmount+iframe+3).name(end-8:end-4),'.mat'];
                    imageTracking=imageTrackings(:,:,iframe);
                    save(filename,'imageTracking');
                end
                
                break;%�˳�ѭ��ִ����һ��field
                
            end
        end
        
    end
    
    loopIndex=loopIndex+1;
    
end

end
%% xyshift correction
function [imageStack,bestPositionAccumulation] = fluoImageCorrection(imageStack)
% ���ӫ��ͼ���ʱ�����е�xy����
% �������Ϊ��ֵͼ�񣬺��Ľ��������Ѿ��Ż���
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
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);%ƫ�����[-15,15];ys2020.01.03
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
    image2New=imdilate(image2,se);  % ����imdilateʵ��ƽ��
    sumImage=image1 & image2New;    % �����߼�����ĳ˷��൱��&
    correlationMatrix(i)=sum(sum(sumImage)); % ֱ�Ӷ��߼���������ٶȱȽϿ�
end
bestPosition=[x(correlationMatrix==max(correlationMatrix)),y(correlationMatrix==max(correlationMatrix))];
bestPosition=bestPosition(1,:);
end
function image=imageCorrectionWithBestPosition(image,bestPosition)
% ��֪Ư��������еĽ���
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
%imageStackΪ��ֵͼ�� 2020.01.03ys Ѱ��ͼ���и�ķ�Χ
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
