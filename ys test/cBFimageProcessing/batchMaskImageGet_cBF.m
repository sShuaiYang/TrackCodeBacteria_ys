function batchMaskImageGet_cBF()
dirFile='F:\2020-04-27 PAO1_IP31_100x 1agar Temp30 ys';
fieldList=dir(dirFile);
for iField=1:64%length(nameList)-2
    t0=clock;
    
    disp('start image processing and xy shift correction ')
    dirField=[dirFile,'\',fieldList(iField+2).name];
    dirTrack=[dirField,'\Tracking'];
    mkdir(dirTrack);
    dirMaskCheck= [dirField,'\maskCheck'];
    mkdir(dirMaskCheck);
    dirMaskSave=strcat(dirField,'\tiff2matlab');
    mkdir(dirMaskSave);
    dirImage=[dirField,'\','BF'];
    imageList=dir(dirImage);
    %去除Thumbs.db文件影响文件读取
    if isfile(strcat(dirImage,'\Thumbs.db'))
        imageList = thumbsFileRemove(imageList);
    end
      
    disp(fieldList(iField+2).name);
  
    %     imageNum=length(imageList)-3;
    imageNum=90;
    checkGap=5;%每隔5frame将图片储存，进行图像处理和校正的check%建议每5min存一张
    
    if imageNum<=400
        disp('imageStack 0001')
        tic;
        tempImages = zeros(2048,2048,imageNum,'uint8');% cBF image type uint8
        parfor i=1:imageNum %i=1:length(imageList)-3
            %         imageStack(:,:,i)=import_tiff_stack( strcat(dirImage,'\',imageList(i+3).name) );
            tempImages(:,:,i) = import_tiff_stack( strcat(dirImage,'\',imageList(i+3).name) );            
        end
        [imageTrackings,~] = myImageProcessing_cBFys(tempImages,'w');
        %xy shift correction
        [imageTrackings,~] = fluoImageCorrection(imageTrackings);
        % save maskImages
        for iframe=1:imageNum
            filename = [dirTrack,'\','imageTracking',imageList(iframe+3).name(end-8:end-4),'.mat'];
            imageTracking = imageTrackings(:,:,iframe);
            save(filename,'imageTracking');            
        end
        j = 1;
        maskImages = false(size(imageTrackings,1),size(imageTrackings,2),size(imageTrackings,3)*2);
        maskImages(:,:,1:2:end) = imageTrackings;
        maskImages(:,:,2:2:end) = imageTrackings;
        maskImages=logical(maskImages); 
%         bigImageStackSplitSave(dirMaskSave,maskImages,200);
        save([dirMaskSave,'\',num2str(j,'%05.f')],'maskImages','-v7.3')   
        
        maskCheckStack=imageTrackings(:,:,1:checkGap:end); 
        save([dirMaskCheck,'\','maskCheckStack'],'maskCheckStack');
        clear tempImages imageTrackings imageTracking maskImages 
        toc;
    else
        tic;
        n=0;
        j=0;
        stackGap = 200;% 保存的frame间隔,每200frame进行处理
        maskCheckStack=[];
        checkGap = 5;%每隔5frame将图片储存，进行图像处理和校正的check
        for i=1:imageNum
            n=n+1;            
            tempImages(:,:,n)=import_tiff_stack( strcat(dirImage,'\',imageList(i+3).name) );
            if n==stackGap||(j==floor(imageNum/stackGap)&&n==mod(imageNum,stackGap))
                disp (['imageStack',num2str((j+1),'%05.f')]);
                
                [imageTrackings,~]=myImageProcessing_cBFys(tempImages,'w');
                %xy shift correction
                if j==0
                    [imageTrackings,~]=fluoImageCorrection(imageTrackings);
                else
                    imageTrackings(:,:,2:end+1)=imageTrackings;
                    imageTrackings(:,:,1)=labelImage;
                    [imageTrackings,~]=fluoImageCorrection(imageTrackings);
                    imageTrackings=imageTrackings(:,:,2:end);
                end
                labelImage=imageTrackings(:,:,end);
                % save maskImages 
                for iframe=1:n
                    filename=[dirTrack,'\','imageTracking',imageList(j*stackGap+iframe+3).name(end-8:end-4),'.mat'];
                    imageTracking=imageTrackings(:,:,iframe);
                    save(filename,'imageTracking');
                end
                j=j+1;
                
                maskImages=false(size(imageTrackings,1),size(imageTrackings,2),size(imageTrackings,3)*2);
                maskImages(:,:,1:2:end)=imageTrackings;
                maskImages(:,:,2:2:end)=imageTrackings;
                maskImages=logical(maskImages); 
                save([dirMaskSave,'\',num2str(j,'%05.f')],'maskImages');  
               
                n=0;
                maskCheckStack=cat(3,maskCheckStack,imageTrackings(:,:,1:checkGap:end));  
                clear imageTrackings imageTracking tempImages maskImages
            end
        end
        maskCheckStack=logical(maskCheckStack);        
        %         bigImageStackSplitSave(dirMaskCheck,maskCheckStack);
        save([dirMaskCheck,'\','maskCheckStack'],'maskCheckStack');
        
        toc;
    end   
    
    load([dirImage,'\frameInfo.mat'])
    frameInfo=frameInfo(1:imageNum,:);
    save([dirTrack,'\','frameInfo.mat'],'frameInfo'); 
    
    frameInfoCheck=frameInfo(1:checkGap:imageNum,:);    
    cellNum=zeros(2,size(maskCheckStack,3));
    for i=1:size(maskCheckStack,3)
        stats = regionprops(maskCheckStack(:,:,i),'Centroid');
        cellNum(1,i)=numel(stats);
        cellNum(2,i)=etime(frameInfoCheck(i,1:6),frameInfoCheck(1,1:6))/60;   %min
    end
    figure,plot(cellNum(2,:),cellNum(1,:))
    saveas(gcf,[dirMaskCheck,'\','cellNum','.fig']);
    saveas(gcf,[dirMaskCheck,'\','cellNum','.tif']);
    close all;
    clear maskCheckStack
    
    disp(['耗时:',num2str(etime(clock,t0)),'秒']);
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
function  bigImageStackSplitSave(dirSave,imageStack,stackGap)
   % save maskImages
    n=0;
    j=0;
%     stackGap=200;
    for i=1:size(imageStack,3)
        n=n+1;
        smallStack(:,:,n)=imageStack(:,:,i);
        if n==stackGap||(j==floor(size(imageStack,3)/stackGap)&&n==mod(size(imageStack,3),stackGap))
            j=j+1;
            maskImages=smallStack;%重命名文件
            save([dirSave,'\',num2str(j,'%05.f')],'maskImages');
            n=0;
            clear smallStack maskImages
        end
    end
end
%%
function  fileList = thumbsFileRemove(fileList)
templogic = true (1, numel (fileList));
for iFile = 1: numel (fileList)
    
    if strcmp (fileList(iFile).name, 'Thumbs.db')
        templogic(iFile) = false;
        fileList = fileList(templogic);
        return
    end
    
end
end