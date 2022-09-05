function [bioTree,dirFile]=batchTreeTrackingJob_withImageCorrelation(picType,dirFile)
global xyShift
xyShift=[0,0];
% dirFile=superBatchConvert(picType);
dirImage=strcat(dirFile,'\tiff2matlab');   %put the folder of the save results
dirResultSave=strcat(dirFile,'\bioTreeResult');   %put the folder of the save results
mkdir(dirFile,'bioTreeResult');
addpath(dirImage);
clc;
nameList=dir(dirImage);
for i=1:(length(nameList)-2)
    gainOneSmallTree(dirImage,i,dirResultSave,picType)
end
rmpath(dirImage)
rmpath(dirResultSave)
disp('combineSmallTreetoBigTree')
bioTree=combineSmallTreetoBigTree(i,dirResultSave);
end
function gainOneSmallTree(dirImage,i,dirResultSave,picType)
global xyShift
nameList=dir(dirImage);
jobNumber=strcat('start Job', num2str(i));
folderName=strcat('t',num2str(i));
mkdir(dirResultSave,folderName);
fileName=nameList(i+2).name;
disp(jobNumber);
disp(strcat('load file',fileName));
tic;disp(strcat('1.load:',fileName,'...'));
imageStack=myLoad(i,fileName,dirResultSave);toc;
frameShift=preStackSize(i, nameList(i+1).name);
clear tempImage;
disp('image Processing')
% imageStack=imageStack(647:647+999,1372:1372+999,:);
[maskImages,imageProcessingInfo]=myImageProcessing(imageStack,picType);
[maskImages,xyShift]=fluoImageCorrection(maskImages,xyShift);
[~,bioTree]=treeTrackingJob(imageStack,maskImages,1,size(maskImages,3),frameShift);
if i==1
    [~,backGroundPara]=backGroundCorrection([],[],'16bit');
    bioTree{1}.backGroundPara=backGroundPara;
    bioTree{1}.imageSize=[size(maskImages,1),size(maskImages,2)];
    bioTree{1}.imageProcessingInfo=imageProcessingInfo;
end
disp('save BioTree');
saveFile1=strcat(dirResultSave,'\',folderName,'\bioTree',num2str(i));
tic;bioTreeSave(bioTree,saveFile1);toc;
end
function imageStack=myLoad(iJob,fileName,dirResultSave)
if iJob==1
    tempImage=load(fileName);
    imageStack=tempImage.imageStack;
else
    tempImage=load(fileName);
    imageStack=tempImage.imageStack;
    addpath(dirResultSave);
    preLastImage=load('lastImage');
    imageStack=cat(3,preLastImage.lastImage,imageStack);
end
lastImage= imageStack(:,:,end);
saveFile=strcat(dirResultSave,'\lastImage');
save(saveFile,'lastImage');
end
function frameShift=preStackSize(iJob, fileName)
if iJob==1
    frameShift=0;
else
    vars = whos('-file', fileName);
    stackSize=vars.size(3);
    frameShift=(iJob-1)*stackSize-1;
end
end
function bioTreeSave(bioTree,saveFile)
bytesTemp=whos('bioTree');
bytes=bytesTemp.bytes;
treeNum=fix(bytes/2000000000)+1;
savePoint=findSavePoint(bioTree,treeNum);
startFrame=1;
for iTree=1:treeNum
    if treeNum== 1
        save(saveFile,'bioTree');
        return;
    end
    bioTreeStack=bioTree(startFrame:savePoint(iTree));
    saveFile1=strcat(saveFile,'-',num2str(treeNum));
    save(saveFile1,'bioTreeStack');
    startFrame=savePoint(iTree)+1;
    bioTreeStack=[];
end
end
function savePoint=findSavePoint(bioTree,treeNum)
savePoint=[];
if treeNum==1
    return;
end
startFrame=1;
for iframe=100:100:size(bioTree,2)
bioTreeTmep=bioTree(startFrame:iframe)  ;  
bytesTemp=whos('bioTreeTmep');
bytes=bytesTemp.bytes;
pointNum=fix(bytes/2000000000);
if pointNum==1
    savePoint=[savePoint,iframe-100];
    startFrame=iframe;
    if size(savePoint,2)==treeNum-1
        savePoint=[savePoint,size(bioTree,2)];
        return;
    end
end
end
end
function [imageStack,xyShift] = fluoImageCorrection(imageStack,xyShift)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
imageSize=size(imageStack);
if min(imageSize(1:2))>=1500;
    imageStackNew=imageStack(600:1100,600:1100,:);
else
    imageStackNew=imageStack;
end
image=imageStackNew(:,:,1);
cc=bwconncomp(image);
if cc.NumObjects<=1
    return
end
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);
end
[imageStack,xyShift]=imageCorrectionWithBestPosition(imageStack,bestPosition,xyShift);
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
function [image,xyShift]=imageCorrectionWithBestPosition(image,bestPosition,xyShift)
% 已知漂移量后进行的较正
bestPositionAccumulation=zeros(size(bestPosition));
bestPositionAccumulation(1,:)=xyShift;
for i=2:size(bestPosition,1)
    bestPositionAccumulation(i,:)=bestPosition(i,:)+bestPositionAccumulation(i-1,:);
end
for i=1:size(bestPosition,1)
    se=translate(strel(1),bestPositionAccumulation(i,:));
    image(:,:,i)=imdilate(image(:,:,i),se);
end
xyShift=bestPositionAccumulation(end,:);
end