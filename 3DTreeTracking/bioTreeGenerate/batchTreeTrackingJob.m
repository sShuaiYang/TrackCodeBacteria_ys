function [bioTree,dirFile]=batchTreeTrackingJob(picType,dirFile)
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
function gainOneSmallTree(dirImage,i,dirResultSave,~)
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
maskImages=imageStack;
imageProcessingInfo.cropInfo=true(512,512);
% [maskImages,imageProcessingInfo]=myImageProcessing(imageStack,picType);
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
function bioTree=nodeReductionforBioTree(bioTree)
[~,nodeList]=findNode(bioTree);
nodeListNew=[];
i=0;
while ~isequal(nodeList,nodeListNew)
    i=i+1;
    if i~=1
        nodeList=nodeListNew;
    end
    disp(strcat('round_',num2str(i)))
    disp('type1__')
    tic;bioTree=type1NodeReduction(bioTree);toc
    disp('type_In__')
    tic;bioTree=type_InNodeReduction(bioTree);toc;
    disp('type_Out__')
    tic;bioTree=type_OutNodeReduction(bioTree);toc;
    disp('type1Deep__')
    tic;bioTree=type1NodeReduction_deep(bioTree);toc;
    disp('type2__')
    tic;bioTree=type2NodeReduction(bioTree);toc;
    disp('type3__')
    tic;bioTree=type3NodeReduction(bioTree);toc;
    disp('type4__')
    tic;bioTree=type4NodeReduction(bioTree);toc;
    [~,nodeListNew]=findNode(bioTree);
end
end
function allTreeSave(bioTree,saveFile)
saveFile=strcat(saveFile,'\bioTeeStack');
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
    saveFile1=strcat(saveFile,num2str(iTree));
    save(saveFile1,'bioTreeStack');
    startFrame=savePoint(iTree)+1;
    bioTreeStack=[];
end
end