function multiFieldFluoImageBioTreeGeneration(dirFile)
disp('MultiFiled fluo image bioTree generation');
fieldList =  dir([dirFile,filesep,'field*']);
% cd(dirFile);
for iField= 1:length(fieldList)
    t0=clock;
    if ~strcmp(fieldList(iField).name(1:5),'field')
        continue
    end
    disp(fieldList(iField).name);
    disp(compose("Start time %d-%d-%d %.2d:%.2d:%02.0f.",t0));
    
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    [bioTree] = batchTreeTrackingJob(dirField);
    %     bioTree = bioTreeMeasure(bioTree,0,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
    bioTree = bioTreeDigHoleAndMeasure(bioTree,1,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
    bioTree = bioTreeTwoPointTracking(bioTree);
    %     [bioTree,branchList,~,~] = myBiograph_new2(bioTree);
    %     [bioTree,~,~] = divisionFinder(bioTree,branchList);7
    
    dirFinalbioTreeSave=strcat(dirField,'\bioTreeResult\finalTree');
    mkdir(dirFinalbioTreeSave)
    save (strcat(dirFinalbioTreeSave,'\bioTree.mat'),'bioTree','-v7.3');
    clear bioTree
    t1 = clock;
    disp(compose("End time %d-%d-%d %.2d:%.2d:%02.0f.",t1))
    disp(['共耗时:',num2str(etime(t1,t0)),'秒']);
end
end
%%
function bioTree = batchTreeTrackingJob(dirFile)
dirImage=strcat(dirFile,'\tiff2matlab');   %put the folder of the save results
dirResultSave=strcat(dirFile,'\bioTreeResult');   %put the folder of the save results
mkdir(dirFile,'bioTreeResult');
addpath(dirImage);
addpath(dirResultSave);
% clc;
nameList=dir(dirImage);
for i=1:(length(nameList)-2)
    if i==1
        imageSize=gainOneSmallTree(dirImage,i,dirResultSave);
    else
        imageSize=gainOneSmallTree(dirImage,i,dirResultSave,imageSize);
    end
end
rmpath(dirImage)
rmpath(dirResultSave)

disp('combineSmallTreetoBigTree')
bioTree=combineSmallTreetoBigTree(length(nameList)-2,dirResultSave);
disp('Lucky, succeed in combining small tree')

% disp('Now we begin bioTreeMeasure')
bioTree=bioTreeMeasure(bioTree,0,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
% mkdir(dirResultSave,'\allTreeAfterProcesing');
mkdir(dirResultSave,'allTree');
allTreeSave(bioTree,strcat(dirResultSave,'\allTree'),'bioTreeStack');
end
%%
function imageSize = gainOneSmallTree(dirImage,i,dirResultSave,imageSize)
nameList=dir(dirImage);
jobNumber=strcat('start Job', num2str(i));
folderName=strcat('t',num2str(i));
mkdir(dirResultSave,folderName);
fileName=nameList(i+2).name;
disp(jobNumber);
disp(strcat('load file',fileName));
tic;disp(strcat('1.load:',fileName,'...'));
maskImages=myLoad(i,fileName,dirResultSave);toc;
frameShift=preStackSize(i, nameList(i+1).name);
clear tempImage;
disp('image Processing')
% maskImages=myImageProcessingBlack(imageStack,0);
% maskImages=myImageProcessing(imageStack,0);

[~,bioTree]=treeTrackingJob(maskImages,maskImages,1,size(maskImages,3),frameShift);
if nargin==3
    bioTree{1}.imageSize=[size(maskImages,1),size(maskImages,2)];
    imageSize=bioTree{1}.imageSize;
end
disp('save BioTree');
saveFile1=strcat(dirResultSave,'\',folderName,'\bioTree',num2str(i));
bytesTemp=whos('bioTree');
bytes=bytesTemp.bytes;
if bytes>2000000000
    tic;save(saveFile1,'bioTree','-v7.3');toc;
else
    tic;save(saveFile1,'bioTree');toc;
end
clear maskImages
end
function maskImages=myLoad(iJob,fileName,dirResultSave)
if iJob==1
    tempImage=load(fileName);
    maskImages=tempImage.imageTrackings;
else
    tempImage=load(fileName);
    maskImages=tempImage.imageTrackings;
    addpath(dirResultSave);
    preLastImage=load('lastImage');
    maskImages=cat(3,preLastImage.lastImage,maskImages);
end
lastImage=maskImages(:,:,end);
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

function treeNum=allTreeSave(bioTree,saveFile,name)
saveFile=strcat(saveFile,'\',name);
bytesTemp=whos('bioTree');
bioTree1=bioTree(1);
bytesTemp1=whos('bioTree1');
bytes=bytesTemp.bytes-bytesTemp1.bytes;
treeNum=fix(bytes/1900000000)+2;
savePoint=1;
savePoint=findSavePoint(bioTree,treeNum,savePoint);
startFrame=1;
for iTree=1:treeNum
    if treeNum==1
        save(saveFile,'bioTree');
        return;
    end
    bioTreeStack=bioTree(startFrame:savePoint(iTree));
    saveFile1=strcat(saveFile,num2str(iTree));
    if savePoint(iTree)==1
        save(saveFile1,'bioTreeStack','-v7.3');
    else
        save(saveFile1,'bioTreeStack');
    end
    startFrame=savePoint(iTree)+1;
    bioTreeStack=[];
end
end
function savePoint=findSavePoint(bioTree,treeNum,savePoint)
if treeNum==1
    return;
end
if ~isempty(savePoint)
    startFrame=2;
else
    startFrame=1;
end
if size(savePoint,2)==treeNum-1
    savePoint=[savePoint,size(bioTree,2)];
    return;
end
for i=1:size(bioTree,2)
    bioTreeCell=bioTree(i);
    bytesTemp=whos('bioTreeCell');
    bytesCell(i)=bytesTemp.bytes;
end
if mod(size(bioTree,2),100)~=0
    iframeExtra=size(bioTree,2);
else
    iframeExtra=[];
end
for iframe=[100:100:size(bioTree,2),iframeExtra]
    bioTreeTmep=bioTree(startFrame:iframe);
    bytesTemp=whos('bioTreeTmep');
    bytes=bytesTemp.bytes;
    pointNum=fix(bytes/1900000000);
    if pointNum==1
        smallbytes=bytesCell(startFrame:iframe);
        bytesNum=1:numel(smallbytes);
        for i=2:numel(smallbytes)
            smallbytes(i)=smallbytes(i)+smallbytes(i-1);
        end
        smallbytes=smallbytes/1000/1000/1000;
        bytesNum=max(bytesNum(smallbytes<=1.9));
        savePoint=[savePoint,startFrame+bytesNum-1];
        startFrame=startFrame+bytesNum;
        if size(savePoint,2)==treeNum-1
            savePoint=[savePoint,size(bioTree,2)];
            return;
        end
    end
end
end