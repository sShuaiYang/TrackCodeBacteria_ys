function bioTree=longTimeProtocal(picType,dirFile)
for i=1
    % for i=1:size(nameList,1)
    clc
    %     dirFile=nameList{i,1};
    %         disp(nameList{i})
    if isempty(dirFile)
        dirFile=uigetdir();
    end
    dataOut=scaleInforead(strcat(dirFile,'\information\imageInfo.txt'));
    [bioTree,dirFile]=generateLongTimeBioTreeAndMakeMovie(picType,dirFile,i,dataOut);
        longTimeStatistic(bioTree,dirFile);
    %     allData=gainAllData(bioTree);
    %     allData=detailsForBacteriaGrowthandMoving(allData,50,strcat(dirFile,'\longTimeResult'));
    %     save(strcat(dirFile,'\longTimeResult\allData'),'allData');
    step=200;
    clusterStatistic(bioTree,dirFile,step,bioTree{1}.scaleInfo);
    %     [~,correctRateAll]=bioTreeCorrectRate(bioTree);
    %     disp(correctRateAll)
    % %     clear bioTree
    %     clear allData
%     close all
end
end
function [bioTree,dirFile]=generateLongTimeBioTreeAndMakeMovie(picType,dirFile,iStack,scaleInfo)
[bioTree,dirFile]=batchTreeTrackingJob(picType,dirFile);
bioTree=autoTidyBioTree(bioTree);
nodeSortPre=[];
[nodeSort,~]=findNode(bioTree);
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,1,25);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
dirForSeg=strcat(dirFile,'\bioTreeResult\allTree\forSeg');
mkdir(dirForSeg);
allTreeSave(bioTree,dirForSeg,'forSeg');
mkdir(strcat(dirFile,'\bioTreeResult\allTree'))
cd(strcat(dirFile,'\bioTreeResult\allTree'))
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type_InNodeReduction(bioTree,1);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,2,25);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,3,40);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,1);
    bioTree=mixMatchReducion(bioTree,5);
    bioTree=type3NodeReduction(bioTree,15);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,3,40);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=mixMatchReducion(bioTree,10);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
    bioTree=autoTidyBioTree(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
bioTree=mergeNode(bioTree);
nodeSortPre=[];
while ~isequal(nodeSortPre,nodeSort)
    nodeSortPre=nodeSort;
    bioTree=type1NodeReduction(bioTree,3,40);
    bioTree=type_InNodeReduction(bioTree,2);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=mixMatchReducion(bioTree,10);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
    [nodeSort,~]=findNode(bioTree);
end
bioTree=leafRootRefineBioTree(bioTree,0);
bioTree=bioTreeOptimize(bioTree);
bioTree{1}.scaleInfo=scaleInfo;
dirAfterSeg=strcat(dirFile,'\bioTreeResult\allTree\afterSeg');
mkdir(dirAfterSeg);
treeNum=allTreeSave(bioTree,dirAfterSeg,'afterSeg');
clear bioTree
dirFinalTree=strcat(dirFile,'\bioTreeResult\allTree\finalTree');
mkdir(dirFinalTree)
for i=1:treeNum
    load(strcat(dirAfterSeg,'\afterSeg',num2str(i)));
    if i==1
        imageSize=bioTreeStack{1}.imageSize;
    end
    bioTreeStack=bioTreeDigHoleAndMeasure(bioTreeStack,1,imageSize(1),imageSize(2));
    allFinalTreeSave(bioTreeStack,dirFinalTree,strcat('_',num2str(i),'_'));
end
measureResult=dir(dirFinalTree);
measureResultNum=numel(measureResult)-2;
for i=1:measureResultNum
    load(strcat(dirFinalTree,'\',strcat('_',num2str(i),'_')));
    if i==1
        bioTree=bioTreeStack;
    else
        bioTree(end+1:end+size(bioTreeStack,2))=bioTreeStack;
    end
end
if iStack==5 || iStack==7 || iStack==9 || iStack==10
    load(strcat(dirFile,'\forbiddenImage.mat'))
    bioTree=bioTreeCorrection(bioTree,forbiddenImage);
else
    bioTree=bioTreeCorrection(bioTree);
end
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
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
function statisticLongTimeBacteriaWithFastProtocal(bioTree,dirFile)
[bioTree,allDataOri]=runaLL(bioTree);
save(strcat(dirFile,'\allDataOri'),'allDataOri');
correctionList=autoGetCorrectionList(bioTree,allDataOri);
dataCorrection(allDataOri,correctionList,[],[],dirFile,[]);
getFeatureParameter([],dirFile);
end
function longTimeStatistic(bioTree,dirFile)
longTimeStatisticForTree(bioTree,dirFile);
createLongTimeDivisionGraph(bioTree,dirFile);
close all
end
function clusterStatistic(bioTree,dirFile,step,dataOut)
% dirSave=strcat(dirFile,'\clusterAnalyse');
% mkdir(dirSave)
% cd(dirSave)
dirSave1=strcat(dirFile,'\longTimeResult');
mkdir(dirSave1)
cd(dirSave1)
bacteriaFrameInfo=getEachBacteriaInFrame(bioTree);
save('frameInfo','bacteriaFrameInfo');
bioTreeBasicInfo(bioTree,bacteriaFrameInfo,dataOut.timeInterval,dataOut.scaleInfo,dirFile);

allResult=dataAnalysisForLongTimeBacteria(bioTree);
saveas(gcf,'figure1.fig');
saveas(gcf,'figure1.tif');
close all

resultMap=traceStatisticMap(bioTree,bacteriaFrameInfo,20);
imshow(resultMap)
colormap(jet);
saveas(gcf,'traceInfo.tif');
close all

[velocityResult,lengthTime]=lengthAndVelocityVsTime(bioTree,bacteriaFrameInfo);
allResult.velocityResult=velocityResult;
allResult.lengthTime=lengthTime;
lengthAndVelocity(allResult.velocityResult,allResult.lengthTime);
save('allResult','allResult');
saveas(gcf,'length&velocity.fig')
saveas(gcf,'lengthAndvelocity.tif')
close all
% newClusterStatistic(bioTree,dirFile,bacteriaFrameInfo);
% 
% clusterTree=accumulateNetwork(bacteriaFrameInfo,bioTree,step);
% clusterTree=testClusterTree(clusterTree,dirSave,step,bioTree,bacteriaFrameInfo);
% save('frameInfo&clusterTree','bacteriaFrameInfo','clusterTree');
% [coverageTime,slopeTime]=degreeDistributionVsBacteriaCoverage(clusterTree,200,dirFile,bacteriaFrameInfo,'uw');
% save('degree&coverage','coverageTime','slopeTime');
% [coverageTime,slopeTime]=degreeDistributionVsBacteriaCoverage(clusterTree,200,dirFile,bacteriaFrameInfo,'w');
% save('weighted coverage','coverageTime','slopeTime');
% load(strcat(dirFile,'\frameInfo&clusterTree.mat'))
% dirDegreeGyrationFile=strcat(dirFile,'\logMapdegreeGyration');
% mkdir(dirDegreeGyrationFile)
% cd(dirDegreeGyrationFile)
% statisticDegreeGyration(bioTree,clusterTree,step,dirDegreeGyrationFile);
end
function lengthAndVelocity(velocityResult,lengthTime)
load('\\192.168.1.163\d\慢扫数据库\2014-02-11 F1 longTime(oilwater) Jzy\allResult.mat');
figure
hold on
subplot(1,2,2),plot(velocityResult(:,1),velocityResult(:,2),'r');
hold on;plot(allResult.velocityResult(:,1),allResult.velocityResult(:,2),'b');
xlabel('Time(min)');ylabel('velocity(um/s)');title('ave velocity vs Time');
subplot(1,2,1),plot(lengthTime(:,1),lengthTime(:,2),'r');
hold on;plot(allResult.lengthTime(:,1),allResult.lengthTime(:,2),'b');
xlabel('Time(min)');ylabel('length(um)');title('ave length vs Time');
end
function newClusterStatistic(bioTree,dirFile,bacteriaFrameInfo)
dirSave=strcat(dirFile,'\newCluster');
mkdir(dirSave)
cd(dirSave)
if nargin==2
    load('bacteriaFrameInfo');
else
    save('bacteriaFrameInfo','bacteriaFrameInfo');
end
getFigureMap(bioTree,bacteriaFrameInfo);
saveas(gcf,'twoWayProfile.fig');
[h1,h2,result]=directLinkMatrix(bacteriaFrameInfo,bioTree);
save('directResult','result')
saveas(h1,'mixBranchMap.fig')
saveas(h2,'RgMap.fig');
end
function treeNum=allFinalTreeSave(bioTree,saveFile,name)
saveFile=strcat(saveFile,'\',name);
bytesTemp=whos('bioTree');
treeNum=fix(bytesTemp.bytes/1900000000)+1;
savePoint=1;
savePoint=findFinalSavePoint(bioTree,treeNum,savePoint);
startFrame=1;
for iTree=1:treeNum
    if treeNum==1
        bioTreeStack=bioTree;
        bytesTmep=whos('bioTreeStack');
        if bytesTmep.bytes<2000000000
            save(saveFile,'bioTreeStack');
        else
            save(saveFile,'bioTreeStack','-v7.3');
        end
        return;
    end
    bioTreeStack=bioTree(startFrame:savePoint(iTree+1));
    saveFile1=strcat(saveFile,num2str(iTree));
    bytesTmep=whos('bioTreeStack');
    if bytesTmep.bytes<2000000000
        save(saveFile1,'bioTreeStack');
    else
        save(saveFile1,'bioTreeStack','-v7.3');
    end
    startFrame=savePoint(iTree+1)+1;
    bioTreeStack=[];
end
end
function savePoint=findFinalSavePoint(bioTree,treeNum,savePoint)
startFrame=1;
for iframe=100:100:size(bioTree,2)
    bioTreeTmep=bioTree(startFrame:iframe);
    bytesTemp=whos('bioTreeTmep');
    bytes=bytesTemp.bytes;
    pointNum=fix(bytes/1900000000);
    if pointNum==1
        savePoint=[savePoint,iframe-100];
        startFrame=iframe-99;
        if size(savePoint,2)==treeNum
            savePoint=[savePoint,size(bioTree,2)];
            return;
        end
    end
end
end
function degreeRadAll=statisticDegreeGyration(bioTree,clusterTree,step,dirFile)
beginOne=step+1-mod(200,step)+200;
timeInterval=[beginOne:step:size(clusterTree,2),size(clusterTree,2)];
n=0;
dirFile1=strcat(dirFile,'\tifFile');
dirFile2=strcat(dirFile,'\figFile');
mkdir(dirFile1)
mkdir(dirFile2)
for i=timeInterval
    n=n+1;
    degreeRad=logDegreeVsRadiusMap(bioTree,clusterTree,i);
    degreeRadAll.data{n}=degreeRad;
    degreeRadAll.num(n)=i;
    h=gca;
    saveas(h,strcat(dirFile1,'\',num2str(i),'.tif'));
    saveas(h,strcat(dirFile2,'\',num2str(i),'.fig'));
end
save(strcat(dirFile,'\degreeRadData'),'degreeRadAll');
end
function dataOut=scaleInforead(filein)
% input filein is the URL of the information.txt
fidin=fopen(filein,'r');
nline=0;
while ~feof(fidin) % 判断是否为文件末尾
    tline=fgetl(fidin); % 从文件读行
    nline=nline+1;
    if nline==5
        % find the line *** x : 1721 * 0.061728 : um ***
        info=textscan(tline,'%s%s%f%s%f');
        scaleInfo=info{5};
        xSize=info{3};
    end
    if nline==6
        info=textscan(tline,'%s%s%f%s%f');
        ySize=info{3};
    end
    %find the line *** Repeat T - 20000 times (5 sec) ***
    isRepeatLine=strfind(tline,'Repeat');
    if ~isempty(isRepeatLine)
        info=textscan(tline,'%s %s %s %f %s %s %s');
        if strcmp(info{7}{1}(1:end-1),'sec')
            timeInterval=str2num(info{6}{1}(2:end));
        end
        if strcmp(info{7}{1}(1:end-1),'ms')
           timeInterval=str2num(info{6}{1}(2:end))/1000;
        end
    end
end
fclose(fidin);
dataOut.rcSize=[ySize,xSize];
dataOut.scaleInfo=scaleInfo;
dataOut.timeInterval=timeInterval;
end