function bioTree=stepByStepTrackingCode(trackingCase)
dirFile=uigetdir();
dirOriginalImage=[dirFile,'\Tracking'];
dirTree=[dirFile,'\bioTreeResult'];
mkdir(dirTree)
% 针对路径下的original image，首先读取单张的大小，估测并等待得到smallTree，而后一起进行分析
% 首先读取一张照片，如果没有，则一直等待。
disp('waiting')
imageAmount=190;  % 2G以下系统可以运算的最多张数

% 让系统一直等待，除非相邻两张照片之间的间隔大于半个小时，则停止计算
waitingTime=0;
imageIndex=0;
seriousIndex=1;
while waitingTime<=30*60
    tic;nameList=dir(dirOriginalImage);
    newImage=numel(nameList)-3-imageIndex;
    wait1=toc;
    if newImage==0;
        pause(60)
        waitingTime=wait1+60+waitingTime;    % 如果没有新的Image,等待60s，并进行下一次循环
        clc
        disp(['waitingTime=',num2str(waitingTime)]);
    else
        waitingTime=0;
        if newImage>=imageAmount-imageIndex      % 若新的照片数目大于最大规定数目，则强制转化为等于最大规定数目，剩下的下次再计算
            newImage=imageAmount-imageIndex;
        end
        for iImage=1:newImage
            imageIndex=imageIndex+1;    % 如果有新的图像，则将其添加到imageStack
            maskImages(:,:,imageIndex)=load([dirOriginalImage,'\',nameList(imageIndex+3).name]);
        end
        if imageIndex==imageAmount   % 当imageIndex到达固定的张数时，计算maskImage,bioTree，并存储
            gainOneSmallTree(maskImages,seriousIndex,dirTree,imageAmount);
            imageIndex=0;
            seriousIndex=seriousIndex+1;
        end
    end
end
% 等待结束后，将剩余的imageStack进行计算。这样就得到了一系列小Tree
gainOneSmallTree(maskImages,seriousIndex,dirTree,imageAmount);

% 等所有实验结束后，combineSmallTreetoBigTree，并进行分割处理。
disp('combineSmallTreetoBigTree')
bioTree=combineSmallTreetoBigTree(seriousIndex,dirTree);

% 接下来进行图像的处理
% 如果是慢扫，明场追踪，则需要用topological Tracking
% 如果是荧光，则用新的图像分析方法
if strcmp(trackingCase,'flowCell')
    % flow cell tracking 程序
    bioTree=autoTidyBioTree(bioTree);
    nodeSortPre=[];
    [nodeSort,~]=findNode(bioTree);
    while ~isequal(nodeSortPre,nodeSort)
        nodeSortPre=nodeSort;
        bioTree=type1NodeReduction(bioTree,3,40);
        bioTree=type2NodeReduction(bioTree,2);
        bioTree=type3NodeReduction(bioTree,30);
        bioTree=type4NodeReduction(bioTree);
        [nodeSort,~]=findNode(bioTree);
    end
end
if strcmp(tracekingCase,'agarose')
    bioTree=autoTidyBioTree(bioTree);
    bioTree=agaroseBottomAutoSeg(bioTree);
    bioTree=type2NodeReduction(bioTree,2);
    bioTree=type3NodeReduction(bioTree,30);
    bioTree=type4NodeReduction(bioTree);
end

% bioTreeMeasure
bioTree=bioTreeDigHoleAndMeasure(bioTree,1,imageSize(1),imageSize(2));

% detect cutNum for each branch
bioTree=detectCutNumForBranch(bioTree);

testResult=[dirFile,'\testResult'];
mkdir(testResult)
branchList=bioTree{1}.branchList;
for iCoreBranch=1:size(branchList,1)
    if branchList(iCoreBranch,4)>8
        [linkMatrix,centroidInfo,~,leafNum,~,~]=generateOneBranchPhytree(bioTree,iCoreBranch,branchList(iCoreBranch,:),'normal');
        plotLinkTree(linkMatrix,centroidInfo,1:leafNum,[1,0,0],'normal',branchList(iCoreBranch,5));
        saveas(gcf,[testResult,num2str(iCoreBranch,'%03.0f'),'.tif'])
    end
end
end

%%  gainOneSmallTree
function gainOneSmallTree(maskImages,i,dirResultSave,imageAmount)
jobNumber=strcat('start Job', num2str(i));
folderName=strcat('t',num2str(i));
mkdir(dirResultSave,folderName);
disp(jobNumber);
if iJob~=1
    addpath(dirResultSave);
    preLastImage=load('lastImage');
    imageStack=cat(3,preLastImage.lastImage,maskImages);
end
lastImage=maskImages(:,:,end);
saveFile=strcat(dirResultSave,'\lastImage');
save(saveFile,'lastImage');
frameShift=preStackSize(i, imageAmount);
clear tempImage;
imageProcessingInfo.cropInfo=true(512,512);
maskImages=imageStack;
if strcmp(trackingCase,'')
[~,bioTree]=treeTrackingJob(imageStack,maskImages,1,size(maskImages,3),frameShift);
if i==1
    [~,backGroundPara]=backGroundCorrection([],[],'16bit');
    bioTree{1}.backGroundPara=backGroundPara;
    bioTree{1}.imageSize=[size(maskImages,1),size(maskImages,2)];
    bioTree{1}.imageProcessingInfo=imageProcessingInfo;
end
saveFile1=strcat(dirResultSave,'\',folderName,'\bioTree',num2str(i,'%02.0f'));
save(saveFile1,'bioTree')
end
end
function frameShift=preStackSize(iJob, stackSize)
if iJob==1
    frameShift=0;
else
    frameShift=(iJob-1)*stackSize-1;
end
end
%% myBioGraph_new2无root信息
function [bioTree,branchList,unCertainList,emptyNodeList]=myBiograph_new2(bioTree)
[bioTree,branchList]=stratRootSearching(bioTree);
bioTree=initClustering(bioTree,branchList);
CertainList=branchList;
while true
    [bioTree,unCertainList,floodNum]=nodeFlood(bioTree,CertainList,branchList);
    if isempty(unCertainList)
        break;
    end
    [bioTree,CertainList]=branchFinder(bioTree,unCertainList);
end
[bioTree,branchList]=getAllinList(bioTree,branchList);
if ~isempty(branchList)
    [emptyNodeList,isError]=emptyFinder(bioTree,branchList);
else
    emptyNodeList=[];
end
end
function [bioTree,branchList]=stratRootSearching(bioTree)
branchList=[];
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            trueIn=0;
            for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
                if bioTree{iframe}.node{iNode}.In{iIn}.isNode==true
                    trueIn=trueIn+1;
                end
            end
            if trueIn==0
                branchList=[branchList;[iframe,iNode,1]];
            end
        end
    end
end
end
function bioTree=initClustering(bioTree,branchList)
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            bioTree{iframe}.root{iroot}.xPos=[];
            bioTree{iframe}.root{iroot}.treePos=[];
            if bioTree{iframe}.root{iroot}.is2Node==true
                bioTree{iframe}.root{iroot}.branchIndex=[];
                bioTree{iframe}.root{iroot}.isBranch=false;
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for iLeaf=1:size(bioTree{iframe}.leavies,2)
            bioTree{iframe}.leavies{iLeaf}.xPos=[];
            bioTree{iframe}.leavies{iLeaf}.treePos=[];
            %             disp(iframe)
            %             disp(iLeaf)
            if bioTree{iframe}.leavies{iLeaf}.is2Node==true
                bioTree{iframe}.leavies{iLeaf}.branchIndex=[];
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            bioTree{iframe}.node{iNode}.xPos=[];
            bioTree{iframe}.node{iNode}.treePos=[];
            bioTree{iframe}.node{iNode}.branchIndex=[];
            bioTree{iframe}.node{iNode}.isHyperNode=false;
            for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
                bioTree{iframe}.node{iNode}.In{iIn}.branchIndex=[];
            end
        end
    end
end
for iList=1:size(branchList,1)
    branchInfo=branchList(iList,:);
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot=[];
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allLeaf=[];
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode=[];
    bioTree{branchInfo(1)}.node{branchInfo(2)}.branchIndex=iList;
    for iIn=1:size(bioTree{branchInfo(1)}.node{branchInfo(2)}.In,2)
        bioTree{branchInfo(1)}.node{branchInfo(2)}.In{iIn}.branchIndex=iList;
    end
end
end
function [bioTree,unCertainList,floodNum]=nodeFlood(bioTree,nodeList,branchList)
unCertainList=[];
floodNum=0;
for iList=1:size(nodeList,1)
    isStart=true;
    allRoot=[];
    allLeaf=[];
    nodeInfo=nodeList(iList,:);
    branchIndex=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex;
    if branchIndex==258
        p=1;
    end
    nodeListS=nodeInfo;
    searching=true;
    while searching
        [bioTree,subNodeList,unCertlist,rootList,leafList,searching]=subNodeSearhing(bioTree,nodeListS,branchIndex,branchList,isStart);
        nodeListS=subNodeList;
        unCertainList=[unCertainList;unCertlist];
        allRoot=[allRoot;rootList];
        allLeaf=[allLeaf;leafList];
        isStart=false;
    end
    branchInfo=branchList(branchIndex,:);
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot=[bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot;(allRoot)];
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allLeaf=[bioTree{branchInfo(1)}.node{branchInfo(2)}.allLeaf;(allLeaf)];
    allNode=bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode;
    floodNum=floodNum+size(allNode,1)+size(allRoot,1)+size(allLeaf,1);
end
end
function [bioTree,subNodeList,unCertlist,rootList,leafList,searching]=subNodeSearhing(bioTree,nodeList,branchIndex,branchList,isStart)
subNodeList=[];
unCertlist=[];
rootList=[];
leafList=[];
for iList=1:size(nodeList,1)
    nodeInfo=nodeList(iList,:);
    [isCert,unCertNode]=isCertain(bioTree,nodeInfo,isStart);
    unCertlist=[unCertlist;unCertNode];
    if isCert==true || isStart==true
        branchInfo=branchList(branchIndex,:);
        bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode=[bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode;nodeInfo];
        if isempty( bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex)
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex=branchIndex;
        end
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==true
                subNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                bioTree{subNodeInfo(1)}.node{subNodeInfo(2)}.In{subNodeInfo(3)}.branchIndex=branchIndex;
                if isempty(bioTree{subNodeInfo(1)}.node{subNodeInfo(2)}.branchIndex)
                    subNodeList=[subNodeList;subNodeInfo];
                end
            end
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==false
                subLeafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                bioTree{subLeafInfo(1)}.leavies{subLeafInfo(2)}.branchIndex=branchIndex;
                leafList=[leafList;subLeafInfo];
            end
        end
        for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
                rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
                rootList=[rootList;rootInfo];
                bioTree{rootInfo(1)}.root{rootInfo(2)}.branchIndex=branchIndex;
            end
        end
    end
end
subNodeList=listRuduction(subNodeList);
if ~isempty(subNodeList)
    searching=true;
else
    searching=false;
end
end
function [isCert,unCertNode]=isCertain(bioTree,nodeInfo,isStart)
indexIn=[];
unCertNode=[];
if isStart==true;
    isCert=true;
    return;
end
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1
    isCert=true;
    return;
end
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if  bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        if isempty(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.branchIndex)
            isCert=false;
            return;
        end
    end
    indexIn=[indexIn;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.branchIndex];
end
indexIn=unique(indexIn(:,1)');
if  size(indexIn,2)==1
    isCert=true;
    return;
else
    isCert=false;
    unCertNode=nodeInfo;
    return;
end
end
function [bioTree,CertainList]=branchFinder(bioTree,unCertainList)
for iList=1:size(unCertainList,1)
    rootList=[];
    leafList=[];
    finder=[];
    bList=[];
    nodeInfo=unCertainList(iList,:);
    for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
            parentNodeInfo2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
            bIndex=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.branchIndex;
            %             EdgeTemp=bioTree{parentNodeInfo2(1)}.node{parentNodeInfo2(2)}.Out{parentNodeInfo2(3)}.traceInfo.measurment{end};
            pixelIdxList=bioTree{parentNodeInfo2(1)}.node{parentNodeInfo2(2)}.Out{parentNodeInfo2(3)}.traceInfo.pixelIdxList{end};
            imageSize=bioTree{1}.imageSize;
            [~,BWImage]=idx2Xy(pixelIdxList,imageSize);
            BWImage=imfill(BWImage,'holes');
            edgeSize=numel(BWImage(BWImage==1));
            %edgeSize=getEdgeSize(EdgeTemp);
            finder=[finder;[bIndex,edgeSize]];
        end
    end
    branchIndex=whoisBigger(finder);
    for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
            rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
            rootList=[rootList;rootInfo];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.branchIndex=branchIndex;
        end
    end
    for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==false
            subLeafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
            bioTree{subLeafInfo(1)}.leavies{subLeafInfo(2)}.branchIndex=branchIndex;
            leafList=[leafList;subLeafInfo];
        end
    end
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex=branchIndex;
    bList=[bList;branchIndex];
end
CertainList=unCertainList;
end
function [bioTree,branchList]=getAllinList(bioTree,branchList)
sizeList=[];
for iList=1:size(branchList,1)
    branchInfo=branchList(iList,:);
    allRoot=listRuduction(bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot);
    allLeaf=listRuduction(bioTree{branchInfo(1)}.node{branchInfo(2)}.allLeaf);
    allNode=listRuduction(bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode);
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot=allRoot;
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allLeaf=allLeaf;
    bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode=allNode;
    sizeBracnch=size(allRoot,1)+size(allLeaf,1)+size(allNode,1);
    sizeList=[sizeList;sizeBracnch];
end
branchList=[branchList,sizeList];
end
function [emptyNodeList,isError]=emptyFinder(bioTree,branchList)
nCount=0;
emptyNodeList=[];
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                if isempty(bioTree{iframe}.root{iroot}.branchIndex)
                    emptyNodeList=[emptyNodeList;[iframe,iroot]];
                end
                nCount=nCount+1;
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for iLeaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{iLeaf}.is2Node==true
                if isempty(bioTree{iframe}.leavies{iLeaf}.branchIndex)
                    emptyNodeList=[emptyNodeList;[iframe,iLeaf]];
                end
                nCount=nCount+1;
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if isempty(bioTree{iframe}.node{iNode}.branchIndex)
                emptyNodeList=[emptyNodeList;[iframe,iNode]];
            end
            nCount=nCount+1;
        end
    end
end
if isempty(emptyNodeList) && sum(branchList(:,4))==nCount
    isError=false;
    disp('Lucky! There is No Error');
else
    isError=true;
    disp('No !! Clustering Error');
end
end
function nodeListafter=listRuduction(nodeListbefore)
nodeListafter=[];
while ~isempty(nodeListbefore)
    nodeListafter=[nodeListafter;nodeListbefore(1,:)];
    tempListFrame=nodeListbefore(:,1)-nodeListbefore(1,1);
    tempListNode=nodeListbefore(:,2)-nodeListbefore(1,2);
    nodeListbefore=nodeListbefore(tempListFrame~=0|tempListNode~=0,:);
end
end
function branchIndex=whoisBigger(finder)
uIndex=unique(finder(:,1)');
cEdgeSize=[];
for i=1:size(uIndex,2)
    cEdgeSize=[cEdgeSize;[uIndex(i),sum(finder(finder(:,1)==uIndex(i),2))]];
end
[~,mIndex]=max(cEdgeSize(:,2));
branchIndex=cEdgeSize(mIndex,1);
end
function egdeSize=getEdgeSize(EgdeTemp)
fillArea=[];
for i=1:size(EgdeTemp,1)
    fillArea=[fillArea,EgdeTemp(i).FilledArea];
end
egdeSize=sum(fillArea);
end
%% flowcell tracking auto segementation

%%% type1 node reduction
function [bioTree,usefulNode]=type1NodeReduction(bioTree,strength,finalOneToOneMatchThreShold)
% ##defalt finalOneToOneMatchThreShold=25 could be reach 40(in ASolution l120)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type1NodeinFrame(bioTree{iframe},iframe,bioTree);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,usefulNode]=fullNodeType1Tracking(bioTree,nodeList(iList,:),strength,finalOneToOneMatchThreShold);
                nodeList(iList,4)=usefulNode;
                if usefulNode==0;
                    bioTree=type1NodeLinker(bioTree,nodeInfoIn,traceInfo);
                end
                nodeCount=nodeCount+1;
            end
            cantDivideNode=nodeList(:,4);
            nodeList(cantDivideNode==1,:)=[];
            if ~isempty(nodeList)
                bioTree=removeType1Node(bioTree,nodeList);
            end
        end
    end
end
end
function nodeList=type1NodeinFrame(bioTreeFrame,iframe,bioTree)
nodeList=[];
for iNode=1:size(bioTreeFrame.node,2)
    if size(bioTreeFrame.node{iNode}.In,2)>=1 && size(bioTreeFrame.node{iNode}.Out,2)==1
        %         if bioTreeFrame.node{iNode}.branchIndex~=12
        %             break
        %         end
        nodeInfo=[iframe,iNode];
        usefulNode=0;
        pixelIdxListIn=getInputMask(bioTree,nodeInfo);
        for iIn=1:size(pixelIdxListIn,2)
            [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
            CC=bwconncomp(BWImage);
            regionNum=CC.NumObjects;
            if regionNum>=2
                usefulNode=1;
                break
            end
        end
        if usefulNode==0
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function bioTree=type1NodeLinker(bioTree,nodeInfoIn,traceInfo)
newNodeIn=0;
for iIn=1:size(nodeInfoIn,2)
    if isempty(traceInfo{iIn}.isBreakNode)
        if nodeInfoIn{iIn}.isNode==false
            rootInfo=nodeInfoIn{iIn}.rootInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 leafIndex=size(bioTree{leafInfo(1)}.leavies,2);
                %                 leafInfo(2)=leafIndex+1;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
                if isempty(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList)
                    bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{1}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
                end
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
            end
        end
        if nodeInfoIn{iIn}.isNode==true
            nodeInfo=nodeInfoIn{iIn}.nodeInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 leafIndex=size(bioTree{leafInfo(1)}.leavies,2);
                %                 leafInfo(2)=leafIndex+1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{end};
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo= leafInfo;
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
            end
        end
    end
    if traceInfo{iIn}.isBreakNode==1
        newNodeIn=newNodeIn+1;
        if newNodeIn==1
            breakNodeInfo=[traceInfo{iIn}.breakNodeInfo(1),size(bioTree{traceInfo{iIn}.breakNodeInfo(1)}.node,2)+1];
        end
        if nodeInfoIn{iIn}.isNode==false
            rootInfo=nodeInfoIn{iIn}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=true;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[breakNodeInfo,newNodeIn];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=[];
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.isNode=false;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.rootInfo=rootInfo;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
        else
            if nodeInfoIn{iIn}.isNode==true
                preNodeInfo=nodeInfoIn{iIn}.nodeInfo;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.is2Node=true;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.nodeInfo=[breakNodeInfo,newNodeIn];
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.leafInfo=[];
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.isNode=1;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.In{newNodeIn}.nodeInfo=preNodeInfo;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
            end
        end
        if  traceInfo{iIn}.is2Node==false
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.is2Node=false;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.leafInfo=leafInfo;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.traceInfo.pixelIdxList=traceInfo{iIn}.afterBreakTrace.pixelIdxList;
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.nodeInfo=[];
            bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[breakNodeInfo,1];
        else
            if traceInfo{iIn}.is2Node==true
                nodeInfo=traceInfo{iIn}.nodeInfo;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.is2Node=true;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.nodeInfo=nodeInfo;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.traceInfo.pixelIdxList=traceInfo{iIn}.afterBreakTrace.pixelIdxList;
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.Out{1}.leafInfo=[];
                bioTree{breakNodeInfo(1)}.node{breakNodeInfo(2)}.state=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=true;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[breakNodeInfo,1];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[];
            end
        end
    end
end
end
function nodeInfoIn=getInputInfo(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        nodeInfoIn{iIn}.isNode=true;
        nodeInfoIn{iIn}.nodeInfo= nodeInfo_pre;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        nodeInfoIn{iIn}.isNode=false;
        nodeInfoIn{iIn}.rootInfo= rootInfo;
    end
end
end
function bioTree=removeType1Node(bioTree,nodeList)
newList=[];
countNode=0;
iframe=nodeList(1,1);
for iList=1:size(nodeList,1)
    nodeInfo=nodeList(iList,:);
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}=[];
end
for iNode=1:size(bioTree{iframe}.node,2)
    if ~isempty(bioTree{iframe}.node{iNode})
        countNode=countNode+1;
        newList=[newList;[iNode,countNode]];
    end
end

if ~isempty(newList)
    for iList=1:size(newList,1)
        newNode{newList(iList,2)}=bioTree{iframe}.node{newList(iList,1)};
    end
    bioTree{iframe}.node=newNode;
else
    bioTree{iframe}.node=[];
end
if ~isempty(bioTree{iframe}.node)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==true
                nodeInfopre=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                bioTree{nodeInfopre(1)}.node{nodeInfopre(2)}.Out{nodeInfopre(3)}.nodeInfo=[iframe,iNode,iIn];
            end
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                rootInfopre=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
                bioTree{rootInfopre(1)}.root{rootInfopre(2)}.nodeInfo=[iframe,iNode,iIn];
            end
        end
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                nodeInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                bioTree{nodeInfonext(1)}.node{nodeInfonext(2)}.In{nodeInfonext(3)}.nodeInfo=[iframe,iNode,iOut];
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                leafInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                bioTree{leafInfonext(1)}.leavies{leafInfonext(2)}.nodeInfo=[iframe,iNode,iOut];
            end
        end
    end
end
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
% BWImageGain=imfill(BWImageGain,'holes');
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
        end
    end
end
end

%%% type2 node reduction
function bioTree=type2NodeReduction(bioTree,strength)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type2NodeinFrame(bioTree{iframe},iframe,bioTree);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,canDivideorNot]=fullNodeType2Tracking(bioTree,nodeList(iList,:),strength);
                nodeList(iList,4)=canDivideorNot;
                if canDivideorNot==1
                    bioTree=type2NodeLinker(bioTree,nodeInfoIn,traceInfo);
                end
                nodeCount=nodeCount+1;
            end
            canDivideNode=nodeList(:,4);
            nodeList(canDivideNode==0,:)=[];
            if ~isempty(nodeList)
                bioTree=removeType2Node(bioTree,nodeList);
            end
        end
    end
end
end
function nodeList=type2NodeinFrame(bioTreeFrame,iframe,bioTree)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if size(bioTreeFrame.node{iNode}.In,2)==size(bioTreeFrame.node{iNode}.Out,2) && size(bioTreeFrame.node{iNode}.In,2)<=8
            nodeInfo=[iframe,iNode];
            usefulNode=0;
            pixelIdxListIn=getInputMask(bioTree,nodeInfo);
            for iIn=1:size(pixelIdxListIn,2)
                [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
                CC=bwconncomp(BWImage);
                regionNum=CC.NumObjects;
                if regionNum>=2
                    usefulNode=1;
                    break
                end
            end
            if usefulNode==0
                nodeList=[nodeList;[iframe,iNode]];
            end
        end
    end
end
end
function bioTree=type2NodeLinker(bioTree,nodeInfoIn,traceInfo)
for iIn=1:size(nodeInfoIn,2)
    if nodeInfoIn{iIn}.isNode==false
        rootInfo=nodeInfoIn{iIn}.rootInfo;
        if traceInfo{iIn}.is2Node==false;
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
            %             if leafInfo(1)~=rootInfo(1)
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            %             end
            bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
        end
        if traceInfo{iIn}.is2Node==true;
            nodeInfoT=traceInfo{iIn}.nodeInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
        end
    end
    if nodeInfoIn{iIn}.isNode==true
        nodeInfo=nodeInfoIn{iIn}.nodeInfo;
        if traceInfo{iIn}.is2Node==false;
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
            %             if leafInfo(1)~=nodeInfo(1)
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            %             end
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo= leafInfo;
        end
        if traceInfo{iIn}.is2Node==true;
            nodeInfoT=traceInfo{iIn}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
            bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
        end
    end
end
end
function bioTree=removeType2Node(bioTree,nodeList)
newList=[];
countNode=0;
iframe=nodeList(1,1);
for iList=1:size(nodeList,1)
    nodeInfo=nodeList(iList,:);
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}=[];
end
for iNode=1:size(bioTree{iframe}.node,2)
    if ~isempty(bioTree{iframe}.node{iNode})
        countNode=countNode+1;
        newList=[newList;[iNode,countNode]];
    end
end

if ~isempty(newList)
    for iList=1:size(newList,1)
        newNode{newList(iList,2)}=bioTree{iframe}.node{newList(iList,1)};
    end
    bioTree{iframe}.node=newNode;
else
    bioTree{iframe}.node=[];
end
if ~isempty(bioTree{iframe}.node)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==true
                nodeInfopre=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                bioTree{nodeInfopre(1)}.node{nodeInfopre(2)}.Out{nodeInfopre(3)}.nodeInfo=[iframe,iNode,iIn];
            end
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                rootInfopre=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
                bioTree{rootInfopre(1)}.root{rootInfopre(2)}.nodeInfo=[iframe,iNode,iIn];
            end
        end
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                nodeInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                bioTree{nodeInfonext(1)}.node{nodeInfonext(2)}.In{nodeInfonext(3)}.nodeInfo=[iframe,iNode,iOut];
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                leafInfonext=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                bioTree{leafInfonext(1)}.leavies{leafInfonext(2)}.nodeInfo=[iframe,iNode,iOut];
            end
        end
    end
end
end

%%% type3 node reduction
function bioTree=type3NodeReduction(bioTree,threShold)
nodeCount=0;
for iframe=1:size(bioTree,2)
    dispFrame(iframe)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type3NodeinFrame(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,inNum,outNum,canDivideorNot]=fullNodeType3Tracking(bioTree,nodeList(iList,:),threShold);
                if canDivideorNot==1
                    bioTree=type3NodeLinker(bioTree,nodeInfoIn,traceInfo,inNum,nodeList(iList,:));
                    bioTree=removeType3Node(bioTree,nodeList(iList,:),inNum,outNum);
                end
                nodeCount=nodeCount+1;
                %                 disp(nodeCount);
            end
        end
    end
end
end
function nodeList=type3NodeinFrame(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if (size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==2)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==2 && size(bioTreeFrame.node{iNode}.Out,2)>=3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==4)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==5)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==5 && size(bioTreeFrame.node{iNode}.Out,2)==4)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==2)
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function bioTree=type3NodeLinker(bioTree,nodeInfoIn,traceInfo,inNum,nodeInfo)
if ~(numel(inNum)==1 && inNum==0)
    for iIn=1:numel(inNum)
        if nodeInfoIn{inNum(iIn)}.isNode==false
            rootInfo=nodeInfoIn{inNum(iIn)}.rootInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
                %                 if leafInfo(1)~=rootInfo(1)
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 end
                bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
            end
        end
        if nodeInfoIn{inNum(iIn)}.isNode==true
            nodeInfo=nodeInfoIn{inNum(iIn)}.nodeInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
                %                 if leafInfo(1)~=nodeInfo(1)
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                %                 end
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
            end
        end
    end
end
if numel(inNum)==1 && inNum==0
    rootNum=size(bioTree{1,nodeInfo(1)}.root,2);
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node=traceInfo{1,1}.is2Node;
    if bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node==0
        leafInfo=traceInfo{1,1}.leafInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.leafInfo=leafInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.nodeInfo=[];
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.rootPixelDetail=traceInfo{1,1}.traceInfo.pixelIdxList{1};
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=traceInfo{1,1}.traceInfo.pixelIdxList;
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.rootInfo=[nodeInfo(1),rootNum+1];
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.nodeInfo=[];
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.is2Node=0;
    end
    if bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node==1
        newNodeInfo=traceInfo{1,1}.nodeInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.nodeInfo=newNodeInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.leafInfo=[];
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.isNode=0;
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.nodeInfo=[];
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.rootInfo=[nodeInfo(1),rootNum+1];
    end
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.rootPixelDetail=traceInfo{1,1}.traceInfo.pixelIdxList{1};
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=traceInfo{1,1}.traceInfo.pixelIdxList;
end
end
function bioTree=removeType3Node(bioTree,nodeInfo,inNum,outNum)
if ~(numel(inNum)==1 && inNum==0)
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In(inNum)=[];
    for i=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.isNode==0
            rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo(3)=i;
        end
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.isNode==1
            NewnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.nodeInfo;
            bioTree{NewnodeInfo(1)}.node{NewnodeInfo(2)}.Out{1,NewnodeInfo(3)}.nodeInfo(3)=i;
        end
    end
end
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out(outNum)=[];
for i=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.is2Node==0
        leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.leafInfo;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo(3)=i;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.is2Node==1
        NewnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.nodeInfo;
        bioTree{NewnodeInfo(1)}.node{NewnodeInfo(2)}.In{1,NewnodeInfo(3)}.nodeInfo(3)=i;
    end
end
end

%%% type4 node reduction
function bioTree=type4NodeReduction(bioTree)
nodeCount=0;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=type4NodeinFrame(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                %nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                nodeInfo=nodeList(iList,:);
                pixelIdxListIn=getInputMask(bioTree,nodeInfo);
                iInProps=0;
                for iIn=1:size(pixelIdxListIn,2)
                    [~,BWImage]=idx2Xy(pixelIdxListIn{iIn},bioTree{1}.imageSize);
                    CC=bwconncomp(BWImage);
                    regionNum=CC.NumObjects;
                    if regionNum>=2
                        iInProps=1;
                    end
                end
                if iInProps==1
                    break
                end
                imageSize=bioTree{1}.imageSize;
                traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.traceInfo.pixelIdxList;
                pictureOut=false(imageSize);
                pictureOut(traceInfo{1})=true;
                cc=bwconncomp(pictureOut);
                if cc.NumObjects==2
                    outPixelIdxList=cc.PixelIdxList;
                    if numel(traceInfo)>1
                        nodeNum=size(bioTree{nodeInfo(1)+1}.node,2);
                        bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1};
                        bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.traceInfo.pixelIdxList(1)=[];
                        if bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.is2Node==1
                            newNodeInfo=bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.nodeInfo;
                            bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{newNodeInfo(3)}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,1];
                        end
                        if bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.is2Node==0
                            newLeafInfo=bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.Out{1}.leafInfo;
                            bioTree{newLeafInfo(1)}.leavies{newLeafInfo(2)}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,1];
                        end
                        for iOut=1:2
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=[nodeInfo(1)+1,nodeNum+1,iOut];
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList=[];
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                        end
                        for iIn=1:2
                            bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.In{1,iIn}.isNode=1;
                            bioTree{nodeInfo(1)+1}.node{1,nodeNum+1}.In{1,iIn}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                        end
                    end
                    if numel(traceInfo)==1
                        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node==1
                            newNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
                            inNum=size(bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In,2);
                            for iIn=1:2
                                if iIn==2
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.isNode=1;
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.rootInfo=[];
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,iIn+inNum-1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                                else
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.isNode=1;
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.rootInfo=[];
                                    bioTree{newNodeInfo(1)}.node{newNodeInfo(2)}.In{1,newNodeInfo(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                                end
                            end
                            for iOut=1:2
                                if iOut==2
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=[newNodeInfo(1),newNodeInfo(2),iOut+inNum-1];
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                                else
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.is2Node=1;
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.leafInfo=[];
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.nodeInfo=newNodeInfo;
                                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,iOut}.traceInfo.pixelIdxList{1}=outPixelIdxList{iOut};
                                end
                            end
                        end
                    end
                end
                nodeCount=nodeCount+1;
                %             disp(nodeCount)
            end
        end
    end
end
fprintf('\n')
end
function nodeList=type4NodeinFrame(bioTreeFrame,iframe)
nodeList=[];
for iNode=1:size(bioTreeFrame.node,2)
    if size(bioTreeFrame.node{iNode}.In,2)==1 && size(bioTreeFrame.node{iNode}.Out,2)==1
        nodeList=[nodeList;[iframe,iNode]];
    end
end
end

%% detect cutNum for eachBranch
function bioTree=detectCutNumForBranch(bioTree)
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
for iCoreBranch=1:size(branchList,1)
    iBranchInfo=branchList(iCoreBranch,:);
    allNode=bioTree{iBranchInfo(1)}.node{iBranchInfo(2)}.allNode;
    [~,order]=sort(allNode(:,1));
    allNode=allNode(order,:);
    branchList(iCoreBranch,5)=numel(bioTree);
    for iNode=1:size(allNode,1)
        nodeInfo=bioTree{allNode(iNode,1)}.node{allNode(iNode,2)};
        isDivision=isDivsionNode(bioTree,nodeInfo);
        if isDivision==0;
            branchList(iCoreBranch,5)=nodeInfo(1)-1;
            break
        end
    end
end
bioTree{1}.branchList=branchList;
end
function isDivision=isDivsionNode(bioTree,nodeInfo)
isDivision=0;
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2
    frameOut=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out;
    bacteriaSize=[frameOut{1}.traceInfo.measurment{1}.FilledArea,frameOut{2}.traceInfo.measurment{1}.FilledArea];
    if ~(max(bacteriaSize)-min(bacteriaSize)>150 || max(bacteriaSize)>min(bacteriaSize)*2 )
        isDivision=1;
    end
end
end
function plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,type,cutNum)
% generally linkMatrix means the relationship for [N, L], N means all
% node, L means all leaf, centroidInfo here could be changed for the second
% column. The first column means the frame that the node begin\
% color=[1,0,1];
leafNum=numel(leafYCoo);
nodeNum=size(linkMatrix,1)-leafNum;
centroidInfo(nodeNum+1:end,2)=leafYCoo;
leafOrder=1:leafNum;
for i=1:nodeNum
    iLine=linkMatrix(nodeNum+1-i,nodeNum+1-i:end);
    centroid2=centroidInfo(:,2);
    centroid2=centroid2(nodeNum+1-i:end);
    leafPos=(centroid2(iLine~=0));
    centroidInfo(nodeNum+1-i,2)=mean(leafPos(~isnan(leafPos)));
end
% figure;
color=[1,0,0];
% new add map size
figure1=figure;
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.288229166666667 0.74854189336235],'FontSize',20);
% ylim(axes1,[0 118])
% xlim(axes1,[0 12.6])
view(axes1,[90 -90]);
box(axes1,'on');
hold all
% Create xlabel
xlabel('time(h)','VerticalAlignment','bottom','Rotation',90,...
    'HorizontalAlignment','center','FontSize',20);
scrsz=get(0,'ScreenSize');
set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
set(gcf, 'PaperPositionMode', 'auto');

for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1)
        if linkMatrix(i,j)~=0
            if nargin==5 && strcmp(type,'aging') 
                linkTwoPointsAging(centroidInfo(i,:),centroidInfo(j,:),color,nodeNum,i,j,linkMatrix(i,j));
            else
                linkTwoPoints(centroidInfo(i,:),centroidInfo(j,:),color,nodeNum,i,j,linkMatrix(i,j),cutNum);
            end
        end
    end
end
end
function linkTwoPoints(c1,c2,color,nodeNum,i,j,para,cutNum)
c1(1)=c1(1)/2;
c2(1)=c2(1)/2;
cutNum=cutNum/2;
% c1(1)=c1(1)*4;
% c2(1)=c2(1)*4;
% support c1(x)<c2(x),link c1-c3-c2

% debug plus para(linkMatrix(i,j))
% if para==-inf
%     color=[0,0,0];
% else
%     colorAll=colormap(jet(8));
%     color=colorAll(para,:);
% end
% if para~=0
%     para=ceil((para-1000)/1000*2);
%     color=color(para,:);
%     if para<=4
%         color=[0,0,1];
%     else
%         color=[1,0,0];
%     end
% end
if c1(1)>c2(1)
    c=c2;
    c2=c1;
    c1=c;
end
hold on;
if c1(1)>cutNum
    return
end
% link line
% c3=[c1(1),c2(2)];
% line([c1(1),c3(1),c2(1)],[c1(2),c3(2),c2(2)],'Color',color,'LineWidth',1);

% link with ellipse
c3=[c2(1),c1(2)];
axisA=abs(c2(1)-c1(1));
axisB=abs(c2(2)-c1(2));
if c2(1)>cutNum
    x=c1(1):0.01:cutNum;
else
    x=c1(1):0.01:c2(1);
end
if c2(2)<c1(2)
    y=c3(2)-axisB*(1-((x-c3(1))/axisA).^2).^(1/2);
else
    y=c3(2)+axisB*(1-((x-c3(1))/axisA).^2).^(1/2);
end
line(x,y,'Color',color,'LineWidth',1)


maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker','^','MarkerSize',12,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<=cutNum
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        end
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<=cutNum
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end
function linkTwoPointsAging(c1,c2,color,nodeNum,i,j,para)
% c1(1)=c1(1)/1200;
% c2(1)=c2(1)/1200;
% support c1(x)<c2(x),link c1-c3-c2

% debug plus para(linkMatrix(i,j))
if para==-inf
    color=[0,0,0];
else
    colorAll=colormap(jet(10));
    color=colorAll(para,:);
end
    
if c1(1)>c2(1)
    c=c2;
    c2=c1;
    c1=c;
end
c3=[c1(1),c2(2)];
hold on;
line([c1(1),c3(1),c2(1)],[c1(2),c3(2),c2(2)],'Color',color,'LineWidth',1);
maker='.';
size=16;
if i==1
    plot(c1(1),c1(2),'Marker','^','MarkerSize',12,'Color',[1,0,0]);
else
    if j<=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
    end
    if j>=nodeNum
        plot(c1(1),c1(2),'Marker',maker,'MarkerSize',size,'Color',[205,155,29]/255);
        if c2(1)<360
            plot(c2(1),c2(2),'Marker',maker,'MarkerSize',size,'Color',[85,107,47]/255);
        end
    end
end
end
