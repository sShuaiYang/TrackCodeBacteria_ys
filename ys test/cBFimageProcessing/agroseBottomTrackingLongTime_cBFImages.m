function  agroseBottomTrackingLongTime_cBFImages()
%移植longtimeProtocal的code进行处理
dirFile='\\192.168.1.14\e\2020-04-23 PAO1(B226)_IP31_100x 08agar Temp30 ys';
fieldList=dir(dirFile);
for iField = 17:36%numel(fieldList)-2
    t0 = clock;
    if ~(strcmp(fieldList(iField+2).name(1:5),'field') ) && strcmp(fieldList(iField+2).name(end-2:end),'mat')
        continue
    end

    fieldIdx=str2double(fieldList(iField+2).name(end-3:end)); 
    disp(['loading',32,fieldList(iField+2).name]);
    disp(compose("Start time %d-%d-%d %.2d:%.2d:%.0f.",t0));
    
    dirField=[dirFile,'\',fieldList(iField+2).name];
    % flow cell tracking
    disp('treeProcessing')
    t1 = clock;
    bioTree=batchTreeTrackingJob(dirField);
    disp(['bioTree generation 耗时:',num2str(etime(clock,t1))]);

    disp(' longtime bioTree Analysis')
    tic;
    bioTree=longTimeBioTreeAnalysis(bioTree);
    toc;
    
    disp('bioTree Tidy and Reduction')
    tic;
    bioTree=autoTidyBioTree(bioTree);
    toc;
    %%  bioTreeMeasure
    tic;
    bioTree{1}.fieldNum=fieldIdx;
    disp('Now bioTree mearsure');
    bioTree=bioTreeDigHoleAndMeasure(bioTree,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
    toc;
    
    disp('Try detect cutNum for each branch and Two point aging')
    tic;
    try
    bioTree=detectCutNumForBranch(bioTree);
    catch
        warning('Problem using the function detectCutNumForBranch. Skipped.');    
    end
    try
        bioTree=bioTreeTwoPointAging(bioTree);
    catch
        warning('Problem using the function bioTreeTwoPointAging. Skipped.');
    end
    toc;
    
    disp ('genealogicalTreeGeneration')
    pedigreeSave=[dirField,'\pedigree'];
    mkdir(pedigreeSave)
    
    load([dirField,'\Tracking\frameInfo.mat']);
    bioTreeInitialTime = frameInfo(1,1:6);
    frameInfoNew(1:2:numel(bioTree)-1,:) = frameInfo(1:numel(bioTree)/2,:);
    frameInfoNew(2:2:numel(bioTree),:) = frameInfo(1:numel(bioTree)/2,:);
    bioTreeTimer = zeros(1,numel(bioTree));
    for i = 1:numel(bioTree)
        bioTreeTimer(i) = etime(frameInfoNew(i,1:6),bioTreeInitialTime)/60;   %min
    end
    bioTree{1}.bioTreeInitialTime = bioTreeInitialTime;
    bioTree{1}.bioTreeTimer=bioTreeTimer;
    save([dirField,'\','bioTree'],'bioTree','-v7.3');
    
    [~] = genealogicalTreeGeneration(bioTree,pedigreeSave);
    [~,~] = generationGrowthRateCal_agrose_new(bioTree,dirField);
    %     growthDataAll{iField}=growthData;
    clear bioTree 
    t1 = clock;
    disp(compose("End time %d-%d-%d %.2d:%.2d:%.0f.",t1))
    disp(['共耗时:',num2str(etime(t1,t0)),'秒']);
end
% save([dirAllFile,'\growthDataAll.mat'],'growthDataAll');
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
    disp('Lucky! There is No Error');
end
[bioTree,branchList]=addIrrelevantBacteria(bioTree,branchList);
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
function   [bioTree,subNodeList,unCertlist,rootList,leafList,searching]=subNodeSearhing(bioTree,nodeList,branchIndex,branchList,isStart)
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
function [bioTree,branchList]=addIrrelevantBacteria(bioTree,branchList)
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        if bioTree{iframe}.root{iRoot}.is2Node==0
            branchNum=size(branchList,1)+1;
            branchList=[branchList;iframe,iRoot,0,2];
            bioTree{iframe}.root{iRoot}.xPos=[];
            bioTree{iframe}.root{iRoot}.treePos=[];
            bioTree{iframe}.root{iRoot}.branchIndex=branchNum;
            bioTree{iframe}.root{iRoot}.isHyperNode=0;
            bioTree{iframe}.root{iRoot}.allRoot=[iframe,iRoot];
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            bioTree{iframe}.root{iRoot}.allNode=[];
            bioTree{iframe}.root{iRoot}.allLeaf=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.xPos=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.treePos=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.branchIndex=branchNum;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.xPos.isHyperNode=0;
        end
    end
end
end
%%
function [ bioTree ] = bioTreeAllInfoGet_cBF( bioTree,dirFile)
%BIOTREEALLINFOGET Summary of this function goes here
load([dirFile,'\Tracking\frameInfo.mat']);
bioTreeInitialTime=frameInfo(1,1:6);
frameInfoNew(1:2:numel(bioTree)-1,:)=frameInfo(1:numel(bioTree)/2,:);
frameInfoNew(2:2:numel(bioTree),:)=frameInfo(1:numel(bioTree)/2,:);
% frameInfoNew=frameInfo;
% bestPositionAccumulation=bioTree{1}.bestPositionAccumulation;
for i=1:numel(bioTree)
    bioTreeTimer(i)=etime(frameInfoNew(i,1:6),bioTreeInitialTime)/60;   %min
end
bioTree{1}.bioTreeInitialTime=bioTreeInitialTime;
bioTree{1}.bioTreeTimer=bioTreeTimer;

% load([dirFile,'\Tracking\frameInfo.mat'])

bioTreeFrame=[];
nameList=dir([dirFile,'\Tracking']);
for i=1:size(frameInfo,1)
    bioTreeTimer=etime(frameInfo(i,1:6),bioTree{1}.bioTreeInitialTime)/60;
    diffTime=abs(bioTree{1}.bioTreeTimer-bioTreeTimer);
    [~,index]=find(diffTime==min(diffTime));
    %     bioTreeFrame(i)=index(1);
    bioTreeFrame(i)=index(2);
end
[~,checkFrame]=find(diff(bioTreeFrame)==0);
bioTreeFrame(checkFrame+1)=0;
for iPic=1:numel(bioTreeFrame)
    if bioTreeFrame(iPic)==0
        continue
    end
    %         image=load([dirFluoFile,'\image',fluoChannel{iChannel},num2str(iPic,'%05.f'),'.mat']);
    %     tempimage=load([dirFile,'\Tracking,'\image','Tacking',num2str(iPic,'%05.f'),'.mat']);
    %     image=tempimage.imageTracking;
    
    for iframe=1:bioTreeFrame(iPic)
        for iRoot=1:numel(bioTree{iframe}.root)
            traceInfo=bioTree{iframe}.root{iRoot}.traceInfo;
            bioTree{iframe}.root{iRoot}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
            if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
                tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
            end
        end
        for iNode=1:numel(bioTree{iframe}.node)
            for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
                traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo;
                bioTree{iframe}.node{iNode}.Out{iOut}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
                if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
                    tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
                end
            end
        end
    end
end

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
n=0;
while n==0
    try
        [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
        branchList(:,5)=0;
        n=1;
    catch err
        bioTree=bioTreeCut(bioTree,numel(bioTree)-2);
    end
end
wrongBranch=[];
for iCoreBranch=1:size(branchList(branchList(:,3)==1,:),1)
    iBranchInfo=branchList(iCoreBranch,:);
    allNode=bioTree{iBranchInfo(1)}.node{iBranchInfo(2)}.allNode;
    [~,order]=sort(allNode(:,1));
    allNode=allNode(order,:);
    branchList(iCoreBranch,5)=numel(bioTree);
    for iNode=1:size(allNode,1)
        nodeInfo=allNode(iNode,:);
        isDivision=isDivsionNode(bioTree,nodeInfo);
        if isDivision==0;
            branchList(iCoreBranch,5)=nodeInfo(1)-1;
            if branchList(iCoreBranch,5)<branchList(iCoreBranch,1)
                wrongBranch=[wrongBranch;iCoreBranch];
            end
            break
        end
    end
end
branchList(wrongBranch,:)=[];
emptyNum=[];
for iBranch=1:size(branchList,1)
    if branchList(iBranch,5)==0
        if branchList(iBranch,1)>numel(bioTree)
            emptyNum=[emptyNum;iBranch];
        else
            branchList(iBranch,5)=numel(bioTree);
        end
    end
end
branchList(emptyNum,:)=[];
bioTree{1}.branchList=branchList;
end
function isDivision=isDivsionNode(bioTree,nodeInfo)
isDivision=0;
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2
    %     frameOut=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out;
    %     bacteriaSize=[frameOut{1}.traceInfo.measurment{1}.FilledArea,frameOut{2}.traceInfo.measurment{1}.FilledArea];
    %     if ~(max(bacteriaSize)-min(bacteriaSize)>150 || max(bacteriaSize)>min(bacteriaSize)*2 )
    isDivision=1;
    %     end
end
end

%% fluoImage processing
function [maskImage,imageStack]=fluoImageProcessing(imageStack)
% minIntensity=300;
% imageStack=imageStack-200;
for iStack=1:size(imageStack,3)
    image=imageStack(:,:,iStack);
    pixelInfo=image(:);
    pixelInfo=sort(pixelInfo);
    backGround=mean(pixelInfo(1:end/3));
    maxPixel=max(pixelInfo);
    image=image-backGround;
    edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
    gaussianFilter1=fspecial('gaussian',[10, 10],3);
    gaussianFilter2=fspecial('gaussian',[10, 10],2);
    image=uint8(double(image)/double(maxPixel)*255);
    image=imfilter(image,gaussianFilter1); % use Guassian blur filter process
    image=imfilter(image,edgeFilter); %use edgeFilter process
    image=imfilter(image,gaussianFilter2);
    thre=5;
    image=im2bw(image,thre/255);
    image=imclearborder(image);
    image=bwareaopen(image,100,4);
    image=logical(image);
    image=imfill(image,'holes');
    maskImage(:,:,iStack)=image;
    %     maskImage(:,:,iStack)=bwmorph(image,'close');
end
end


function [f1] = plotLinkTree(linkMatrix,centroidInfo,leafYCoo,color,type,cutNum)
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
close all
f1=figure;
set(f1,'visible','off')
axes1 = axes('Parent',f1,...
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
c1(1)=c1(1)/60;
c2(1)=c2(1)/60;
cutNum=cutNum/60;
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
    x=c1(1):0.001:cutNum;
else
    x=c1(1):0.001:c2(1);
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
% function bacNum=getBasicFigure(maskImages,dirFile,mip,fieldNum)
function bacNum=getBasicFigure(maskImages,dirFile,fieldNum)
fluoChannel{1,1}='CyOFP';
fluoChannel{2,1}='sfGFP';
fluoChannel{3,1}='mScarletI';
fluoChannel{4,1}='TDsmURFP';
fluoChannel{5,1}='CyPet';
fluoChannel{6,1}='Venus';
fluoChannel{7,1}='mAmetrine';
if isempty(maskImages)
    for iChannel=1:numel(fluoChannel)
        result=[];
        try
            load([dirFile,'\',fluoChannel{iChannel},'new\','frameInfo.mat'])
        catch err
            continue
        end
        switch iChannel
            case 1
                bacInfo{1}.meanCyOFP=1;
            case 2
                bacInfo{1}.meanGFP=1;
            case 3
                bacInfo{1}.meanmScalet=1;
            case 4
                bacInfo{1}.meanRFP=1;
            case 5
                bacInfo{1}.meanCyPet=1;
            case 6
                bacInfo{1}.meanVenus=1;
            case 7
                bacInfo{1}.meanmAmetrine=1;
        end
    end
    bacInfo{1}.majorLength=1;
    bacInfo{1}.minorLength=1;
    bacInfo{1}.time=1;
    %     if ~isempty(bacInfo)
    %         bacInfo{1}.tag=mip.feildTag.tag{fieldNum};
    %         bacInfo{1}.tagValue=mip.feildTag.tagValue(fieldNum);
    %     end
    save([dirFile,'\bacInfo&tree'],'bacInfo');
    bacNum=[];
    return
end
load([dirFile,'\Tracking\frameInfo.mat']);
frameInfo=frameInfo(1:size(maskImages,3)/2,:);%ys
timeBegin=frameInfo(1,1:6);
for iTime=1:size(frameInfo,1)
    timeAll(iTime)=etime(frameInfo(iTime,1:6),timeBegin)/60;
end
maskImages=maskImages(:,:,2:2:end);
dirMeanIntensityFile=[dirFile,'\meanIntensity'];
mkdir(dirMeanIntensityFile)

% plot bacNum vs time
for iImage=1:size(maskImages,3)
    image=maskImages(:,:,iImage);
    cc=bwconncomp(image);
    bacNum(iImage)=cc.NumObjects;
    cc=regionprops(image,'MajorAxisLength','MinorAxisLength');
    majorLength=[];
    minorLength=[];
    for iSmallBac=1:numel(cc)
        majorLength=[majorLength;cc(iSmallBac).MajorAxisLength];
        minorLength=[minorLength;cc(iSmallBac).MinorAxisLength];
    end
    bacInfo{iImage}.majorLength=majorLength;
    bacInfo{iImage}.minorLength=minorLength;
    bacInfo{iImage}.time=timeAll(iImage);
end
if ~isempty(maskImages)
    plot(timeAll/60,bacNum)
    saveas(gcf,[dirMeanIntensityFile,'\bacNum','.tif'])
    saveas(gcf,[dirMeanIntensityFile,'\bacNum','.fig'])
    close all
end

% plot meanIntensity vs time
for iChannel=1:numel(fluoChannel)
    bioTreeFrame=[];
    result=[];
    try
        load([dirFile,'\',fluoChannel{iChannel},'new\','frameInfo.mat'])
    catch err
        continue
    end
    for i=1:size(frameInfo,1)
        bioTreeTimer=etime(frameInfo(i,1:6),timeBegin)/60;
        diffTime=abs(timeAll-bioTreeTimer);
        [~,index]=find(diffTime==min(diffTime));
        bioTreeFrame(i)=index(1);
    end
    for iPic=1:numel(bioTreeFrame)
        if bioTreeFrame(iPic)==0
            continue
        end
        %         image=load([dirFile,'\',fluoChannel{iChannel},'new','\image',fluoChannel{iChannel},num2str((iPic),'%05.f'),'.mat']);
        %         image=image.imageI;
        image=load([dirFile,'\',fluoChannel{iChannel},'new','\image','Tacking',num2str(iPic,'%05.f'),'.mat']);
        image=image.imageTracking;
        maskI=maskImages(:,:,bioTreeFrame(iPic));
        cc=regionprops(maskI,image,'MeanIntensity');
        meanIntensity=[];
        for iSmallBac=1:numel(cc)
            meanIntensity=[meanIntensity;cc(iSmallBac).MeanIntensity];
        end
        switch iChannel
            case 1
                bacInfo{bioTreeFrame(iPic)}.meanCyOFP=meanIntensity;
            case 2
                bacInfo{bioTreeFrame(iPic)}.meanGFP=meanIntensity;
            case 3
                bacInfo{bioTreeFrame(iPic)}.meanmScalet=meanIntensity;
            case 4
                bacInfo{bioTreeFrame(iPic)}.meanRFP=meanIntensity;
            case 5
                bacInfo{bioTreeFrame(iPic)}.meanCyPet=meanIntensity;
            case 6
                bacInfo{bioTreeFrame(iPic)}.meanVenus=meanIntensity;
            case 7
                bacInfo{bioTreeFrame(iPic)}.meanmAmetrine=meanIntensity;
        end
        result=[result;timeAll(bioTreeFrame(iPic))/60,mean(image(maskI))];
    end
    plot(result(:,1),result(:,2))
    %     if ~isempty(bacInfo)
    %         bacInfo{1}.tag=mip.feildTag.tag{fieldNum};
    %         bacInfo{1}.tagValue=mip.feildTag.tagValue(fieldNum);
    %     end
    save([dirMeanIntensityFile,'\',fluoChannel{iChannel}],'result')
    saveas(gcf,[dirMeanIntensityFile,'\',fluoChannel{iChannel},'.tif'])
    saveas(gcf,[dirMeanIntensityFile,'\',fluoChannel{iChannel},'.fig'])
    close all
end
fluoChannel{8,1}='RedMask';
fluoChannel{9,1}='BlueMask';
fluoChannel{10,1}='GreenMask';
for iChannel=8:10
    bioTreeFrame=[];
    result=[];
    try
        load([dirFile,'\',fluoChannel{iChannel},'\','frameInfo.mat'])
    catch err
        continue
    end
    for i=1:size(frameInfo,1)-1
        deltaTime(i)=abs(etime(frameInfo(i+1,1:6),frameInfo(i,1:6)));
    end
    deltaTime(end+1)=deltaTime(end);
    frameInfo(:,15)=frameInfo(:,15).*frameInfo(:,9)./deltaTime'/1000;
    for i=1:size(frameInfo,1)
        bioTreeTimer=etime(frameInfo(i,1:6),timeBegin)/60;
        diffTime=abs(timeAll-bioTreeTimer);
        [~,index]=find(diffTime==min(diffTime));
        bioTreeFrame(i)=index(1);
    end
    for iPic=1:numel(bioTreeFrame)
        if bioTreeFrame(iPic)==0
            continue
        end
        image=import_tiff_stack([dirFile,'\',fluoChannel{iChannel}(1:end-4),'Control\image',fluoChannel{iChannel}(1:end-4),'Control',num2str(1,'%05.f'),'.tif']);
        maskI=maskImages(:,:,bioTreeFrame(iPic));
        cc=regionprops(maskI,image,'MeanIntensity');
        meanIntensity=[];
        for iSmallBac=1:numel(cc)
            meanIntensity=[meanIntensity;cc(iSmallBac).MeanIntensity];
        end
        switch iChannel
            case 8
                bacInfo{bioTreeFrame(iPic)}.redControl=meanIntensity*16/frameInfo(iPic,14)*frameInfo(iPic,15);
            case 9
                bacInfo{bioTreeFrame(iPic)}.blueControl=meanIntensity*16/frameInfo(iPic,14)*frameInfo(iPic,15);
            case 10
                bacInfo{bioTreeFrame(iPic)}.greenControl=meanIntensity*16/frameInfo(iPic,14)*frameInfo(iPic,15);
        end
    end
end
save([dirFile,'\bacInfo&tree'],'bacInfo');
end
function [imageStack,bestPositionAccumulation] = fluoImageCorrection(imageStack)
% 针对荧光图像的时间序列的xy矫正
% 输入变量为二值图像，核心矫正函数已经优化过
imageSize=size(imageStack);
if min(imageSize(1:2))>=1500
    imageStackNew=imageStack(1024-500:1024+500,1024-500:1024+500,:);
else
    imageStackNew=imageStack;
end
imageStackNew1(:,:,2:size(imageStackNew,3)+1)=imageStackNew;
parfor i=2:size(imageStack,3)
    bestPosition(i,:)=caculateCrossCorrelationForImage(imageStackNew1(:,:,i),imageStackNew(:,:,i),15);
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

%% bioTreeTwoPointAging
function bioTree=bioTreeTwoPointAging(bioTree)
bioTree=bioTreeTwoPointTracking(bioTree);
bioTree=bioTreePointMatch(bioTree);
end
function bioTree=bioTreeTwoPointTracking(bioTree)
bioTree=getTwoPointPosition(bioTree);
end
function bioTree=getTwoPointPosition(bioTree)
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if ~isempty(bioTree{iframe}.root{iroot}.traceInfo.measurment)
                for iMeasurment=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment{1},1)
                    centroid=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).Centroid;
                    orientation=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).Orientation;
                    majorAxisLength=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).MajorAxisLength;
                    eccentricity=bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).Eccentricity;
                    [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity);
                    bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).p1Position=p1Position;
                    bioTree{iframe}.root{iroot}.traceInfo.measurment{1}(iMeasurment).p2Position=p2Position;
                end
                for iTrace=2:size(bioTree{iframe}.root{iroot}.traceInfo.measurment,2)
                    for iMeasurment=1:size(bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace},1)
                        centroid_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).Centroid;
                        orientation_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).Orientation;
                        majorAxisLength_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).MajorAxisLength;
                        eccentricity_next=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).Eccentricity;
                        [p1Position_next,p2Position_next]=twoPointPosition(centroid_next,orientation_next,majorAxisLength_next,eccentricity_next);
                        p1Position_pre= bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}(1).p1Position;
                        p2Position_pre= bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}(1).p2Position;
                        [p1Position,p2Position]=twoPointConnect(p1Position_pre,p1Position_next,p2Position_pre,p2Position_next);
                        bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).p1Position=p1Position;
                        bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}(iMeasurment).p2Position=p2Position;
                    end
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                if ~isempty(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment)
                    for iMeasurment=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1},1)
                        centroid= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).Centroid;
                        orientation= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).Orientation;
                        majorAxisLength= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).MajorAxisLength;
                        eccentricity= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).Eccentricity;
                        [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity);
                        bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).p1Position=p1Position;
                        bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}(iMeasurment).p2Position=p2Position;
                    end
                    for iTrace=2:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment,2)
                        for iMeasurment=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace},1)
                            centroid_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).Centroid;
                            orientation_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).Orientation;
                            majorAxisLength_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).MajorAxisLength;
                            eccentricity_next= bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).Eccentricity;
                            [p1Position_next,p2Position_next]=twoPointPosition(centroid_next,orientation_next,majorAxisLength_next,eccentricity_next);
                            p1Position_pre=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace-1}(1).p1Position;
                            p2Position_pre=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace-1}(1).p2Position;
                            [p1Position,p2Position]=twoPointConnect(p1Position_pre,p1Position_next,p2Position_pre,p2Position_next);
                            bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).p1Position=p1Position;
                            bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(iMeasurment).p2Position=p2Position;
                        end
                    end
                end
            end
        end
    end
end
end
function [p1Position,p2Position]=twoPointPosition(centroid,orientation,majorAxisLength,eccentricity)
O_sign=sign(orientation);
if orientation==0
    O_sign=1;
end
O_cos=cos(abs(pi*orientation/180));
O_sin=sin(abs(pi*orientation/180));
Axis_length=0.5*majorAxisLength;
x=centroid(1);
y=centroid(2);
p1_x=x-O_sign*Axis_length*O_cos*eccentricity;
p1_y=y+Axis_length*O_sin*eccentricity;
p2_x=x+O_sign*Axis_length*O_cos*eccentricity;
p2_y=y-Axis_length*O_sin*eccentricity;
p1Position=[p1_x,p1_y];
p2Position=[p2_x,p2_y];
end
function  [p1Position,p2Position]=twoPointConnect(p1Position_pre,p1Position_next,p2Position_pre,p2Position_next)
vectorP1P2_pre=p1Position_pre-p2Position_pre;
vectorP1P2_next=p1Position_next-p2Position_next;
if sum(vectorP1P2_pre.*vectorP1P2_next)>0
    p1Position=p1Position_next;
    p2Position=p2Position_next;
    return;
end
if sum(vectorP1P2_pre.*vectorP1P2_next)<0
    p1Position=p2Position_next;
    p2Position=p1Position_next;
    return;
end
if sum(vectorP1P2_pre.*vectorP1P2_next)==0
    D1=sum((p1Position_pre-p1Position_next).^2);
    D2=sum((p1Position_pre-p2Position_next).^2);
    if D1<=D2
        p1Position=p1Position_next;
        p2Position=p2Position_next;
        return;
    else
        p1Position=p2Position_next;
        p2Position=p1Position_next;
        return;
    end
end
end
function bioTree=bioTreePointMatch(bioTree)
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        bioTree{iframe}.root{iRoot}.isOld1=1;
        bioTree{iframe}.root{iRoot}.isOld2=1;
    end
end
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            bioTree{iframe}.node{iNode}.Out{iOut}.isOld1=-inf;
            bioTree{iframe}.node{iNode}.Out{iOut}.isOld2=-inf;
        end
    end
end
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        if size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)~=1
            preIsNode=bioTree{iframe}.node{iNode}.In{1}.isNode;
            if preIsNode==0
                outInfo=[];
                preRoot=bioTree{iframe}.node{iNode}.In{1}.rootInfo;
                if iframe~=preRoot(1)
                    inInfo=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}.p1Position,bioTree{preRoot(1)}.root{preRoot(2)}.isOld1;...
                        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}.p2Position,bioTree{preRoot(1)}.root{preRoot(2)}.isOld2];
                else
                    inInfo=[0,0,1;0,0,1];
                end
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    outInfo(iOut*2-1:2*iOut,:)=[bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p1Position;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p2Position];
                end
                distMatrix=pdist2(inInfo(:,1:2),outInfo);
                orderNum=1:size(outInfo,1);
                matchNum=orderNum(distMatrix(1,:)==min(distMatrix(1,:)));
                matchNum=matchNum(1);
                outNum=fix((matchNum+1)/2);
                leaveNum=matchNum-(outNum-1)*2;
                if leaveNum==1
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preRoot(1)}.root{preRoot(2)}.isOld1+1;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                end
                if leaveNum==2
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preRoot(1)}.root{preRoot(2)}.isOld1+1;
                end
                distMatrix(:,matchNum)=50;
                matchNum=orderNum(distMatrix(2,:)==min(distMatrix(2,:)));
                matchNum=matchNum(1);
                outNum=fix((matchNum+1)/2);
                leaveNum=matchNum-(outNum-1)*2;
                if leaveNum==1
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preRoot(1)}.root{preRoot(2)}.isOld2+1;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                end
                if leaveNum==2
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                    bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preRoot(1)}.root{preRoot(2)}.isOld2+1;
                end
            end
            if preIsNode==1
                outInfo=[];
                preNode=bioTree{iframe}.node{iNode}.In{1}.nodeInfo;
                inInfo=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}.p1Position;bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}.p2Position];
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    outInfo(iOut*2-1:2*iOut,:)=[bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p1Position;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{1}.p2Position];
                end
                distMatrix=pdist2(inInfo(:,1:2),outInfo);
                orderNum=1:size(outInfo,1);
                if ~isempty(bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld1)
                    matchNum=orderNum(distMatrix(1,:)==min(distMatrix(1,:)));
                    matchNum=matchNum(1);
                    outNum=fix((matchNum+1)/2);
                    leaveNum=matchNum-(outNum-1)*2;
                    if leaveNum==1
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld1+1;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                    end
                    if leaveNum==2
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld1+1;
                    end
                    distMatrix(:,matchNum)=50;
                    matchNum=orderNum(distMatrix(2,:)==min(distMatrix(2,:)));
                    matchNum=matchNum(1);
                    outNum=fix((matchNum+1)/2);
                    leaveNum=matchNum-(outNum-1)*2;
                    if leaveNum==1
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld2+1;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=0;
                    end
                    if leaveNum==2
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld1=0;
                        bioTree{iframe}.node{iNode}.Out{outNum}.isOld2=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.isOld2+1;
                    end
                end
            end
        end
    end
end
end

% montageImageInitialized
function montageImageInitialized(dirFile,mip)
% montageXY=input('montage X num,Y num amd gapValue [n,m,gap] = ');
montageXY=[12,12,-0.3];
picSizeX=round(2048*abs(montageXY(3))/2);
picSizeY=2048-picSizeX;

fluo_calibration=mip.Calib.fluo_protein_calibration1_6;
fluoChannel=mip.Calib.fluo_protein_info;
fieldNameList=dir(dirFile);

iField=3;
fluoHere=[];
n=0;
dirFieldFile=[dirFile,'\',fieldNameList(iField+2).name];
for iFluo=1:numel(fluoChannel)
    try load([dirFieldFile,'\',fluoChannel{iFluo},'\frameInfo.mat']);
        n=n+1;
        fluoChannelChoose(n)=iFluo;
        fluoHere=[fluoHere,fluoChannel(iFluo)];
        mkdir([dirFieldFile,'\',fluoChannel{iFluo},'new']);
        copyfile([dirFieldFile,'\',fluoChannel{iFluo},'\frameInfo.mat'],[dirFieldFile,'\',fluoChannel{iFluo},'new\frameInfo.mat'])
        switch iFluo
            case 1
                backGround(:,:,n)=mip.Calib.LumencorBlueField;
            case 2
                backGround(:,:,n)=mip.Calib.LumencorTealField;
            case 3
                backGround(:,:,n)=mip.Calib.LumencorVioletField;
            case 4
                backGround(:,:,n)=mip.Calib.LumencorCyanField;
            case 5
                backGround(:,:,n)=mip.Calib.LumencorCyanField;
            case 6
                backGround(:,:,n)=mip.Calib.LumencorGreenField;
            case 7
                backGround(:,:,n)=mip.Calib.LumencorRedField;
        end
    end
end
for iFluo=1:numel(fluoHere)
    fieldNum=0;
    imageAll=[];
    for iField=1:numel(fieldNameList)-2
        dirFieldFile=[dirFile,'\',fieldNameList(iField+2).name];
        if numel(fieldNameList(iField+2).name)>=9 && strcmp(fieldNameList(iField+2).name(end-8:end-4),'feild')
            fieldNum=fieldNum+1;
            try
                for i=1:numel(dir([dirFieldFile,'\',fluoHere{1}]))-3
                    load([dirFieldFile,'\',fluoHere{iFluo},'\frameInfo.mat']);
                    image=import_tiff_stack([dirFieldFile,'\',fluoHere{iFluo},'\image',fluoHere{iFluo},num2str(i,'%05.f'),'.tif']);
                    if iFluo==1
                        image=double(image)-150;
                    else
                        image=double(image)-150;
                    end
                    image=image./backGround(:,:,iFluo);
                    image=image(picSizeX:picSizeY,picSizeX:picSizeY);
                    imageAll(:,:,fieldNum)=uint16(image);
                    imageAll=uint16(imageAll);
                end
            end
        end
    end
    imageSize=picSizeY-picSizeX+1;
    image_final=uint16(zeros(imageSize*montageXY(1),imageSize*montageXY(2)));
    indexY=0:imageSize:imageSize*montageXY(1);
    indexX=0:imageSize:imageSize*montageXY(2);
    for i=1:montageXY(1)
        for j=1:montageXY(2)
            montageNum=(i-1)*montageXY(2)+j;
            image_final(indexY(i)+1:indexY(i+1),indexX(j)+1:indexX(j+1))=imrotate(imageAll(:,:,montageNum),90);
        end
    end
    imwrite(image_final,[dirFile,'\feild0001\',fluoHere{iFluo},'.tif']);
    save([dirFile,'\feild0001\',fluoHere{iFluo}],'image_final')
end
end
function fluoImageInitialized(dirFile,mip)
if strcmp(mip.Calib.Olympus,'1x')
    fluo_calibration=mip.Calib.fluo_protein_calibration1;
else
    fluo_calibration=mip.Calib.fluo_protein_calibration1_6;
end
fluoChannel=mip.Calib.fluo_protein_info;
fieldNameList=dir(dirFile);
sita=144.9;
load([dirFile,'\LumencorField.mat']);
for iField=1:numel(fieldNameList)-2
    fieldNum=0;
    dirFieldFile=[dirFile,'\',fieldNameList(iField+2).name];
    if numel(fieldNameList(iField+2).name)>=9 && strcmp(fieldNameList(iField+2).name(end-8:end-4),'feild')
        fieldNum=fieldNum+1;
        if fieldNum==1
            fluoHere=[];
            n=0;
            for iFluo=1:numel(fluoChannel)
                try load([dirFieldFile,'\',fluoChannel{iFluo},'\frameInfo.mat']);
                    n=n+1;
                    fluoChannelChoose(n)=iFluo;
                    fluoHere=[fluoHere,fluoChannel(iFluo)];
                    mkdir([dirFieldFile,'\',fluoChannel{iFluo},'new']);
                    copyfile([dirFieldFile,'\',fluoChannel{iFluo},'\frameInfo.mat'],[dirFieldFile,'\',fluoChannel{iFluo},'new\frameInfo.mat'])
                    switch iFluo
                        case 1
                            backGround(:,:,n)=LumencorField.BlueField;
                        case 2
                            backGround(:,:,n)=LumencorField.TealField;
                        case 3
                            backGround(:,:,n)=LumencorField.VioletField;
                        case 4
                            backGround(:,:,n)=LumencorField.CyanField;
                        case 5
                            backGround(:,:,n)=LumencorField.CyanField;
                        case 6
                            backGround(:,:,n)=LumencorField.GreenField;
                        case 7
                            backGround(:,:,n)=LumencorField.RedField;
                    end
                end
            end
        end
        trackingFluoChannelReGet(dirFieldFile,fluoHere);
        for i=1:numel(dir([dirFieldFile,'\',fluoHere{1}]))-3
            try
                for iFluo=1:numel(fluoHere)
                    load([dirFieldFile,'\',fluoHere{iFluo},'\frameInfo.mat']);
                    image=import_tiff_stack([dirFieldFile,'\',fluoHere{iFluo},'\image',fluoHere{iFluo},num2str(i,'%05.f'),'.tif']);
                    image=double(image)-100;
                    image=image./backGround(:,:,iFluo);
                    image=image/frameInfo(i,14)*frameInfo(i,13)/frameInfo(i,9);
                    image=image(:);
                    newMatrix(iFluo,:)=double(image);
                end
                image=[];
                caliMatrix=fluo_calibration(fluoChannelChoose,fluoChannelChoose);
                newMatrix=inv(caliMatrix)'*newMatrix;
                for iFluo=1:numel(fluoHere)
                    imageI=reshape(newMatrix(iFluo,:),2048,2048)*sita/1000;
                    save([dirFieldFile,'\',fluoHere{iFluo},'new\image',fluoHere{iFluo},num2str(i,'%05.f'),'.mat'],'imageI')
                end
            end
        end
    end
end
end
%% getTreegeneration maker
function tree=treeSizeAndGenerationMarker(tree,calTime)
fluoProtein{1}='CyOFP';
fluoProtein{2}='GFP';
fluoProtein{3}='mScalet';
fluoProtein{4}='RFP';
fluoProtein{5,1}='CyPet';
fluoProtein{6,1}='Venus';
fluoProtein{7,1}='mAmetrine';
for iTree=1:numel(tree)
    try
        smallTree=tree(iTree);
        linkInfo=smallTree.linkInfo;
        tree(iTree).linkInfo(1).linkGeneration=1;
        linkGeneration{1}=1;
        linkGeneration{2}=1;
        for iLink=2:numel(linkInfo)
            targetInfo=smallTree.linkColumn(smallTree.linkRow==smallTree.linkRow(iLink));
            linkTagetNum=numel(targetInfo);
            index=find(targetInfo==smallTree.linkColumn(iLink));
            linkGeneration{smallTree.linkColumn(iLink)}=[linkGeneration{smallTree.linkRow(iLink)};index];
            linkInfo(iLink).linkGeneration=linkGeneration{smallTree.linkColumn(iLink)};
            tree(iTree).linkInfo(iLink).linkGeneration=numel(linkInfo(iLink).linkGeneration);
        end
    catch err
        linkInfo=smallTree.linkInfo;
        for iLink=1:numel(linkInfo)
            tree(iTree).linkInfo(iLink).linkGeneration=NaN;
        end
    end
    for iProtein=1:numel(fluoProtein)
        try
            bacInfo=getLeafInfo(linkInfo,tree(iTree).timer(end),fluoProtein{iProtein});
        catch err
            bacInfo=[];
        end
        if isempty(bacInfo)
            result=[];
        else
            try
                result=calculateCorrelation(bacInfo);
            catch err
                result=[];
            end
        end
        switch iProtein
            case 1
                tree(iTree).CyOFPcorr=result;
            case 2
                tree(iTree).GFPcorr=result;
            case 3
                tree(iTree).mScaletcorr=result;
            case 4
                tree(iTree).RFPcorr=result;
            case 5
                tree(iTree).CyPetcorr=result;
            case 6
                tree(iTree).Venuscorr=result;
            case 7
                tree(iTree).mAmetrinecorr=result;
        end
    end
end
end
function bacInfo=getLeafInfo(linkInfo,givenTime,fluoProtein)
n=0;
bacInfo=[];
for i=2:numel(linkInfo);
    linkTimer=linkInfo(i).linkTimer;
    diffTime=abs(linkTimer-givenTime);
    indexDiff=find(diffTime<=2);
    if isempty(indexDiff)
        continue
    end
    switch fluoProtein
        case 'CyOFP'
            result=linkInfo(i).CyOFP(indexDiff(2));
        case 'GFP'
            result=linkInfo(i).GFP(indexDiff(2));
        case 'mScalet'
            result=linkInfo(i).mScalet(indexDiff(2));
        case 'RFP'
            result=linkInfo(i).RFP(indexDiff(2));
        case 'CyPet'
            result=linkInfo(i).CyPet(indexDiff(2));
        case 'Venus'
            result=linkInfo(i).Venus(indexDiff(2));
        case 'mAmetrine'
            result=linkInfo(i).mAmetrine(indexDiff(2));
    end
    if result~=0;
        n=n+1;
        bacInfo{n}.generateMark=linkInfo(i).linkGeneration;
        bacInfo{n}.meanIntensity=result;
    end
end
end
function result=calculateCorrelation(bacInfo)
allTable=[];
for i=1:numel(bacInfo)
    bacI.generateMark=bacInfo{i}.generateMark;
    bacI.meanIntensity=bacInfo{i}.meanIntensity;
    for j=1:numel(bacInfo)
        if j<=i
            bacJ.generateMark=bacInfo{j}.generateMark;
            bacJ.meanIntensity=bacInfo{j}.meanIntensity;
            dis=calDis(bacI.generateMark,bacJ.generateMark);
            allTable=[allTable;dis,bacI.meanIntensity,bacJ.meanIntensity];
        end
    end
end
for i=0:max(allTable(:,1))
    corrResult=corrcoef(allTable(allTable(:,1)==i,2),allTable(allTable(:,1)==i,3));
    result(i+1,1)=corrResult(1,2);
end
end
function dis=calDis(mark1,mark2)
% 分裂后两个子代之间的距离为1
sameGeneration=0;
for i=1:min(numel(mark1),numel(mark2))
    if mark1(i)==mark2(i)
        sameGeneration=i;
    else
        break
    end
end
dis=numel(mark1)-sameGeneration+numel(mark2)-sameGeneration-1;
if dis==-1
    dis=0;
end
end

function trackingFluoChannelReGet(dirFieldFile,fluoHere)
for i=1:numel(fluoHere)
    fileNum(i)=numel(dir([dirFieldFile,'\',fluoHere{i}]));
end
if max(fileNum)==min(fileNum)
    return
end
trackingNum=find(fileNum==max(fileNum));
fluoNum=find(fileNum==min(fileNum));
fluoFrameInfo=load([dirFieldFile,'\',fluoHere{fluoNum(1)},'\frameInfo.mat']);
fluoFrameInfo=fluoFrameInfo.frameInfo;
trackingFrameInfo=load([dirFieldFile,'\',fluoHere{trackingNum(1)},'\frameInfo.mat']);
trackingFrameInfo=trackingFrameInfo.frameInfo;
movefile([dirFieldFile,'\',fluoHere{trackingNum(1)}],[dirFieldFile,'\',fluoHere{trackingNum(1)},'Pre']);
mkdir([dirFieldFile,'\',fluoHere{trackingNum(1)}]);
for i=1:size(fluoFrameInfo,1)
    timeI=fluoFrameInfo(i,1:6);
    trackingTime=trackingFrameInfo(:,1:6);
    for j=1:size(trackingTime,1)
        deltaTime(j)=etime(trackingTime(j,:),timeI);
    end
    focusTime=find(abs(deltaTime)==min(abs(deltaTime)));
    copyfile([dirFieldFile,'\',fluoHere{trackingNum(1)},'Pre\image',fluoHere{trackingNum(1)},num2str(focusTime,'%05.f'),'.tif'],[dirFieldFile,'\',fluoHere{trackingNum(1)},'\image',fluoHere{trackingNum(1)},num2str(i,'%05.f'),'.tif']);
    frameInfo(i,:)=trackingFrameInfo(focusTime,:);
end
save([dirFieldFile,'\',fluoHere{trackingNum(1)},'\frameInfo.mat'],'frameInfo');
end
function createAllData(dirFile)
%GETALLDATAPREPARED Summary of this function goes here
%   Detailed explanation goes here
allData.maskResult=[];
allData.maskResultTag=[];
allData.maskResultTagValue=[];
allData.trackingResult=[];
allData.trackingResultTag=[];
allData.trackingResultTagValue=[];
allData.correlationDataCYOFP=[];
allData.correlationDataGFP=[];
allData.correlationDatamScalet=[];
allData.correlationDataRFP=[];
allData.correlationDataCyPet=[];
allData.correlationDataVenus=[];
allData.correlationDatamAmetrine=[];
treeSizeAll=[];
% allData里面含两类数据，一个是从mask得到的bacInfo
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP,第八列为CyPet,第九列为Venus,第十列为Ametrine一一对应，没有数据的地方用NaN表示
% 第一列为time(min),第二列为majorLength,第三列为minorLength,第四列为CyOFP,第五列为GFP,第六列为mScalet,第七列为RFP，一一对应，没有数据的地方用NaN表示，第八列为Orientation,第九列为FilledArea,
% 第十列为growthRate,第十一列为fpProduceGFP，第十二列为treeSize,第十三列为generation,第十四列为CyPet,第十五列为Venus,第十六列为Ametrine，第十七列为fpProduceCyOFP,第十八列为fpProducemScalet
% 第十九列为fpProduceRFP，第二十列为fpProduceCyPet，第二十一列为fpProduceVenus,第二十二列为fpProduceAmetrine,第二十三列为redControl,第二十四列为blueControl,第二十五列为greenControl
data1Size=13;
data2Size=25;
treeAll=[];
a=load([dirFile,'\bacInfo&tree.mat']);
bacInfo=a.bacInfo;
% scale=a.tree(1).scale;
scale=0.065;
for iTime=1:numel(bacInfo)
    maskResult=[];
    maskResultTag=[];
    maskResultTagValue=[];
    dataLength=numel(bacInfo{iTime}.majorLength);
    for iData=1:data1Size
        switch iData
            case 1
                maskResult(1:dataLength,iData)=bacInfo{iTime}.time;
            case 2
                maskResult(1:dataLength,iData)=bacInfo{iTime}.majorLength*scale;
            case 3
                maskResult(1:dataLength,iData)=bacInfo{iTime}.minorLength*scale;
            case 4
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanCyOFP;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 5
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanGFP;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 6
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanmScalet;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 7
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanRFP;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 8
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanCyPet;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 9
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.Venus;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 10
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.mAmetrine;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 11
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.redControl;
                catch err
                    maskResult(1:dataLength,iData)=0;
                end
            case 12
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.blueControl;
                catch err
                    maskResult(1:dataLength,iData)=0;
                end
            case 13
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.greenControl;
                catch err
                    maskResult(1:dataLength,iData)=0;
                end
        end
    end
    %     for iResult=1:size(maskResult,1)
    %         maskResultTag{iResult,1}=bacInfo{1}.tag;
    %         maskResultTagValue(iResult,1)=bacInfo{1}.tagValue;
    %     end
    %     allData.maskResult=[allData.maskResult;maskResult];
    %     allData.maskResultTagValue=[allData.maskResultTagValue;maskResultTagValue];
    %     allData.maskResultTag=[allData.maskResultTag;maskResultTag];
end
bacTree=a.tree;
treeAll=[treeAll,bacTree];
for iTime=1:numel(bacTree)
    leafNum=bacTree(iTime).leafNum;
    treeSizeAll=[treeSizeAll;ceil(log(leafNum)/log(2)+1)];
    if ~isempty(bacTree(iTime).CyOFPcorr)
        allData.correlationDataCYOFP=[allData.correlationDataCYOFP;cat(2,(0:numel(bacTree(iTime).CyOFPcorr)-1)',abs(bacTree(iTime).CyOFPcorr),ones(size(bacTree(iTime).CyOFPcorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    if ~isempty(bacTree(iTime).GFPcorr)
        allData.correlationDataGFP=[allData.correlationDataGFP;cat(2,(0:numel(bacTree(iTime).GFPcorr)-1)',abs(bacTree(iTime).GFPcorr),ones(size(bacTree(iTime).GFPcorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    if ~isempty(bacTree(iTime).mScaletcorr)
        allData.correlationDatamScalet=[allData.correlationDatamScalet;cat(2,(0:numel(bacTree(iTime).mScaletcorr)-1)',abs(bacTree(iTime).mScaletcorr),ones(size(bacTree(iTime).mScaletcorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    if ~isempty(bacTree(iTime).RFPcorr)
        allData.correlationDataRFP=[allData.correlationDataRFP;cat(2,(0:numel(bacTree(iTime).RFPcorr)-1)',abs(bacTree(iTime).RFPcorr),ones(size(bacTree(iTime).RFPcorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    if ~isempty(bacTree(iTime).CyPetcorr)
        allData.correlationDataCyPet=[allData.correlationDataCyPet;cat(2,(0:numel(bacTree(iTime).CyPetcorr)-1)',abs(bacTree(iTime).CyPetcorr),ones(size(bacTree(iTime).CyPetcorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    if ~isempty(bacTree(iTime).Venuscorr)
        allData.correlationDataVenus=[allData.correlationDataVenus;cat(2,(0:numel(bacTree(iTime).Venuscorr)-1)',abs(bacTree(iTime).Venuscorr),ones(size(bacTree(iTime).Venuscorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    if ~isempty(bacTree(iTime).mAmetrinecorr)
        allData.correlationDatamAmetrine=[allData.correlationDatamAmetrine;cat(2,(0:numel(bacTree(iTime).mAmetrinecorr)-1)',abs(bacTree(iTime).mAmetrinecorr),ones(size(bacTree(iTime).mAmetrinecorr,1),1)*ceil(log(leafNum)/log(2)+1))];
    end
    for iLink=1:numel(bacTree(iTime).linkInfo)
        for iBac=1:numel(bacTree(iTime).linkInfo(iLink).linkTimer)/2
            trackingResult=[];
            trackingResultTag=[];
            trackingResultTagValue=[];
            for iData=1:data2Size
                switch iData
                    case 1
                        trackingResult(iData)=bacTree(iTime).linkInfo(iLink).linkTimer(2*iBac);
                    case 2
                        trackingResult(iData)=bacTree(iTime).linkInfo(iLink).measurment{2*iBac}(1).MajorAxisLength*scale;
                    case 3
                        trackingResult(iData)=bacTree(iTime).linkInfo(iLink).measurment{2*iBac}(1).MinorAxisLength*scale;
                    case 4
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).CyOFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 5
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).GFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 6
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).mScalet(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 7
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).RFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 8
                        trackingResult(iData)=bacTree(iTime).linkInfo(iLink).measurment{2*iBac}.Orientation;
                    case 9
                        trackingResult(iData)=bacTree(iTime).linkInfo(iLink).measurment{2*iBac}.FilledArea;
                    case 10
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).growthRate(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 11
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProduceGFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 12
                        trackingResult(iData)=ceil(log(leafNum)/log(2)+1);
                    case 13
                        trackingResult(iData)=bacTree(iTime).linkInfo(iLink).linkGeneration;
                    case 14
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).CyOFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 15
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).Venus(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 16
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).mAmetrine(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 17
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProduceCyOFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 18
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProducemScalet(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 19
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProduceRFP(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 20
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProduceCyPet(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 21
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProduceVenus(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 22
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).fpProducemAmetrine(2*iBac);
                        catch err
                            trackingResult(iData)=NaN;
                        end
                    case 23
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).redControl(2*iBac);
                        catch err
                            trackingResult(iData)=0;
                        end
                    case 24
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).blueControl(2*iBac);
                        catch err
                            trackingResult(iData)=0;
                        end
                    case 25
                        try
                            trackingResult(iData)=bacTree(iTime).linkInfo(iLink).greenControl(2*iBac);
                        catch err
                            trackingResult(iData)=0;
                        end
                end
            end
            %             for iResult=1:size(trackingResult,1)
            %                 trackingResultTag{iResult,1}=bacInfo{1}.tag;
            %                 trackingResultTagValue(iResult,1)=bacInfo{1}.tagValue;
            %             end
            %             allData.trackingResult=[allData.trackingResult;trackingResult];
            %             allData.trackingResultTagValue=[allData.trackingResultTagValue;trackingResultTagValue];
            %             allData.trackingResultTag=[allData.trackingResultTag;trackingResultTag];
        end
    end
end
allData.treeAll=treeAll;
allData.treeSize=treeSizeAll;
save([dirFile,'\allData.mat'],'allData')
end
function createAllData_brieflyMask(dirFile,mip)
%GETALLDATAPREPARED Summary of this function goes here
%   Detailed explanation goes here
allData.maskResult=[];
allData.maskResultTag=[];
allData.maskResultTagValue=[];
allData.trackingResult=[];
allData.trackingResultTag=[];
allData.trackingResultTagValue=[];
allData.correlationDataCYOFP=[];
allData.correlationDataGFP=[];
allData.correlationDatamScalet=[];
allData.correlationDataRFP=[];
allData.correlationDataCyPet=[];
allData.correlationDataVenus=[];
allData.correlationDatamAmetrine=[];
treeSizeAll=[];
data1Size=13;
data2Size=25;
treeAll=[];
a=load([dirFile,'\bacInfo&tree.mat']);
bacInfo=a.bacInfo;
scale=mip.Calib.scale;
for iTime=1:numel(bacInfo)
    maskResult=[];
    maskResultTag=[];
    maskResultTagValue=[];
    dataLength=numel(bacInfo{iTime}.majorLength);
    for iData=1:data1Size
        switch iData
            case 1
                maskResult(1:dataLength,iData)=bacInfo{iTime}.time;
            case 2
                maskResult(1:dataLength,iData)=bacInfo{iTime}.majorLength*scale;
            case 3
                maskResult(1:dataLength,iData)=bacInfo{iTime}.minorLength*scale;
            case 4
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanCyOFP;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 5
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanGFP;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 6
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanmScalet;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 7
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanRFP;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 8
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.meanCyPet;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 9
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.Venus;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 10
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.mAmetrine;
                catch err
                    maskResult(1:dataLength,iData)=NaN;
                end
            case 11
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.redControl;
                catch err
                    maskResult(1:dataLength,iData)=0;
                end
            case 12
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.blueControl;
                catch err
                    maskResult(1:dataLength,iData)=0;
                end
            case 13
                try
                    maskResult(1:dataLength,iData)=bacInfo{iTime}.greenControl;
                catch err
                    maskResult(1:dataLength,iData)=0;
                end
        end
    end
    for iResult=1:size(maskResult,1)
        maskResultTag{iResult,1}=bacInfo{1}.tag;
        maskResultTagValue(iResult,1)=bacInfo{1}.tagValue;
    end
    allData.maskResult=[allData.maskResult;maskResult];
    allData.maskResultTagValue=[allData.maskResultTagValue;maskResultTagValue];
    allData.maskResultTag=[allData.maskResultTag;maskResultTag];
end
allData.trackingResult=[];
allData.trackingResultTagValue=[];
allData.trackingResultTag=[];
allData.treeAll=[];
allData.treeSize=[];
save([dirFile,'\allData.mat'],'allData')
end
function bioTree=bioTreeDigHoleAndMeasure(bioTree,xSize,ySize)
parfor iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        bioTree{iframe}.root=DigRoot(bioTree{iframe}.root,xSize,ySize);
    end
    if ~isempty(bioTree{iframe}.node)
        bioTree{iframe}.node=DigNode(bioTree{iframe}.node,xSize,ySize);
    end
    if ~isempty(bioTree{iframe}.leavies)
        bioTree{iframe}.leavies=DigLeaf(bioTree{iframe}.leavies,xSize,ySize);
    end
end
end
function bioRoot=DigRoot(bioRoot,xSize,ySize)
imageSize=[xSize,ySize];
for iRoot=1:size(bioRoot,2)
    pixelIdxList=bioRoot{iRoot}.rootPixelDetail;
    [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
    bioRoot{iRoot}.rootPixelDetail=pixelIdxList;
    bioRoot{iRoot}.rootMeasurment=pos;
    bioRoot{iRoot}.traceInfo.measurment=[];
    for iTrace=1:size(bioRoot{iRoot}.traceInfo.pixelIdxList,2)
        pixelIdxList=bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace};
        [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
        bioRoot{iRoot}.traceInfo.measurment{iTrace}=pos;
        bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace}=pixelIdxList;
    end
end
end
function bioLeaf=DigLeaf(bioLeaf,xSize,ySize)
imageSize=[xSize,ySize];
try
    bioLeaf{iLeaf}.leaviesPixelDetail;
    for iLeaf=1:size(bioLeaf,2)
        pixelIdxList=bioLeaf{iLeaf}.leaviesPixelDetail;
        [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
        bioLeaf{iLeaf}.leaviesPixelDetail=pixelIdxList;
        bioLeaf{iLeaf}.leafMeasurment=pos;
        
    end
catch err
end
end
function bioNode=DigNode(bioNode,xSize,ySize)
imageSize=[xSize,ySize];
for iNode=1:size(bioNode,2)
    for iNodeOut=1:size(bioNode{iNode}.Out,2)
        bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment=[];
        for iNodeTrace=1:size(bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
            pixelIdxList=bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
            [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
            bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}=pos;
            bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}=pixelIdxList;
        end
    end
end
end
function [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize)
[xyMin,BWImage]=idx2Xy(pixelIdxList,imageSize);
pos=regionprops(BWImage,'FilledArea','Centroid','Eccentricity','MajorAxisLength','Orientation','MinorAxisLength');
% BWImage=imfill(BWImage,'holes');  % 原来是屏蔽的
for i=1:size(pos,1)
    pos(i).Centroid(1)=pos(i).Centroid(1)+xyMin(2)-1;
    pos(i).Centroid(2)=pos(i).Centroid(2)+xyMin(1)-1;
end
pixelIdxList=xy2Idx(xyMin,BWImage,imageSize);
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function bioTree=longTimeBioTreeAnalysis(bioTree)
bioTree=autoTidyBioTree(bioTree);
disp('Node Reduction')
nodeSortPre=[];
[nodeSort,~]=findNode(bioTree);
try
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
end
try
    bioTree=mergeNode(bioTree);
catch
    warning('Problem using the function mergeNode. Skipped.');
end
nodeSortPre=[];
try
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
end
disp('leafRootRefineBioTree')
try
    bioTree=leafRootRefineBioTree(bioTree,0);
catch 
     warning('Problem using the function leafRootRefineBioTree. Skipped.');
end
disp('bioTreeOptimize')
try
    bioTree=bioTreeOptimize(bioTree);
catch
    warning('Problem using the function bioTreeOptimize. Skipped.');
end
end
function bioTree=batchTreeTrackingJob(dirFile)
% dirTracking=[dirFile,'\Tracking'];
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

n=0;
attempTimes=0;
while n==0
    try
        bioTree=combineSmallTreetoBigTree(length(nameList)-2,dirResultSave);
        n=1;
    catch
        warning('Problem using the function combineSmallTreetoBigTree.Trying again.');        
        % % 重新再生成一遍small bioTree
        addpath(dirImage);
        addpath(dirResultSave);
        for i=1:(length(nameList)-2)
            if i==1
                imageSize=gainOneSmallTree(dirImage,i,dirResultSave);
            else
                imageSize=gainOneSmallTree(dirImage,i,dirResultSave,imageSize);
            end
        end
        rmpath(dirImage)
        rmpath(dirResultSave)

        % set repeat times
        attempTimes=attempTimes+1;
        if attempTimes>3
            msg = ['failedTimes=',num2str(attempTimes),'.' 'Error occurred in combining small tree .'];
            error(msg)
        end
    end
end
disp('Lucky, succeed in combining small tree')

% disp('Now we begin bioTreeMeasure')
% bioTree=bioTreeMeasure(bioTree,0,bioTree{1}.imageSize(1),bioTree{1}.imageSize(2));
% mkdir(dirResultSave,'\allTreeAfterProcesing');
mkdir(dirResultSave,'allTree');
allTreeSave(bioTree,strcat(dirResultSave,'\allTree'),'bioTreeStack');
end
%%
function imageSize=gainOneSmallTree(dirImage,i,dirResultSave,imageSize)
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
% tic;save(saveFile1,'bioTree');toc;
end
function maskImages=myLoad(iJob,fileName,dirResultSave)
if iJob==1
    tempImage=load(fileName);
    maskImages=tempImage.maskImages;
else
    tempImage=load(fileName);
    maskImages=tempImage.maskImages;
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
% function savePoint=findSavePoint(bioTree,treeNum,savePoint)
% if treeNum==1
%     return;
% end
% if ~isempty(savePoint)
%     startFrame=2;
% else
%     startFrame=1;
% end
% if size(savePoint,2)==treeNum-1
%     savePoint=[savePoint,size(bioTree,2)];
%     return;
% end
% for iframe=100:100:size(bioTree,2)
%     bioTreeTmep=bioTree(startFrame:iframe);
%     bytesTemp=whos('bioTreeTmep');
%     bytes=bytesTemp.bytes;
%     pointNum=fix(bytes/2000000000);
%     if pointNum==1
%         savePoint=[savePoint,iframe-100];
%         startFrame=iframe;
%         if size(savePoint,2)==treeNum-1
%             savePoint=[savePoint,size(bioTree,2)];
%             return;
%         end
%     end
% end
% end
% function allTreeSave(bioTree,saveFile)
% saveFile=strcat(saveFile,'\bioTreeStack');
% bytesTemp=whos('bioTree');
% bioTree1=bioTree(1);
% bytesTemp1=whos('bioTree1');
% if bytesTemp1.bytes>2000000000
%     bytes=bytesTemp.bytes-bytesTemp1.bytes;
%     treeNum=fix(bytes/2000000000)+2;
%     savePoint=1;
% else
%     bytes=bytesTemp.bytes;
%     treeNum=fix(bytes/2000000000)+1;
%     savePoint=[];
% end
% savePoint=findSavePoint(bioTree,treeNum,savePoint);
% startFrame=1;
% for iTree=1:treeNum
%     if treeNum==1
%         save(saveFile,'bioTree');
%         return;
%     end
%     bioTreeStack=bioTree(startFrame:savePoint(iTree));
%     saveFile1=strcat(saveFile,num2str(iTree));
%     if savePoint(iTree)==1
%         save(saveFile1,'bioTreeStack','-v7.3');
%     else
%         save(saveFile1,'bioTreeStack');
%     end
%     startFrame=savePoint(iTree)+1;
%     bioTreeStack=[];
% end
% end
function bioTree=combineSmallTreetoBigTree(smallTreeNumber,dirBioTree) % this function cnc be use to connect the bioTree in diffrent time point
bioTree=bioTreeLoader(smallTreeNumber,dirBioTree);
bioTree=leafRootRefineBioTree(bioTree,0);
bioTree=bioTreeSizeReduction(bioTree,0);
% bioTree=nodeEdgeReduction(bioTree,0);
% allTreeSave(bioTree);
end
%%
function bioTree=bioTreeSizeReduction(bioTree,frameShift)
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            %             if bioTree{iframe}.root{iroot}.is2Node==true
            if bioTree{iframe}.root{iroot}.is2Node==true&&~isempty(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList)
                bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList(end)=[];
                %                 bioTree{iframe}.root{iroot}.traceInfo.measurment(end)=[];
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(end)=[];
                    %                     bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment(end)=[];
                end
            end
        end
    end
end
end

%%
function bioTreeCombine=bioTreeLoader(smallTreeNumber,dirBioTree)
% dirBioTree=uigetdir();
cd(dirBioTree);
bioTreeCombine=[];
for i=1:smallTreeNumber
    disp(strcat('Load small Tree',num2str(i)));
    folderName=strcat('t',num2str(i));
    cd(folderName);
    nameList=dir();
    if size(nameList,1)==3
        fileName=strcat('bioTree',num2str(i));
        treeTemp=load(fileName);
        bioTree=treeTemp.bioTree;
        if i==1
            bioTreeCombine=bioTree;
            xSize=bioTree{1}.imageSize(1);
            ySize=bioTree{1}.imageSize(2);
        else
            connectPoint=size(bioTreeCombine,2);
            bioTreeCombine=bioTreeConnecterNew(bioTreeCombine,bioTree,connectPoint,xSize,ySize);
        end
        cd(dirBioTree);
    end
    if size(nameList,1)>=4
        treeIndex=size(nameList,1)-2;
        fileName1=strcat('bioTree',num2str(i),'-',num2str(treeIndex));
        treeTemp1=load(fileName1);
        bioTree=treeTemp1.bioTree;
        preTreeSize=size(treeTemp1.bioTreeStack,2);
        for j=2:size(nameList,1)-2
            treeIndex=size(nameList,1)-j-1;
            fileName1=strcat('bioTree',num2str(i),'-',num2str(treeIndex));
            treeTemp2=load(fileName1);
            bioTree(1:preTreeSize)=treeTemp2.bioTreeStack;
            preTreeSize=size(bioTree,2);
        end
        if i==1
            bioTreeCombine=bioTree;
            xSize=bioTree{1}.imageSize(1);
            ySize=bioTree{1}.imageSize(2);
        else
            connectPoint=size(bioTreeCombine,2);
            bioTreeCombine=bioTreeConnecterNew(bioTreeCombine,bioTree,connectPoint,xSize,ySize);
        end
        cd(dirBioTree);
    end
end
end
function bioTree1=bioTreeConnecterNew(bioTree1,bioTree2,connectPoint,xSize,ySize)
%link the leavies in i frame @bioTree1 and roots in the i+1 frame @bioTree2. (one to one link)
i=1;
%need input correct the iamge size in this function
while i<=size(bioTree1{connectPoint}.leavies,2)
    ileaf=bioTree1{connectPoint}.leavies{i}.leaviesPixelDetail;
    for j=1:size(bioTree2{connectPoint}.root,2)
        jroot=bioTree2{connectPoint}.root{j}.rootPixelDetail;
        if isequal(ileaf,jroot)
            [bioTree1,bioTree2]=linearLinker(bioTree1,bioTree2,connectPoint,i,j);
            i=i-1;
            break;
        end
    end
    i=i+1;
end

if size(bioTree1{connectPoint}.leavies,2)>0
    nodeLeafInTree1=bioTree1{connectPoint}.leavies;
    nodeRootInTree2=bioTree2{connectPoint}.root;
    kNode=createNodeMask(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,xSize,ySize);
    
    for iNode=1:size(kNode,2)
        [nodeIn,nodeOut]=searchingNode(nodeLeafInTree1,nodeRootInTree2,kNode,iNode);
        [bioTree1,bioTree2]=nodeLinker(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,nodeIn,nodeOut,iNode);
    end
end

bioTree1=connectTwoTree(bioTree1,bioTree2,connectPoint);
end
function [bioTree1,bioTree2] =linearLinker(bioTree1,bioTree2,connectPoint,i,j)

if  (bioTree1{connectPoint}.leavies{i}.is2Node==false)&&(bioTree2{connectPoint}.root{j}.is2Node==false)
    rootInTree1=bioTree1{connectPoint}.leavies{i}.rootInfo;
    leafInTree2=bioTree2{connectPoint}.root{j}.leafInfo;
    sizeTrace1=size(bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    %     bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=leafInTree2;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=rootInTree1;
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

if  (bioTree1{connectPoint}.leavies{i}.is2Node==false)&&(bioTree2{connectPoint}.root{j}.is2Node==true)
    rootInTree1=bioTree1{connectPoint}.leavies{i}.rootInfo;
    nodeInTree2=bioTree2{connectPoint}.root{j}.nodeInfo;
    sizeTrace1=size(bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    %     bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.is2Node=true;
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=[];
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.nodeInfo=nodeInTree2;
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=rootInTree1;
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

if  (bioTree1{connectPoint}.leavies{i}.is2Node==true)&&(bioTree2{connectPoint}.root{j}.is2Node==false)
    nodeInTree1=bioTree1{connectPoint}.leavies{i}.nodeInfo;
    leafInTree2=bioTree2{connectPoint}.root{j}.leafInfo;
    sizeTrace1=size(bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    %     bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=leafInTree2;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.is2Node=true;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.nodeInfo=nodeInTree1;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=[];
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

if (bioTree1{connectPoint}.leavies{i}.is2Node==true)&&(bioTree2{connectPoint}.root{j}.is2Node==true)
    nodeInTree1=bioTree1{connectPoint}.leavies{i}.nodeInfo;
    nodeInTree2=bioTree2{connectPoint}.root{j}.nodeInfo;
    sizeTrace1=size(bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    %     bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.is2Node=true;
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=[];
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.nodeInfo=nodeInTree2;
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.isNode=true;
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=[];
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.nodeInfo=nodeInTree1;
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

end
function kNode=createNodeMask(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,xSize,ySize)
allSize=xSize*ySize;
nS1=false(xSize,ySize);
nS2=false(xSize,ySize);
nS3=false(xSize,ySize);

for ileaf=1:size(nodeLeafInTree1,2)
    pixelIdxList2=bioTree1{connectPoint}.leavies{ileaf}.leaviesPixelDetail;
    nS2(pixelIdxList2)=1;
    clear pixelIdxList2;
    
    if nodeLeafInTree1{ileaf}.is2Node==false
        rootInTree1=nodeLeafInTree1{ileaf}.rootInfo;
        if size(bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2)>1
            pixelIdxList1=bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList{end-1};
            nS1(pixelIdxList1)=1;
            clear pixelIdxList1;
        end
        continue;
    end
    if nodeLeafInTree1{ileaf}.is2Node==true
        nodeInTree1=nodeLeafInTree1{ileaf}.nodeInfo;
        if size(bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2)>1
            pixelIdxList1=bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList{end-1};
            nS1(pixelIdxList1)=1;
            clear pixelIdxList1;
        end
        continue;
    end
end

for iroot=1:size(nodeRootInTree2,2)
    
    if size(bioTree2{connectPoint}.root{iroot}.traceInfo.pixelIdxList,2)>1
        pixelIdxList3=bioTree2{connectPoint}.root{iroot}.traceInfo.pixelIdxList{2};
        nS3(pixelIdxList3)=1;
        clear pixelIdxList3;
    end
    
end

nodeStacks(:,:,1)=nS1;
nodeStacks(:,:,2)=nS2;
nodeStacks(:,:,3)=nS3;
CCnode=bwconncomp(nodeStacks,26);

for iObject=1:CCnode.NumObjects
    pixelFrame=fix(CCnode.PixelIdxList{iObject}./allSize);
    kNode{iObject}=CCnode.PixelIdxList{iObject}(pixelFrame==1)-allSize;
end
end
function [nodeIn,nodeOut]=searchingNode(nodeLeafInTree1,nodeRootInTree2,kNode,iNode)
nodeIn=[];
nodeOut=[];
for iNodeIn=1:size(nodeLeafInTree1,2)
    if ismember(nodeLeafInTree1{iNodeIn}.leaviesPixelDetail,kNode{iNode})
        nodeIn=[nodeIn,iNodeIn];
    end
end

for iNodeOut=1:size(nodeRootInTree2,2)
    if ismember(nodeRootInTree2{iNodeOut}.rootPixelDetail,kNode{iNode})
        nodeOut=[nodeOut,iNodeOut];
    end
end

end
function [bioTree1,bioTree2]=nodeLinker(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,nodeIn,nodeOut,iNode)

for iNodeIn=1:size(nodeIn,2)
    if  nodeLeafInTree1{nodeIn(iNodeIn)}.is2Node==false
        rootInTree1=nodeLeafInTree1{nodeIn(iNodeIn)}.rootInfo;
        bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.is2Node=true;
        bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.nodeInfo=[connectPoint,iNode,iNodeIn];
        bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=[];
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.isNode=false;
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.rootInfo=rootInTree1;
        clear rootInTree1;
        continue;
    end
    
    if  nodeLeafInTree1{nodeIn(iNodeIn)}.is2Node==true
        nodeInTree1=nodeLeafInTree1{nodeIn(iNodeIn)}.nodeInfo;
        bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.is2Node=true;
        bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=[];
        bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.nodeInfo=[connectPoint,iNode,iNodeIn];
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.isNode=true;
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.nodeInfo=nodeInTree1;
        clear nodeInTree1;
        continue;
    end
end

for iNodeOut=1:size(nodeOut,2)
    if nodeRootInTree2{nodeOut(iNodeOut)}.is2Node==false
        leafInTree2=nodeRootInTree2{nodeOut(iNodeOut)}.leafInfo;
        bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.is2Node=true;
        bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=[];
        bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.nodeInfo=[connectPoint,iNode,iNodeOut];
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.is2Node=false;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.leafInfo=leafInTree2;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.traceInfo=nodeRootInTree2{nodeOut(iNodeOut)}.traceInfo;
        clear leafInTree2;
        continue;
    end
    
    if nodeRootInTree2{nodeOut(iNodeOut)}.is2Node==true
        nodeInTree2=nodeRootInTree2{nodeOut(iNodeOut)}.nodeInfo;
        bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.isNode=true;
        bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=[];
        bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.nodeInfo=[connectPoint,iNode,iNodeOut];
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.is2Node=true;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.nodeInfo=nodeInTree2;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.traceInfo=nodeRootInTree2{nodeOut(iNodeOut)}.traceInfo;
        clear nodeInTree2;
        continue;
    end
    
end

end
function bioTree1=connectTwoTree(bioTree1,bioTree2,connectPoint)
bioTree1{connectPoint}.leavies=[];
bioTree2{connectPoint}.root=[];
bioTree1{connectPoint}.leavies=bioTree2{connectPoint}.leavies;
bioTree1{connectPoint}.node=bioTree2{connectPoint}.node;
bioTree1(connectPoint+1:1:size(bioTree2,2))=bioTree2(connectPoint+1:1:size(bioTree2,2));
end
