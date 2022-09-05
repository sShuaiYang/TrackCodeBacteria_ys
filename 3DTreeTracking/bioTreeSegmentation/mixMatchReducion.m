%% main function 
% find node, do the match, create new node, correct the relationShip
function bioTree=mixMatchReducion(bioTree,threShold)
% in<out
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=mixMatchNodeinFrame(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                matchResult=matchProcess(bioTree,nodeList(iList,:),threShold,'in>out');
                if ~isempty(matchResult)
                    bioTree=nodeCreation(bioTree,nodeList(iList,:),matchResult);
                    bioTree=nodeCorrection(bioTree,nodeList(iList,:));
                    nodeList(iList+1:end,2)=nodeList(iList+1:end,2)-1;
                end
            end
        end
    end
end
fprintf('\n')
bioTree=unbalanceNodeReduction(bioTree,threShold);
% in>out
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=mixMatchNodeinFrame2(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                matchResult=matchProcess(bioTree,nodeList(iList,:),threShold,'in<out');
                if ~isempty(matchResult)
                    bioTree=nodeCreation2(bioTree,nodeList(iList,:),matchResult);
                    bioTree=nodeCorrection(bioTree,nodeList(iList,:));
                    nodeList(iList+1:end,2)=nodeList(iList+1:end,2)-1;
                end
            end
        end
    end
end
end
%% function 1
function nodeList=mixMatchNodeinFrame(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if size(bioTreeFrame.node{iNode}.In,2) >= size(bioTreeFrame.node{iNode}.Out,2) && size(bioTreeFrame.node{iNode}.Out,2)>=3 ||...
                size(bioTreeFrame.node{iNode}.In,2)>=4 && size(bioTreeFrame.node{iNode}.Out,2)==2
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function matchResult=matchProcess(bioTree,nodeInfo,threShold,type)
pictureSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
matchResult=individualToMixMatch(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,pictureSize,threShold,type);
end
function bioTree=nodeCreation(bioTree,nodeInfo,matchResult)
pre_Node_Info=bioTree{nodeInfo(1)}.node{nodeInfo(2)};
bioTree{nodeInfo(1)}.node(nodeInfo(2))=[];
nodeNumNew=size(bioTree{nodeInfo(1)}.node,2);
realNode=0;
for iOut=1:size(pre_Node_Info.Out,2)
    inNum=numel(matchResult{iOut,1});
    if inNum>=2
        realNode=realNode+1;
        newNum=nodeNumNew+realNode;
        for iIn=1:inNum
            bioTree{nodeInfo(1)}.node{newNum}.In(iIn)=pre_Node_Info.In(matchResult{iOut,1}(iIn));
        end
        bioTree{nodeInfo(1)}.node{newNum}.Out(1)=pre_Node_Info.Out(iOut);
    end
    if inNum==1
        preIsNode=pre_Node_Info.In{matchResult{iOut,1}(1)}.isNode;
        if preIsNode==1
            preNodeInfo=pre_Node_Info.In{matchResult{iOut,1}(1)}.nodeInfo;
            bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.is2Node=pre_Node_Info.Out{iOut}.is2Node;
            if bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.is2Node==1
                nextNode=pre_Node_Info.Out{iOut}.nodeInfo;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.leafInfo=[];
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.nodeInfo=nextNode;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=preNodeInfo;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
            else
                nextLeaf=pre_Node_Info.Out{iOut}.leafInfo;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.leafInfo=nextLeaf;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=preNodeInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
            end
            bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList,pre_Node_Info.Out{iOut}.traceInfo.pixelIdxList];
        end
        if preIsNode==0
            preRootInfo=pre_Node_Info.In{matchResult{iOut,1}(1)}.rootInfo;
            bioTree{preRootInfo(1)}.root{preRootInfo(2)}.is2Node=pre_Node_Info.Out{iOut}.is2Node;
            if  bioTree{preRootInfo(1)}.root{preRootInfo(2)}.is2Node==1
                nextNode=pre_Node_Info.Out{iOut}.nodeInfo;
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.leafInfo=[];
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.nodeInfo=nextNode;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=preRootInfo;
            else
                nextLeaf=pre_Node_Info.Out{iOut}.leafInfo;
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.leafInfo=nextLeaf;
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=preRootInfo;
            end
            bioTree{preRootInfo(1)}.root{preRootInfo(2)}.traceInfo.pixelIdxList=[bioTree{preRootInfo(1)}.root{preRootInfo(2)}.traceInfo.pixelIdxList,pre_Node_Info.Out{iOut}.traceInfo.pixelIdxList];
        end
    end
end
end
function bioTree=nodeCorrection(bioTree,nodeInfo)
iframe=nodeInfo(1);
for iNode=1:size(bioTree{iframe}.node,2)
    for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
        if bioTree{iframe}.node{iNode}.In{iIn}.isNode==1
            preNode=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[iframe,iNode,iIn];
        end
        if bioTree{iframe}.node{iNode}.In{iIn}.isNode==0
            preRoot=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[iframe,iNode,iIn];
        end
    end
    for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
        if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==1
            nextNode=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[iframe,iNode,iOut];
        end
        if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==0
            nextLeaf=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[iframe,iNode,iOut];
        end
    end
end
end
%% match process
function matchResult=individualToMixMatch(pixelIdxListIn,bioTreeOut,pictureSize,diff,type)
if strcmp(type,'in>out')
    regionInfo=getSumRegion(pixelIdxListIn,bioTreeOut,pictureSize);
    for iIn=1:numel(pixelIdxListIn)
        [xyMin,regionImage,isbreak]=basicTestForIndividual(pixelIdxListIn{iIn},pictureSize,regionInfo);
        if isbreak==1
            matchResult=[];
            return
        else
            in{iIn}.xyMin=xyMin;
            in{iIn}.BWImage=regionImage;
        end
    end
    for iOut=1:numel(bioTreeOut)
        [regionNum,xyMin,regionImage]=findRegionNum(bioTreeOut{iOut}.traceInfo.pixelIdxList{1},pictureSize,regionInfo);
        if regionNum>=2
            matchResult=[];
            return
        else
            out{iOut}.xyMin=xyMin;
            out{iOut}.BWImage=regionImage{1};
        end
    end
end
if strcmp(type,'in<out')
    regionInfo=getSumRegion(pixelIdxListIn,bioTreeOut,pictureSize);
    for iIn=1:numel(pixelIdxListIn)
        [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListIn{iIn},pictureSize,regionInfo);
        if regionNum>=2
            matchResult=[];
            return
        else
            in{iIn}.xyMin=xyMin;
            in{iIn}.BWImage=regionImage{1};
        end
    end
    for iOut=1:numel(bioTreeOut)
        [xyMin,regionImage,isbreak]=basicTestForIndividual(bioTreeOut{iOut}.traceInfo.pixelIdxList{1},pictureSize,regionInfo);
        if isbreak==1
            matchResult=[];
            return
        else
            out{iOut}.xyMin=xyMin;
            out{iOut}.BWImage=regionImage;
        end
    end
    trans=in;
    in=out;
    out=trans;
end
matchResult=matchMix(in,out,diff);
end
function regionInfo=getSumRegion(pixelIdxListIn,bioTreeOut,pictureSize)
% get the min region for input and output images
AllList=[];
for i=1:numel(pixelIdxListIn)
    AllList=[AllList;pixelIdxListIn{i}];
end
for i=1:numel(bioTreeOut)
    AllList=[AllList;bioTreeOut{i}.traceInfo.pixelIdxList{1}];
end
xSize=pictureSize(1);
yresult=ceil(AllList/xSize);
xresult=AllList-(yresult-1)*xSize;
regionInfo.xMin=min(xresult);
regionInfo.xMax=max(xresult);
regionInfo.yMin=min(yresult);
regionInfo.yMax=max(yresult);
end
function [xyMin,regionImage,isbreak]=basicTestForIndividual(iIn,pictureSize,regionInfo)
% detect regionNum
[regionNum,xyMin,regionImage]=findRegionNum(iIn,pictureSize,regionInfo);
if regionNum>=2
    isbreak=1;
    return
end
regionImage=regionImage{1};
% detect whether the region is convex polygon
% areaDiffThr=0.3;
stats=regionprops(regionImage,'FilledArea','MajorAxisLength','MinorAxisLength');
% if abs(stats.MajorAxisLength*stats.MinorAxisLength*pi/4-stats.FilledArea)/stats.FilledArea>=areaDiffThr
%     isbreak=1;
%     return
% end     %% 觉得不需要，本质在于配对并验证，不需要管是不是几个元素在一起，单连通还是多连通

% detect whether the region is too big
if stats.FilledArea>=800
    isbreak=1;
    return
end
isbreak=0;
end

%% easy match and test the result
function matchResult=matchMix(in,out,diff)
selectionNum=1:numel(in);
for iOut=1:numel(out)
    matchResult{iOut,1}=[];
    deleteI=[];
    for iSel=1:numel(selectionNum)
        target=in{selectionNum(iSel)};
        overlapImage=target.BWImage & out{iOut}.BWImage;
        if double(numel(overlapImage(overlapImage==1)))/numel(target.BWImage(target.BWImage==1))>=0.5
            matchResult{iOut,1}=[matchResult{iOut,1},selectionNum(iSel)];
            deleteI=[deleteI,iSel];
        end
    end
    if isempty(matchResult{iOut,1})
        matchResult=[];
        return
    end
    selectionNum(deleteI)=[];
end
matchNum=0;
for i=1:numel(out)
    matchNum=matchNum+numel(matchResult{i,1});
end
if matchNum~=numel(in)
    matchResult=[];
    return
end
matchResult=matchResultCheckOut(matchResult,in,out,diff);
end
function matchResult=matchResultCheckOut(matchResult,in,out,diff)
for i=1:numel(out)
    inInfo=[];
    outInfo=[];
    if numel(matchResult{i,1})==1
        statsIn=regionprops(in{matchResult{i,1}}.BWImage,'MajorAxisLength','MinorAxisLength','Centroid');
        inInfo(1,:)=[statsIn.MajorAxisLength,statsIn.MinorAxisLength,statsIn.Centroid(1),statsIn.Centroid(2)];
        statsOut=regionprops(out{i}.BWImage,'MajorAxisLength','MinorAxisLength','Centroid');
        outInfo(1,:)=[statsOut.MajorAxisLength,statsOut.MinorAxisLength,statsOut.Centroid(1),statsOut.Centroid(2)];
        dist=pdist2(inInfo,outInfo);
        if dist>=diff
            matchResult=[];
            return
        end
    end
    if numel(matchResult{i,1})>=2
        mergeImage=false(size(in{1}.BWImage));
        for iSIn=1:numel(matchResult{i,1})
            mergeImage=mergeImage | in{matchResult{i,1}(iSIn)}.BWImage;
        end
        inArea=numel(mergeImage(mergeImage==1));
        outArea=numel(out{i}.BWImage(out{i}.BWImage==1));
        if iSIn<=2
            areaThre=150;
        else
            areaThre=250;
        end
        if abs(inArea-outArea)>=areaThre
            matchResult=[];
            return
        end
    end
end
end
%% basic functions
function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
function [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize,regionInfo)
% find how many regions there are in a input/output
if isempty(pixelIdxList)
    regionNum=0;
    xyMin=zeros(1,2);
    regionImage{1}=[];
else
    [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize,regionInfo);
    CC=bwconncomp(BWImage);
    regionNum=CC.NumObjects;
    for i=1:regionNum
        pixelIdxList2=CC.PixelIdxList{i};
        BW=false(CC.ImageSize);
        BW(pixelIdxList2)=1;
        regionImage{i}=BW;
    end
end
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize,regionInfo)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=regionInfo.xMin;
xMax=regionInfo.xMax;
yMin=regionInfo.yMin;
yMax=regionInfo.yMax;
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
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
%% unbalance node reduction(function 2)
% an unusual situation, two bacteria that is nearly to divide, 
% 1,[1,1] -- [1,1] ,1, which could not be match by other code
function bioTree=unbalanceNodeReduction(bioTree,threShold)
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        nodeList=unbalanceNodeinFrame(bioTree{iframe},iframe);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                matchResult=unbalanceMatchProcess(bioTree,nodeList(iList,:),threShold);
                if ~isempty(matchResult)
                    bioTree=nodeCreation(bioTree,nodeList(iList,:),matchResult);
                    bioTree=nodeCorrection(bioTree,nodeList(iList,:));
                    nodeList(iList+1:end,2)=nodeList(iList+1:end,2)-1;
                end
            end
        end
    end
end
end
function nodeList=unbalanceNodeinFrame(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==3
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end    
function matchResult=unbalanceMatchProcess(bioTree,nodeInfo,threShold)
pictureSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
matchResult=unbalanceMixMatch(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,pictureSize,threShold);
end
function matchResult=unbalanceMixMatch(pixelIdxListIn,bioTreeOut,pictureSize,diff)
regionInfo=getSumRegion(pixelIdxListIn,bioTreeOut,pictureSize);
for iIn=1:numel(pixelIdxListIn)
    [xyMin,regionImage,isbreak]=basicTestForIndividual(pixelIdxListIn{iIn},pictureSize,regionInfo);
    if isbreak==1
        matchResult=[];
        return
    else
        in{iIn}.xyMin=xyMin;
        in{iIn}.BWImage=regionImage;
    end
end
for iOut=1:numel(bioTreeOut)
    [regionNum,xyMin,regionImage]=findRegionNum(bioTreeOut{iOut}.traceInfo.pixelIdxList{1},pictureSize,regionInfo);
    if regionNum>=2
        matchResult=[];
        return
    else
        out{iOut}.xyMin=xyMin;
        out{iOut}.BWImage=regionImage{1};
    end
end
matchResult=unbalanceMatchMix(in,out,diff);
end  
function matchResult=unbalanceMatchMix(in,out,diff)
area=[];
for iIn=1:numel(in);
    area=[area;numel(in{iIn}.BWImage(in{iIn}.BWImage==1))];
end
for iOut=1:numel(out);
    area=[area;numel(in{iOut}.BWImage(out{iOut}.BWImage==1))];
end
number=1:3;
matchResult=[];
if numel(area>=(max(area)-20))<=2
    maxArea=find(area==max(area));
    maxArea=maxArea(1);
    if maxArea<=3
        matchResult{1,1}=maxArea;
        matchResult{1,2}=[];
        matchResult{2,1}=number(number~=maxArea);
        matchResult{2,2}=[];
        for iOut=1:3
            bigIn=in{maxArea}.BWImage;
            overlapImage=bigIn & out{iOut}.BWImage;
            if double(numel(overlapImage(overlapImage==1)))/numel(out{iOut}.BWImage(out{iOut}.BWImage==1))>=0.6
                matchResult{1,2}=[matchResult{1,2},iOut];
            end
        end
        if numel(matchResult{1,2})==2
            matchResult{2,2}=number(~ismember(number,matchResult{1,2}));
        else
            matchResult=[];
            return
        end
    end
    if maxArea>=4
        matchResult{1,2}=maxArea-3;
        matchResult{2,2}=number(number~=(maxArea-3));
        matchResult{1,1}=[];
        matchResult{2,1}=[];
        for iIn=1:3
            bigOut=out{maxArea-3}.BWImage;
            overlapImage=bigOut & in{iIn}.BWImage;
            if double(numel(overlapImage(overlapImage==1)))/numel(in{iIn}.BWImage(in{iIn}.BWImage==1))>=0.6
                matchResult{1,1}=[matchResult{1,1},iIn];
            end
        end
        if numel(matchResult{1,1})==2
            matchResult{1,2}=number(~ismember(number,matchResult{1,1}));
        else
            matchResult=[];
            return
        end
    end
    matchResult=testUnbalanceMatchResult(matchResult,in,out,diff);
end
end
function matchResult=testUnbalanceMatchResult(matchResult,in,out,diff)
if numel(matchResult{1,1})==1
    matchResult=matchResult(end:-1:1,:);
end
% 2-1
mergeImage=false(size(in{1}.BWImage));
for iIn=1:2
    mergeImage=mergeImage | in{matchResult{1,1}(iIn)}.BWImage;
end
overlapImage=mergeImage & out{matchResult{1,2}(iOut)}.BWImage;

for i=1:numel(out)
    inInfo=[];
    outInfo=[];
    if numel(matchResult{i,1})==1
        statsIn=regionprops(in{matchResult{i,1}}.BWImage,'MajorAxisLength','MinorAxisLength','Centroid');
        inInfo(1,:)=[statsIn.MajorAxisLength,statsIn.MinorAxisLength,statsIn.Centroid(1),statsIn.Centroid(2)];
        statsOut=regionprops(out{1}.BWImage,'MajorAxisLength','MinorAxisLength','Centroid');
        outInfo(1,:)=[statsIn.MajorAxisLength,statsIn.MinorAxisLength,statsIn.Centroid(1),statsIn.Centroid(2)];
        dist=pdist2(inInfo,outInfo);
        if dist>=diff
            matchResult=[];
            return
        end
    end
    if numel(matchResult{i,1})>=2
        mergeImage=false(size(in{1}.BWImage));
        for iSIn=1:numel(matchResult{i,1})
            mergeImage=mergeImage | in{matchResult{i,1}(iSIn)}.BWImage;
        end
        inArea=numel(mergeImage(mergeImage==1));
        outArea=numel(out{i}.BWImage(out{i}.BWImage==1));
        if iSIn<=2
            areaThre=150;
        else
            areaThre=250;
        end
        if abs(inArea-outArea)>=areaThre
            matchResult=[];
            return
        end
    end
end
end

%% function 3
function nodeList=mixMatchNodeinFrame2(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if size(bioTreeFrame.node{iNode}.Out,2) > size(bioTreeFrame.node{iNode}.In,2) && size(bioTreeFrame.node{iNode}.In,2)>=3 ||...
                size(bioTreeFrame.node{iNode}.Out,2)>=4 && size(bioTreeFrame.node{iNode}.In,2)==2
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function bioTree=nodeCreation2(bioTree,nodeInfo,matchResult)
pre_Node_Info=bioTree{nodeInfo(1)}.node{nodeInfo(2)};
bioTree{nodeInfo(1)}.node(nodeInfo(2))=[];
nodeNumNew=size(bioTree{nodeInfo(1)}.node,2);
realNode=0;
for iIn=1:size(pre_Node_Info.In,2)
    outNum=numel(matchResult{iIn,1});
    if outNum>=2
        realNode=realNode+1;
        newNum=nodeNumNew+realNode;
        for iOut=1:outNum
            bioTree{nodeInfo(1)}.node{newNum}.Out(iOut)=pre_Node_Info.Out(matchResult{iIn,1}(iOut));
        end
        bioTree{nodeInfo(1)}.node{newNum}.In(1)=pre_Node_Info.In(iIn);
    end
    if outNum==1
        preIsNode=pre_Node_Info.In{iIn}.isNode;
        if preIsNode==1
            preNodeInfo=pre_Node_Info.In{iIn}.nodeInfo;
            bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.is2Node=pre_Node_Info.Out{matchResult{iIn,1}(1)}.is2Node;
            if bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.is2Node==1
                nextNode=pre_Node_Info.Out{matchResult{iIn,1}(1)}.nodeInfo;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.leafInfo=[];
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.nodeInfo=nextNode;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=preNodeInfo;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
            else
                nextLeaf=pre_Node_Info.Out{matchResult{iIn,1}(1)}.leafInfo;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.leafInfo=nextLeaf;
                bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=preNodeInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
            end
            bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out{preNodeInfo(3)}.traceInfo.pixelIdxList,pre_Node_Info.Out{matchResult{iIn,1}(1)}.traceInfo.pixelIdxList];
        end
        if preIsNode==0
            preRootInfo=pre_Node_Info.In{iIn}.rootInfo;
            bioTree{preRootInfo(1)}.root{preRootInfo(2)}.is2Node=pre_Node_Info.Out{matchResult{iIn,1}(1)}.is2Node;
            if  bioTree{preRootInfo(1)}.root{preRootInfo(2)}.is2Node==1
                nextNode=pre_Node_Info.Out{matchResult{iIn,1}(1)}.nodeInfo;
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.leafInfo=[];
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.nodeInfo=nextNode;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=preRootInfo;
            else
                nextLeaf=pre_Node_Info.Out{matchResult{iIn,1}(1)}.leafInfo;
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.leafInfo=nextLeaf;
                bioTree{preRootInfo(1)}.root{preRootInfo(2)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=preRootInfo;
            end
            bioTree{preRootInfo(1)}.root{preRootInfo(2)}.traceInfo.pixelIdxList=[bioTree{preRootInfo(1)}.root{preRootInfo(2)}.traceInfo.pixelIdxList,pre_Node_Info.Out{matchResult{iIn}(1)}.traceInfo.pixelIdxList];
        end
    end
end
end