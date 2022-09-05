function bioTree=troubleNodeSolution(bioTree)
dirFile='D:\troubleNodeSaveResult';
mkdir(dirFile);
cd(dirFile)
if nargin==1
    [bioTree,~,~,~]=myBiograph_new2(bioTree);
    % find each trouble node pixelIdxInfo
    troubleNodeList=findTroubleNode(bioTree);
    save('troubleNodeList.mat','troubleNodeList')
else
    load('troubleNodeList')
end
% bioTree{1}.imageSize=[188,203];
% load('solution.mat')
% solution for each trouble node
for iNode=1:numel(troubleNodeList)
    clc
    disp([iNode,troubleNodeList{iNode}.nodeInfo])
    solution{iNode}=manualSegmentation(troubleNodeList{iNode}.pixelIdxListIn,troubleNodeList{iNode}.pixelIdxListOut,troubleNodeList{iNode}.is2Node,bioTree{1}.imageSize);
    solution{iNode}.nodeInfo=troubleNodeList{iNode}.nodeInfo;
    save('solution.mat','solution')  % save for every 10 node so that you work would not be broken
end
save('solution.mat','solution')
% change bioTree with solutionResult
bioTree=solveTroubleNode(bioTree,solution);
end
%% find each trouble node
function troubleNodeList=findTroubleNode(bioTree)
i=0;
for iframe=1:size(bioTree,2)
    nodeList=troubleNodeinFrame(bioTree{iframe},iframe);
    for iNode=1:size(nodeList)
        nodeInfo=nodeList(iNode,:);
        i=i+1;
        troubleNodeList{i}.pixelIdxListIn=getInputMask(bioTree,nodeInfo);
        troubleNodeList{i}.nodeInfo=nodeInfo;
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            troubleNodeList{i}.pixelIdxListOut{iOut}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
        end
        is2Node=[];
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            is2Node=[is2Node;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node];
        end
        troubleNodeList{i}.is2Node=is2Node;
    end
end
end
function nodeList=troubleNodeinFrame(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if ~(size(bioTreeFrame.node{iNode}.In,2)==1 && size(bioTreeFrame.node{iNode}.Out,2)==2)
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
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
%% change bioTree with solutionResult
function bioTree=solveTroubleNode(bioTree,solution)
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        bioTree{iframe}.node{iNode}.state=[];
    end
end
nodeList=[];
for iSol=1:numel(solution)
    nodeList=[nodeList;solution{iSol}.nodeInfo];
end
nodeFrame=unique(nodeList(:,1));
for i=1:numel(nodeFrame)
    iframe=nodeFrame(i);
    if iframe>=78
        p=1;
    end
    disp(iframe)
    solutionNum=find(nodeList(:,1)==iframe);
    for iNode=1:numel(solutionNum)
        switch solution{solutionNum(iNode)}.state
            case 'pass'
                continue
            case 'sep'
                [bioTree,changedNode]=bioTreeNodeReplace(bioTree,solution{solutionNum(iNode)});
                for iLine=1:size(changedNode,1)
                    nodeInfo=changedNode(iLine,:);
                    focusNode=find((nodeInfo(1)==nodeList(:,1)&nodeInfo(2)==nodeList(:,2))==1);
                    if ~isempty(focusNode)
                        solution{focusNode}.state='pass';
                    end
                end
            case 'deal'
                bioTree=bioTreeSolveTroubleNode(bioTree,solution{solutionNum(iNode)});
            case 'mustDiv'
                [bioTree,nodeInfo]=bioTreeMustDiv(bioTree,solution{solutionNum(iNode)});
                if ~isempty(nodeInfo)
                    focusNode=find((nodeInfo(1)==nodeList(:,1)&nodeInfo(2)==nodeList(:,2))==1);
                    if ~isempty(focusNode)
                        solution{focusNode}.state='pass';
                    end
                end
            case 'mistake'
                nodeInfo=solution{solutionNum(iNode)}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.state='mistake';
        end
    end
    bioTree=bioTreeLinkCorrection(bioTree,iframe);
end
bioTree=correctForMistake(bioTree);
end
function [bioTree,changedNode]=bioTreeNodeReplace(bioTree,solutionI)
nodeInfo=solutionI.nodeInfo;
realOutNumber=unique(solutionI.outNum);
aimNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)};
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out=[];
thisOut=0;
changedNode=[];
for iRealOut=realOutNumber
    if numel(aimNode.Out{iRealOut}.traceInfo.pixelIdxList)==1   %next frame is a node/a leaf.  if is a node, change the nextNode inInfo, if is a leaf, generate enough leavies.
        is2Node=aimNode.Out{iRealOut}.is2Node;
        if is2Node==0
            nextLeaf=aimNode.Out{iRealOut}.leafInfo;
            matchOut=find(solutionI.outNum==iRealOut);
            for iOut=1:numel(matchOut)
                if iOut==1
                    leafInfo=nextLeaf;
                else
                    leafInfo=[nextLeaf(1),size(bioTree{nextLeaf(1)}.leavies,2)+1];
                end
                thisOut=thisOut+1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.traceInfo.pixelIdxList{1}=solutionI.pixelIdxList{matchOut(iOut)};
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.is2Node=0;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.leafInfo=leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=1;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[nodeInfo,thisOut];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=solutionI.pixelIdxList{matchOut(iOut)};
            end
        end
        if is2Node==1
            nextNode=aimNode.Out{iRealOut}.nodeInfo;
            matchOut=find(solutionI.outNum==iRealOut);
            if numel(matchOut)~=1
                changedNode=[changedNode;nextNode];
            end
            for iOut=1:numel(matchOut)
                realNode=nextNode;
                thisOut=thisOut+1;
                if iOut==1
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.traceInfo.pixelIdxList{1}=solutionI.pixelIdxList{matchOut(iOut)};
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.is2Node=1;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.nodeInfo=realNode;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.leafInfo=[];
                    bioTree{realNode(1)}.node{realNode(2)}.In{realNode(3)}.isNode=1;
                    bioTree{realNode(1)}.node{realNode(2)}.In{realNode(3)}.nodeInfo=[nodeInfo,thisOut];
                    bioTree{realNode(1)}.node{realNode(2)}.In{realNode(3)}.rootInfo=[];
                else
                    nextInNum=size(bioTree{realNode(1)}.node{realNode(2)}.In,2)+1;
                    realNode(3)=nextInNum;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.traceInfo.pixelIdxList{1}=solutionI.pixelIdxList{matchOut(iOut)};
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.is2Node=1;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.nodeInfo=realNode;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.leafInfo=[];
                    bioTree{realNode(1)}.node{realNode(2)}.In{realNode(3)}.isNode=1;
                    bioTree{realNode(1)}.node{realNode(2)}.In{realNode(3)}.nodeInfo=[nodeInfo,thisOut];
                    bioTree{realNode(1)}.node{realNode(2)}.In{realNode(3)}.rootInfo=[];
                end  
            end
        end
    end
    if numel(aimNode.Out{iRealOut}.traceInfo.pixelIdxList)>1
        matchOut=find(solutionI.outNum==iRealOut);
        if numel(matchOut)==1
            thisOut=thisOut+1;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}=aimNode.Out{iRealOut};
            is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.is2Node;
            if is2Node==1
                nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.nodeInfo;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[nodeInfo,thisOut];
            end
            if is2Node==0
                nextLeaf=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.leafInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[nodeInfo,thisOut];
            end
        end
        if numel(matchOut)>1
            nextNodeInfo=[nodeInfo(1)+1,size(bioTree{nodeInfo(1)+1}.node,2)+1];
            bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.state=[];
            for iOut=1:numel(matchOut)
                thisOut=thisOut+1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.traceInfo.pixelIdxList{1}=solutionI.pixelIdxList{matchOut(iOut)};
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.is2Node=1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.nodeInfo=[nextNodeInfo,iOut];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{thisOut}.leafInfo=[];
                bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{iOut}.isNode=1;
                bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{iOut}.nodeInfo=[nodeInfo,thisOut];
                bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{iOut}.rootInfo=[];
            end
            bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{1}=aimNode.Out{iRealOut};
            bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{1}.traceInfo.pixelIdxList(1)=[];
            is2Node=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{1}.is2Node;
            if is2Node==1
                nextNode=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{1}.nodeInfo;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[nextNodeInfo,1];
            end
            if is2Node==0
                nextLeaf=bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.Out{1}.leafInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[nextNodeInfo,1];
            end
        end
    end
end
end
function bioTree=bioTreeSolveTroubleNode(bioTree,solutionI)
nodeInfo=solutionI.nodeInfo;
inSet=[];
outSet=[];
for i=1:numel(solutionI.matchResult)
    index=find(solutionI.matchResult{i}==0);
    inSet=[inSet,solutionI.matchResult{i}(1:index-1)];
    outSet=[outSet,solutionI.matchResult{i}(index+1:end)];
end

%  case for in-empty(generate a leaf for in)
for iIn=1:numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In) %generate leaf
    if ~ismember(iIn,inSet)
        isNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode;
        if isNode==1
            preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=0;    
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[];
            leafInfo=[nodeInfo(1)-1,size(bioTree{nodeInfo(1)-1}.leavies,2)+1];
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=1;    
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=preNode;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList{end};
        end
        if isNode==0
            preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=0;
            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[];
            leafInfo=[nodeInfo(1)-1,size(bioTree{nodeInfo(1)-1}.leavies,2)+1];
            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=0;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=preRoot;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail= bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList{end};
        end
    end
end

%  case for out-empty(generate a root for out)
outNum=numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out);
for iOut=1:outNum
    if ~ismember(iOut,outSet)
        newRootInfo=[nodeInfo(1),size(bioTree{nodeInfo(1)}.root,2)+1];
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==1
            nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.is2Node=1;
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.leafInfo=[];
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.nodeInfo=nextNode;
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.rootPixelDetail=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=newRootInfo;
        end
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==0
            nextLeaf=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.is2Node=0;
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.leafInfo=nextLeaf;
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.nodeInfo=[];
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.rootPixelDetail=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
            bioTree{newRootInfo(1)}.root{newRootInfo(2)}.traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=newRootInfo;
        end
    end
end

% right match for each in and outPair
aimNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)};
bioTree{nodeInfo(1)}.node{nodeInfo(2)}=[];
realNode=0;
for iMatch=1:numel(solutionI.matchResult)
    index=find(solutionI.matchResult{iMatch}==0);
    inIndex=solutionI.matchResult{iMatch}(1:index-1);
    outIndex=solutionI.matchResult{iMatch}(index+1:end);
    if numel(inIndex)==1 && numel(outIndex)==1
        preIsNode=aimNode.In{inIndex}.isNode;
        nextIsNode=aimNode.Out{outIndex}.is2Node;
        if preIsNode==0
            preRoot=aimNode.In{inIndex}.rootInfo;
            if nextIsNode==1
                nextNode=aimNode.Out{outIndex}.nodeInfo;
                bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
                bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
                bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=nextNode;
                bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,aimNode.Out{outIndex}.traceInfo.pixelIdxList];
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=preRoot;
            end
            if nextIsNode==0
                nextLeaf=aimNode.Out{outIndex}.leafInfo;
                bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=0;
                bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=nextLeaf;
                bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[];
                bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,aimNode.Out{outIndex}.traceInfo.pixelIdxList];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=preRoot;
            end
        end
        if preIsNode==1
            preNode=aimNode.In{inIndex}.nodeInfo;
            if nextIsNode==1
                nextNode=aimNode.Out{outIndex}.nodeInfo;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=nextNode;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,aimNode.Out{outIndex}.traceInfo.pixelIdxList];
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=preNode;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
            end
            if nextIsNode==0
                nextLeaf=aimNode.Out{outIndex}.leafInfo;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=0;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=nextLeaf;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nextNode=[];
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,aimNode.Out{outIndex}.traceInfo.pixelIdxList];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=preNode;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
            end
        end
    end
    if numel(inIndex)>1 || numel(outIndex)>1
        realNode=realNode+1;
        if realNode==1
            newNode=[nodeInfo(1),nodeInfo(2)];
            bioTree{newNode(1)}.node{newNode(2)}.state=[];
        else
            newNode=[nodeInfo(1),numel(bioTree{nodeInfo(1)}.node)+1];
            bioTree{newNode(1)}.node{newNode(2)}.state=[];
        end
        for iIn=1:numel(inIndex)
            bioTree{newNode(1)}.node{newNode(2)}.In{iIn}=aimNode.In{inIndex(iIn)};
        end
        for iOut=1:numel(outIndex)
            bioTree{newNode(1)}.node{newNode(2)}.Out{iOut}=aimNode.Out{outIndex(iOut)};
        end
    end
end
end
function [bioTree,changeNode]=bioTreeMustDiv(bioTree,solutionI)
nodeInfo=solutionI.nodeInfo;
if numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In)>=2 && numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out)==1
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.is2Node==1
        changeNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1}.nodeInfo;
    else
        changeNode=[];
    end
    [bioTree,usefulNode]=type1NodeReduction(bioTree,3,40,nodeInfo);
    if usefulNode==1
        changeNode=[];
    end
else if numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In) == numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out)
        bioTree=type2NodeReduction(bioTree,2,nodeInfo);
        changeNode=[];
    end
end
end
function bioTree=bioTreeLinkCorrection(bioTree,iframe)
emptyNode=[];
for iNode=1:numel(bioTree{iframe}.node)
    if isempty(bioTree{iframe}.node{iNode})
        emptyNode=[emptyNode;iNode];
    end
end
bioTree{iframe}.node(emptyNode)=[];
for iNode=1:numel(bioTree{iframe}.node)
    for iIn=1:numel(bioTree{iframe}.node{iNode}.In)
        isNode=bioTree{iframe}.node{iNode}.In{iIn}.isNode;
        if isNode==1
            preNode=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[iframe,iNode,iIn];
        end
        if isNode==0
            preRoot=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[iframe,iNode,iIn];
            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
        end
    end
    for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
        is2Node=bioTree{iframe}.node{iNode}.Out{iOut}.is2Node;
        if is2Node==1
            nextNode=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[iframe,iNode,iOut];
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.leafInfo=[];
        end
        if is2Node==0
            nextLeaf=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[iframe,iNode,iOut];
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.leafInfo=[];
        end
    end
end
end
function bioTree=correctForMistake(bioTree)
a=[];
for iframe=1:size(bioTree,2)
    allNode=[];
    for iNode=1:size(bioTree{iframe}.node,2)
        if strcmp(bioTree{iframe}.node{iNode}.state,'mistake')
            disp(iframe)
            [bioTree,nodeList]=correctMistakeNode(bioTree,[iframe,iNode]);
            allNode=[allNode;nodeList];
        end
    end
    bioTree=bioTreeDeleteNode(bioTree,allNode);
    a=[a;allNode];
end
end
function [nextNodeList,rootList,leaviesList]=detailNodeInfo(bioTree,nodeInfo)
nextNodeList=[];
rootList=[];
leaviesList=[];
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==0
        rootList=[rootList;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo];
    end
end
for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==1
        nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
        if  isempty(nextNodeList) || max(nextNodeList(:,1)==nextNode(1) & nextNodeList(:,2)==nextNode(2))==0
            nextNodeList=[nextNodeList;[nextNode(1),nextNode(2)]];
        end
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==0
        leaviesList=[leaviesList;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo];
    end
end
end
function [bioTree,nextNodeList]=correctMistakeNode(bioTree,nodeInfo)
nextNodeList=[];
inOutInfo=[size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2),size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)];
if min(inOutInfo)>=2  %% n-n-n node
    nextNodeList=nodeInfo;
    rootList=[];
    leafList=[];
    [nextNode,~,leaf]=detailNodeInfo(bioTree,nodeInfo);
    if ~(size(nextNode,1)==1 && isempty(leaf)) 
        if ~(isempty(nextNode) && all(leaf(:,1)==numel(bioTree)))     %% leaf is the final page(no way to segration--must match)       
            nextNodeList=[];
            return
        end
    end
    if ~isempty(nextNode)
        leafList=[leafList;leaf];
        while strcmp(bioTree{nextNode(1)}.node{nextNode(2)}.state,'mistake')
            nextNodeList=[nextNodeList;nextNode(1:2)];
            [nextNode,root,leaf]=detailNodeInfo(bioTree,nextNode);
            rootList=[rootList;root];
            if size(nextNode,1)~=1
                break
            end
        end
    end
    inOutInfoFinal=[size(bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.In,2),size(bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.Out,2)];
    if isempty(leafList) && isempty(rootList) && inOutInfo(1)==inOutInfoFinal(2) 
        pixelIdxListIn=getInputMask(bioTree,nextNodeList(1,:));
        [newTrace,isN2nNode]=n2nNodeMatch(pixelIdxListIn,bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.Out,bioTree{1}.imageSize,40);
        if isN2nNode==1
            bioTree=matchWithNewTrace(bioTree,newTrace,nextNodeList,pixelIdxListIn,inOutInfo(1));
            return
        end
    end
    if isempty(leafList) && isempty(rootList) && inOutInfo(1)==2 && inOutInfoFinal(2)==1  %% n -n-n-1 node (need to find a 1-2 node)
        if (nextNode(1)-nextNodeList(end,1))<=100;
            nextNodeList=[nextNodeList;nextNode];
            inOutInfoFinal=[size(bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.In,2),size(bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.Out,2)];
            if inOutInfoFinal(2)==2 && inOutInfoFinal(1)==1
                pixelIdxListIn=getInputMask(bioTree,nextNodeList(1,:));
                [newTrace,isN2nNode]=n2nNodeMatch(pixelIdxListIn,bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.Out,bioTree{1}.imageSize,40);
                if isN2nNode==1
                    bioTree=matchWithNewTrace(bioTree,newTrace,nextNodeList,pixelIdxListIn,inOutInfo(1));
                    return
                end
            end
        end
    end
    nextNodeList=[];
end
if inOutInfo(1)==2 && inOutInfo(2)==1  % 2-1-2  type (need to find 1-2 node)
    nextNodeList=nodeInfo;
    [nextNode,~,~]=detailNodeInfo(bioTree,nodeInfo);
    nextNodeList=[nextNodeList;nextNode];
    inOutInfoFinal=[size(bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.In,2),size(bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.Out,2)];
    if inOutInfoFinal(1)==1 && inOutInfoFinal(2)==2 && (nextNodeList(end,1)-nextNodeList(1,1))<=100
        pixelIdxListIn=getInputMask(bioTree,nextNodeList(1,:));
        [newTrace,isN2nNode]=n2nNodeMatch(pixelIdxListIn,bioTree{nextNodeList(end,1)}.node{nextNodeList(end,2)}.Out,bioTree{1}.imageSize,40);
        if isN2nNode==1
            bioTree=matchWithNewTrace(bioTree,newTrace,nextNodeList,pixelIdxListIn,inOutInfo(1));
            return
        end
    end
    nextNodeList=[];
end  
end
function bioTree=matchWithNewTrace(bioTree,newTrace,nextNodeList,pixelIdxListIn,num)
for iIn=1:num
    missTrace=[];
    for iMissTrace=1:(nextNodeList(end,1)-nextNodeList(1,1))
        missTrace{1,iMissTrace}=pixelIdxListIn{iIn};
    end
    if bioTree{nextNodeList(1,1)}.node{nextNodeList(1,2)}.In{iIn}.isNode==1
        preNode=bioTree{nextNodeList(1,1)}.node{nextNodeList(1,2)}.In{iIn}.nodeInfo;
        is2Node=newTrace{iIn}.is2Node;
        if is2Node==1
            aimNode=newTrace{iIn}.nodeInfo;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=aimNode;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,missTrace,newTrace{iIn}.traceInfo.pixelIdxList];
            bioTree{aimNode(1)}.node{aimNode(2)}.In{aimNode(3)}.isNode=1;
            bioTree{aimNode(1)}.node{aimNode(2)}.In{aimNode(3)}.nodeInfo=preNode;
            bioTree{aimNode(1)}.node{aimNode(2)}.In{aimNode(3)}.rootInfo=[];
        end
        if is2Node==0;
            aimLeaf=newTrace{iIn}.leafInfo;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=0;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[];
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=aimLeaf;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,missTrace,newTrace{iIn}.traceInfo.pixelIdxList];
            bioTree{aimLeaf(1)}.leavies{aimLeaf(2)}.is2Node=1;
            bioTree{aimLeaf(1)}.leavies{aimLeaf(2)}.nodeInfo=preNode;
            bioTree{aimLeaf(1)}.leavies{aimLeaf(2)}.rootInfo=[];
        end
    end
    if bioTree{nextNodeList(1,1)}.node{nextNodeList(1,2)}.In{iIn}.isNode==0
        preRoot=bioTree{nextNodeList(1,1)}.node{nextNodeList(1,2)}.In{iIn}.rootInfo;
        is2Node=newTrace{iIn}.is2Node;
        if is2Node==1
            aimNode=newTrace{iIn}.nodeInfo;
            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=aimNode;
            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,missTrace,newTrace{iIn}.traceInfo.pixelIdxList];
            bioTree{aimNode(1)}.node{aimNode(2)}.In{aimNode(3)}.isNode=0;
            bioTree{aimNode(1)}.node{aimNode(2)}.In{aimNode(3)}.nodeInfo=[];
            bioTree{aimNode(1)}.node{aimNode(2)}.In{aimNode(3)}.rootInfo=preRoot;
        end
        if is2Node==0
            aimLeaf=newTrace{iIn}.leafInfo;
            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=0;
            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[];
            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=aimLeaf;
            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,missTrace,newTrace{iIn}.traceInfo.pixelIdxList];
            bioTree{aimLeaf(1)}.leavies{aimLeaf(2)}.is2Node=0;
            bioTree{aimLeaf(1)}.leavies{aimLeaf(2)}.nodeInfo=[];
            bioTree{aimLeaf(1)}.leavies{aimLeaf(2)}.rootInfo=preRoot;
        end
    end
end
end
function bioTree=bioTreeDeleteNode(bioTree,allNode)
for i=1:size(allNode,1)
    nodeInfo=allNode(i,:);
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}=[];
end
for i=1:size(allNode,1)
    frameNum=allNode(i,1);
    bioTree=bioTreeLinkCorrection(bioTree,frameNum);
end
end










