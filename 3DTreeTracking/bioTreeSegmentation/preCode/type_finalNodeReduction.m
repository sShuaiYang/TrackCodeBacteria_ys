function bioTree=type_finalNodeReduction(bioTree)
nodeCount=0;
for iframe=1:size(bioTree,2)
% if iframe==1
%     3498;
% end
    disp(strcat(num2str(iframe),'_'))
    if ~isempty(bioTree{iframe}.node)
        nodeList=type_finalNodeinFrame(bioTree{iframe},iframe,bioTree);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                [traceInfo,canDivideorNot,bioTree,bacteriaList]=fullNodeType_finalTracking(bioTree,nodeList(iList,:));
                nodeList(iList,4)=canDivideorNot;
                if canDivideorNot==1
                bioTree=type1NodeLinker(bioTree,nodeInfoIn,traceInfo,bacteriaList,nodeList(iList,:));
                end
                nodeCount=nodeCount+1;
                disp(nodeCount);
            end
            canDivideNode=nodeList(:,4);
            nodeList(canDivideNode==0,:)=[];
            if ~isempty(nodeList)
            bioTree=removeType1Node(bioTree,nodeList);
            end
        end
    end
end
end
function nodeList=type_finalNodeinFrame(bioTreeFrame,iframe,bioTree)
nodeList=[];
for iNode=1:size(bioTreeFrame.node,2)
    if size(bioTreeFrame.node{iNode}.In,2)>1
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
function bioTree=type1NodeLinker(bioTree,nodeInfoIn,traceInfo,bacteriaList,nodeInfo)
for iIn=1:size(nodeInfoIn,2)
    if nodeInfoIn{iIn}.isNode==false
        rootInfo=nodeInfoIn{iIn}.rootInfo;
        if traceInfo{iIn}.is2Node==false;
            leafInfo=traceInfo{iIn}.leafInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
            if leafInfo(1)~=rootInfo(1)
            bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            end
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
            if leafInfo(1)~=nodeInfo(1)
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
            end
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
function [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize)
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
xyMin=[xMin,yMin];
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