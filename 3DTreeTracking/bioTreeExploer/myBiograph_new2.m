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
function egdeSize=getEdgeSize(EgdeTemp)
fillArea=[];
for i=1:size(EgdeTemp,1)
    fillArea=[fillArea,EgdeTemp(i).FilledArea];
end
egdeSize=sum(fillArea);
end
function [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=pixelIdxList-(yresult-1)*xSize;
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