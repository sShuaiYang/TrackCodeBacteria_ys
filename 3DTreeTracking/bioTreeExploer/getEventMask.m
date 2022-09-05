function [rootMask,isRootMask,leafMask, isLeafMask,divisionMask,isDivisionMask]=getEventMask(bioTree,allList,frame,stepFrame)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
rootMask=false(xSize,ySize);
leafMask=false(xSize,ySize);
divisionMask=false(xSize,ySize);
isRootMask=false;
isLeafMask=false;
isDivisionMask=false;
allRoot=allList.allRoot;
allLeaf=allList.allLeaf;
allNode=allList.allNode;
coreBranch=[3,4,6];   %% add by jzy map 9.10
for iroot=1:size(allRoot,1)
    rootInfo=allRoot(iroot,1:2);
    if ~ismember(bioTree{rootInfo(1)}.root{rootInfo(2)}.branchIndex,coreBranch) %add by jzy 9.10
        continue                                                                %
    end                                                                         %
    if rootInfo(1)>frame
        continue;
    end
    if rootInfo(1)<=frame &&  rootInfo(1)> frame-stepFrame
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==1
            nextnodeInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
            if nextnodeInfo(1)<=frame
                continue;
            end
            if nextnodeInfo(1)>frame
                pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{frame-rootInfo(1)+1};
                rootMask(pixelIdxList)=1;
                isRootMask=true;
            end
        end
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==0
            leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
            if leafInfo(1)<=frame
                continue;
            end
            if leafInfo(1)>frame
                pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{frame-rootInfo(1)+1};
                rootMask(pixelIdxList)=1;
                isRootMask=true;
            end
        end
    end
end
for iNode=1:size(allNode,1)
    nodeInfo=allNode(iNode,1:2);
    if ~ismember(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex,coreBranch)  % add by jzy 9.10
        continue                                                                 %
    end                                                                          %
    if allNode(iNode,4)==false
        continue;
    end
    if nodeInfo(1)>=frame && nodeInfo(1)<frame+stepFrame
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==true
            parentNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
            pixelIdxList=bioTree{parentNodeInfo(1)}.node{ parentNodeInfo(2)}.Out{parentNodeInfo(3)}.traceInfo.pixelIdxList{end};
            divisionMask(pixelIdxList)=1;
            isDivisionMask=true;
        end
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==false
            rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo;
            pixelIdxList=bioTree{ rootInfo(1)}.root{ rootInfo(2)}.traceInfo.pixelIdxList{end};
            divisionMask(pixelIdxList)=1;
            isDivisionMask=true;
        end
    end
    if nodeInfo(1)<=frame &&  nodeInfo(1)> frame-stepFrame
        %         for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
        %             if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==true
        %                 nextnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
        %                 if nextnodeInfo(1)<=frame
        %                     continue;
        %                 end
        %                 if nextnodeInfo(1)>frame
        %                     pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{frame-nodeInfo(1)+1};
        %                     divisionMask(pixelIdxList)=1;
        %                     isDivisionMask=true;
        %                 end
        %             end
        %         end
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
            divisionMask(pixelIdxList)=1;
            isDivisionMask=true;
        end
    end
end
for iLeaf=1:size(allLeaf,1)
    leafInfo=allLeaf(iLeaf,1:2);
    if ~ismember(bioTree{leafInfo(1)}.leavies{leafInfo(2)}.branchIndex,coreBranch)  % add by jzy 9.10
        continue                                                                 %
    end                                                                          %
    if leafInfo(1)>frame+stepFrame
        continue;
    end
    if leafInfo(1)<frame+stepFrame &&  leafInfo(1)>=frame && leafInfo(1)~=size(bioTree,2)
%         nextnodeInfo=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo;
%         if nextnodeInfo(1)>frame
%             continue;
%         end
%         pixelIdxList=bioTree{nextnodeInfo(1)}.node{nextnodeInfo(2)}.Out{nextnodeInfo(3)}.traceInfo.pixelIdxList{frame-nextnodeInfo(1)+1};
        pixelIdxList=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail;
        leafMask(pixelIdxList)=1;
        isLeafMask=true;
    end
end
if isRootMask==true
    rootMask=imfill(rootMask,'holes');
end
if isLeafMask==true
    leafMask=imfill(leafMask,'holes');
end
if isDivisionMask==true
    divisionMask=imfill(divisionMask,'holes');
end
end