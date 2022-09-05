function [branchMask,isMask]=getBranchMask(bioTree,branchList,branchIndex,frame)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
branchMask=false(xSize,ySize);
isMask=false;
branchInfo=branchList(branchIndex,1:3);
if branchInfo(3)==1
    allRoot=bioTree{branchInfo(1)}.node{branchInfo(2)}.allRoot;
    allNode=bioTree{branchInfo(1)}.node{branchInfo(2)}.allNode;
end
if branchInfo(3)==0
    allRoot=bioTree{branchInfo(1)}.root{branchInfo(2)}.allRoot;
    allNode=[];
end
for iroot=1:size(allRoot,1)
    rootInfo=allRoot(iroot,1:2);
    if rootInfo(1)>frame
        continue;
    end
    if rootInfo(1)<=frame
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==1
            nextnodeInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
            if nextnodeInfo(1)<=frame
                continue;
            end
            if nextnodeInfo(1)>frame
                pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{frame-rootInfo(1)+1};
                branchMask(pixelIdxList)=1;
                isMask=true;
            end
        end
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==0
            leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
            if leafInfo(1)<=frame
                continue;
            end
            if leafInfo(1)>frame
                pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{frame-rootInfo(1)+1};
                branchMask(pixelIdxList)=1;
                isMask=true;
            end
        end
    end
end
for iNode=1:size(allNode,1)
    nodeInfo=allNode(iNode,1:2);
    if nodeInfo(1)>frame
        continue;
    end
    if nodeInfo(1)<=frame
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==true
                nextnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                if nextnodeInfo(1)<=frame
                    continue;
                end
                if nextnodeInfo(1)>frame
                    pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{frame-nodeInfo(1)+1};
                    branchMask(pixelIdxList)=1;
                    isMask=true;
                end
            end
            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==false
                leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                if leafInfo(1)<=frame
                    continue;
                end
                if leafInfo(1)>frame
                    pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{frame-nodeInfo(1)+1};
                    branchMask(pixelIdxList)=1;
                    isMask=true;
                end
            end
        end
    end
end
% if isMask==true
%     branchMask=imfill(branchMask,'holes');
% end
end