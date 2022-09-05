function bioTree=generateBioTree(oriTree,frameNum,imageSize)
for iframe=1:frameNum
    bioTree{iframe}.root=[];
    bioTree{iframe}.node=[];
    bioTree{iframe}.leavies=[];
end
bioTree{1}.imageSize=imageSize;
for iNum=1:size(oriTree,2)
    if ~isempty(oriTree{iNum}.isRoot)
        iframe=oriTree{iNum}.beginFrame;
        iRoot=size(bioTree{iframe}.root,2)+1;
        bioTree{iframe}.root{iRoot}.rootPixelDetail=oriTree{iNum}.pixelIdxList{1};
        bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList=oriTree{iNum}.pixelIdxList';
        if ~isempty(oriTree{iNum}.is2Node)
            bioTree{iframe}.root{iRoot}.is2Node=1;
            bioTree{iframe}.root{iRoot}.leafInfo=[];
            nodeFrame=size(oriTree{iNum}.pixelIdxList,1)+iframe;
            bioTree{iframe}.root{iRoot}.nodeInfo=[nodeFrame,size(bioTree{nodeFrame}.node,2)+1,1];
            nodeInfo=bioTree{iframe}.root{iRoot}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode=0;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo=[iframe,iRoot];
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo=[];
            for i=1:2
                oriTree{oriTree{iNum}.nodeInfo(i)}.bioTreeNode=[nodeInfo(1),nodeInfo(2),i];
            end
        end
        if isempty(oriTree{iNum}.is2Node)
            bioTree{iframe}.root{iRoot}.is2Node=0;
            bioTree{iframe}.root{iRoot}.nodeInfo=[];
            if ~isempty(oriTree{iNum}.hasEnd)
                leafFrame=size(oriTree{iNum}.pixelIdxList,1)+iframe-1;
            else
                leafFrame=frameNum;
            end
            bioTree{iframe}.root{iRoot}.leafInfo=[leafFrame,size(bioTree{leafFrame}.leavies,2)+1];
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=0;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iRoot];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=oriTree{iNum}.pixelIdxList{end};
        end
    end
    if isempty(oriTree{iNum}.isRoot)
        bioTreeNode=oriTree{iNum}.bioTreeNode;
        bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.traceInfo.pixelIdxList=oriTree{iNum}.pixelIdxList';
        if ~isempty(oriTree{iNum}.is2Node)
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.is2Node=1;
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.leafInfo=[];
            nodeFrame=size(oriTree{iNum}.pixelIdxList,1)+bioTreeNode(1);
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.nodeInfo=[nodeFrame,size(bioTree{nodeFrame}.node,2)+1,1];
            nodeInfo=bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode=1;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo=[];
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo=bioTreeNode;
            for i=1:2
                oriTree{oriTree{iNum}.nodeInfo(i)}.bioTreeNode=[nodeInfo(1),nodeInfo(2),i];
            end
        end
        if isempty(oriTree{iNum}.is2Node)
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.is2Node=0;
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.nodeInfo=[];
            if ~isempty(oriTree{iNum}.hasEnd)
                leafFrame=size(oriTree{iNum}.pixelIdxList,1)+bioTreeNode(1)-1;
            else
                leafFrame=frameNum;
            end
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.leafInfo=[leafFrame,size(bioTree{leafFrame}.leavies,2)+1];
            leafInfo=bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=1;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=bioTreeNode;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=oriTree{iNum}.pixelIdxList{end};
        end
    end
end
end