function bioTree=bioTreePreCut(bioTree,cutNum)
% 前截的程序，从cutNum 的frame起建立新的根
for iframe=cutNum:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==1
                preNode=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
                if preNode(1)<cutNum
                    rootSize=size(bioTree{cutNum}.root,2);
                    bioTree{cutNum}.root{rootSize+1}.is2Node=1;
                    bioTree{cutNum}.root{rootSize+1}.nodeInfo=[iframe,iNode,iIn];
                    bioTree{cutNum}.root{rootSize+1}.leafInfo=[];
                    bioTree{cutNum}.root{rootSize+1}.traceInfo.pixelIdxList=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList(cutNum-preNode(1)+1:end);
                    bioTree{cutNum}.root{rootSize+1}.traceInfo.measurment=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment(cutNum-preNode(1)+1:end);
                    bioTree{iframe}.node{iNode}.In{iIn}.isNode=0;
                    bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo=[];
                    bioTree{iframe}.node{iNode}.In{iIn}.rootInfo=[cutNum,rootSize+1];
                end
            end
            if bioTree{iframe}.node{iNode}.In{iIn}.isNode==0
                preRoot=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
                if preRoot(1)<cutNum
                    rootSize=size(bioTree{cutNum}.root,2);
                    bioTree{cutNum}.root{rootSize+1}.is2Node=1;
                    bioTree{cutNum}.root{rootSize+1}.nodeInfo=[iframe,iNode,iIn];
                    bioTree{cutNum}.root{rootSize+1}.leafInfo=[];
                    bioTree{cutNum}.root{rootSize+1}.traceInfo.pixelIdxList=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList(cutNum-preRoot(1)+1:end);
                    bioTree{cutNum}.root{rootSize+1}.traceInfo.measurment=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment(cutNum-preRoot(1)+1:end);
                    bioTree{iframe}.node{iNode}.In{iIn}.isNode=0;
                    bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo=[];
                    bioTree{iframe}.node{iNode}.In{iIn}.rootInfo=[cutNum,rootSize+1];
                end
            end
        end
    end
    for iLeaf=1:size(bioTree{iframe}.leavies,2)
        if bioTree{iframe}.leavies{iLeaf}.is2Node==1
            preNode=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
            if preNode(1)<cutNum
                rootSize=size(bioTree{cutNum}.root,2);
                bioTree{cutNum}.root{rootSize+1}.is2Node=0;
                bioTree{cutNum}.root{rootSize+1}.nodeInfo=[];
                bioTree{cutNum}.root{rootSize+1}.leafInfo=[iframe,iLeaf];
                bioTree{cutNum}.root{rootSize+1}.traceInfo.pixelIdxList=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList(cutNum-preNode(1)+1:end);
                bioTree{cutNum}.root{rootSize+1}.traceInfo.measurment=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment(cutNum-preNode(1)+1:end);
                bioTree{iframe}.leavies{iLeaf}.is2Node=0;
                bioTree{iframe}.leavies{iLeaf}.nodeInfo=[];
                bioTree{iframe}.leavies{iLeaf}.rootInfo=[cutNum,rootSize+1];
            end
        end
        if bioTree{iframe}.leavies{iLeaf}.is2Node==0
            preRoot=bioTree{iframe}.leavies{iLeaf}.rootInfo;
            if preRoot(1)<cutNum
                rootSize=size(bioTree{cutNum}.root,2);
                bioTree{cutNum}.root{rootSize+1}.is2Node=0;
                bioTree{cutNum}.root{rootSize+1}.nodeInfo=[];
                bioTree{cutNum}.root{rootSize+1}.leafInfo=[iframe,iLeaf];
                bioTree{cutNum}.root{rootSize+1}.traceInfo.pixelIdxList=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList(cutNum-preRoot(1)+1:end);
                bioTree{cutNum}.root{rootSize+1}.traceInfo.measurment=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment(cutNum-preRoot(1)+1:end);
                bioTree{iframe}.leavies{iLeaf}.is2Node=0;
                bioTree{iframe}.leavies{iLeaf}.nodeInfo=[];
                bioTree{iframe}.leavies{iLeaf}.rootInfo=[cutNum,rootSize+1];
            end
        end
    end
end
for iframe=1:cutNum-1
    bioTree{iframe}.node=[];
    bioTree{iframe}.root=[];
    bioTree{iframe}.leavies=[];
end
end