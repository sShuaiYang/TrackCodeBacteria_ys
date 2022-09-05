function allData=gainAllData(bioTree)
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,~,allList]=divisionFinder(bioTree,branchList);
divisionNode=allList.allNode(allList.allNode(:,4)==1,:);
allDataNum=0;
limitRegionIdx=getLimitRegion(bioTree{1}.imageSize);
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        allDataNum=allDataNum+1;
        allData{allDataNum}.beginInfo=[iframe,iRoot,0];
        allData{allDataNum}.branchIndex=bioTree{iframe}.root{iRoot}.branchIndex;
        if ismember(bioTree{iframe}.root{iRoot}.rootPixelDetail,limitRegionIdx)
            allData{allDataNum}.isAttaching=0;
        else
            allData{allDataNum}.isAttaching=1;
        end
        rootTraceSize=size(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList,2);
        allData{allDataNum}.isDivision=0;
        allData{allDataNum}.frameNum=[iframe,iframe+rootTraceSize];
        allData{allDataNum}.isDetaching=0;
        if bioTree{iframe}.root{iRoot}.is2Node==0
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            if ~(any(ismember(bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail,limitRegionIdx))==1)
                allData{allDataNum}.isDetaching=1;
            end
        end
        allData{allDataNum}.is2Divi=0;
        if bioTree{iframe}.root{iRoot}.is2Node==1
            nodeInfo=bioTree{iframe}.root{iRoot}.nodeInfo;
            for i=1:size(divisionNode,1)
                if isequal(nodeInfo(1:2),divisionNode(i,1:2))
                    allData{allDataNum}.is2Divi=1;
                    break
                end
            end
        end
        if rootTraceSize==0
            allData{allDataNum}.bacteriaInfo=[];
        else
            for iList=1:rootTraceSize
                bacteriaInfo=bioTree{iframe}.root{iRoot}.traceInfo.measurment{iList}(1);
                allData{allDataNum}.bacteriaInfo(iList,:)=[bacteriaInfo.MajorAxisLength,bacteriaInfo.Centroid(2),bacteriaInfo.Centroid(1),bacteriaInfo.Orientation];
            end
        end
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            allDataNum=allDataNum+1;
            allData{allDataNum}.beginInfo=[iframe,iNode,iOut];
            allData{allDataNum}.branchIndex=bioTree{iframe}.node{iNode}.branchIndex;
            allData{allDataNum}.isAttaching=0;
            allData{allDataNum}.isDivision=0;
            for i=1:size(divisionNode,1)
                if isequal([iframe,iNode],divisionNode(i,1:2))
                    allData{allDataNum}.isDivision=1;
                    break
                end
            end
            nodeTraceSize=size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList,2);
            allData{allDataNum}.frameNum=[iframe,iframe+nodeTraceSize];
            allData{allDataNum}.isDetaching=0;
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==0
                leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                if ~(any(ismember(bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail,limitRegionIdx))==1)
                    allData{allDataNum}.isDetaching=1;
                end
            end
            allData{allDataNum}.is2Divi=0;
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==1
                nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                for i=1:size(divisionNode,1)
                    if isequal(nodeInfo(1:2),divisionNode(i,1:2))
                        allData{allDataNum}.is2Divi=1;
                        break
                    end
                end
            end
            if nodeTraceSize==0
                allData{allDataNum}.bacteriaInfo=[];
            else
                for iList=1:nodeTraceSize
                    bacteriaInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iList}(1);
                    allData{allDataNum}.bacteriaInfo(iList,:)=[bacteriaInfo.MajorAxisLength,bacteriaInfo.Centroid(2),bacteriaInfo.Centroid(1),bacteriaInfo.Orientation];
                end
            end
        end
    end
end
end
function limitRegionIdx=getLimitRegion(imageSize)
limitRegion=false(imageSize);
limitRegion(1:3,:)=true;
limitRegion(end-2:end,:)=true;
limitRegion(:,1:3)=true;
limitRegion(:,end-2:end)=true;
limitRegion(1,:)=false;
limitRegion(end,:)=false;
limitRegion(:,1)=false;
limitRegion(:,end)=false;
cc=bwconncomp(limitRegion);
limitRegionIdx=cc.PixelIdxList{1,1};
end