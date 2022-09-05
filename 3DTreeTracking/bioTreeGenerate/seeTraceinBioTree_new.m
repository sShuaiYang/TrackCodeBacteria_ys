function seeTraceinBioTree_new(bioTree,startFrame,is2Node,info)
xSize=1024;
ySize=1024;
if is2Node==0
    rootNum=info(1);
    bwImageStacks=false(xSize,ySize,1:size(bioTree{startFrame}.root{rootNum}.traceInfo.pixelIdxList,2));
    for iframe=1:size(bioTree{startFrame}.root{rootNum}.traceInfo.pixelIdxList,2)
        bwImageTemp=false(xSize,ySize);
        pixelIdxList=bioTree{startFrame}.root{rootNum}.traceInfo.pixelIdxList{iframe};
        bwImageTemp(pixelIdxList)=1;
        bwImageStacks(:,:,iframe)=bwImageTemp;
        clear bwImageTemp;
    end
    implay(bwImageStacks);
end

if is2Node==1
    nodeNum=info(1);
    bwImageStacks=false(xSize,ySize,1:size(bioTree{startFrame}.node{nodeNum}.Out{info(2)}.traceInfo.pixelIdxList,2));
    for iframe=1:size(bioTree{startFrame}.node{nodeNum}.Out{info(2)}.traceInfo.pixelIdxList,2)
        bwImageTemp=false(xSize,ySize);
        pixelIdxList=bioTree{startFrame}.node{nodeNum}.Out{info(2)}.traceInfo.pixelIdxList{iframe};
        bwImageTemp(pixelIdxList)=1;
        bwImageStacks(:,:,iframe)=bwImageTemp;
        clear bwImageTemp;
    end
    implay(bwImageStacks);
end
end