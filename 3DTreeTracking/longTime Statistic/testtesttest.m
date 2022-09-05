allResult.detachingResult=[];
n=0;
bacPieceInfo=[];
leafPieceInfo=[];
for iframe=1:size(bioTree,2)-1
    n=n+1;
    if n==201
        n=1;
        allResult.detachingResult=[allResult.detachingResult;mean(bacPieceInfo),sum(leafPieceInfo)/mean(bacPieceInfo)];
        bacPieceInfo=[];
        leafPieceInfo=[];
    end
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    bacNum=size(bacteriaInfo,1);
    leafNum=size(bioTree{iframe}.leavies,2);
    for iLeaf=1:leafNum
        leafInfo=bioTree{iframe}.leavies{iLeaf}.leaviesPixelDetail;
        if any(ismember(leafInfo,limitRegionIdx))
            leafNum=leafNum-1;
        end
    end
    bacPieceInfo=[bacPieceInfo;bacNum];
    leafPieceInfo=[leafPieceInfo;leafNum];
end