function fluoResult=aveFrameFluo(bioTree,step,fluoImage)
fluoResult=getProperPixel(bioTree,step);
fluoMatrix=zeros(size(fluoResult.fluoPixelMatrix));
for i=1:size(fluoMatrix,2)
    perFluo=fluoImage(:,:,i);
    for j=1:size(fluoMatrix,1)
        pixelValue=fluoResult.fluoPixelMatrix{j,i};
        if isempty(pixelValue)
            continue
        end
        fluoMatrix(j,i)=mean(perFluo(pixelValue{1,1}));
    end
end
fluoResult.fluoMatrix=fluoMatrix;
end
function fluoResult=getProperPixel(bioTree,step)
fluoMatrix=[];
nLine=0;
branchNum=[];
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,~,~]=divisionFinder(bioTree,branchList);
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        traceSize=size(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList,2);
        if traceSize>=2
            beginNum=fix((iframe-1)/step)+1;
            endNum=fix((iframe+traceSize-2)/step)+1;
            if mod(iframe,step)~=1
                if endNum==beginNum
                    continue
                end
                nLine=nLine+1;
                branchNum=[branchNum;bioTree{iframe}.root{iRoot}.branchIndex];
                for i=beginNum+1:endNum
                    fluoMatrix{nLine,i}=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList((i-1)*step+1-iframe+1);
                end
            end
            if mod(iframe,step)==1
                nLine=nLine+1;
                branchNum=[branchNum;bioTree{iframe}.root{iRoot}.branchIndex];
                for i=beginNum:endNum
                    fluoMatrix{nLine,i}=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList((i-1)*step+1-iframe+1);
                end
            end
        end
    end
end
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            traceSize=size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList,2);
            if traceSize>=2
                beginNum=fix((iframe-1)/step)+1;
                endNum=fix((iframe+traceSize-2)/step)+1;
                if mod(iframe,step)~=1
                    if endNum==beginNum
                        continue
                    end
                    nLine=nLine+1;
                    branchNum=[branchNum;bioTree{iframe}.node{iNode}.branchIndex];
                    for i=beginNum+1:endNum
                        fluoMatrix{nLine,i}=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList((i-1)*step+1-iframe+1);
                    end
                end
                if mod(iframe,step)==1
                    nLine=nLine+1;
                    branchNum=[branchNum;bioTree{iframe}.node{iNode}.branchIndex];
                    for i=beginNum:endNum
                        fluoMatrix{nLine,i}=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList((i-1)*step+1-iframe+1);
                    end
                end
            end
        end
    end
end
fluoResult.fluoPixelMatrix=fluoMatrix;
fluoResult.branchNum=branchNum;
end