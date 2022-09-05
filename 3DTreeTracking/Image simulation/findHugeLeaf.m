function [info,leafInfo]=findHugeLeaf(bioTree)
imageSize=bioTree{1}.imageSize;
bioTreeEnd=bioTree{end};
parfor iLeaf=1:size(bioTreeEnd.leavies,2)
    pixelIdxList=bioTreeEnd.leavies{iLeaf}.leaviesPixelDetail;
    a=false(imageSize);
    a(pixelIdxList)=1;
    stats=regionprops(a,'FilledArea');
    leafInfo(iLeaf,1)=stats.FilledArea;
end
meanLeafInfo=mean(leafInfo)+50;
info=find(leafInfo>meanLeafInfo);
end
