function bioTree=bioTreeCorrection(bioTree,forbiddenImage)
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
searchRegion=20;
searchThreShold=5;
for iframe=1:size(bioTree,2)
    iLeaf=0;
    while iLeaf~=size(bioTree{iframe}.leavies,2)
        iLeaf=iLeaf+1;
        pixelIdxList=bioTree{iframe}.leavies{iLeaf}.leafPixelDetail;
        %         if all(ismember(pixelIdxList,limitRegionIdx)==0)
        [needLink,rootFrame,iRoot]=searchLink(bioTree,iframe,iLeaf,searchRegion,searchThreShold);
        if needLink==1
            bioTree=bioTreeLink(bioTree,iframe,iLeaf,rootFrame,iRoot);
            bioTree=bioTreeResize(bioTree,iframe,iLeaf,rootFrame,iRoot);
            iLeaf=iLeaf-1;
            %             end
        end
    end
end
bioTree{1}.forbiddenImage=[];
if nargin==2
    % especially for bacteria which would form cluster easily
    bioTree=bioTreeReductionAccordingForbiddenImage(bioTree,forbiddenImage);
    bioTree{1}.forbiddenImage=forbiddenImage;
end
bioTree=bioTreeClearBorder(bioTree);
end
function limitRegionIdx=getLimitRegion(cropInfo)
cropInfo1=bwmorph(cropInfo,'remove');
cropInfo1=imdilate(cropInfo1,true(13));
limitRegion=cropInfo1 & cropInfo;
cc=bwconncomp(limitRegion);
limitRegionIdx=cc.PixelIdxList{1,1};
end
function [needLink,rootFrame,iRoot]=searchLink(bioTree,iframe,iLeaf,searchRegion,searchThreShold)
leafCentroid=bioTree{iframe}.leavies{iLeaf}.leafMeasurment.Centroid;
needLink=0;
searchEnd=iframe+searchRegion;
if searchEnd>size(bioTree,2)
    searchEnd=size(bioTree,2);
end
rootFrame=[];
iRoot=[];
for rootFrame=iframe+1:searchEnd
    for iRoot=1:size(bioTree{rootFrame}.root,2)
        targetCentroid=bioTree{rootFrame}.root{iRoot}.rootMeasurment.Centroid;
        targetDist=(sum((leafCentroid-targetCentroid).^2)).^0.5;
        if targetDist<searchThreShold
            needLink=1;
            return
        end
    end
end
end
function bioTree=bioTreeLink(bioTree,iframe,iLeaf,rootFrame,iRoot)
preIsNode=bioTree{iframe}.leavies{iLeaf}.is2Node;
leafPixel=bioTree{iframe}.leavies{iLeaf}.leafPixelDetail;
leafMeasurment=bioTree{iframe}.leavies{iLeaf}.leafMeasurment;
if preIsNode==0
    preRoot=bioTree{iframe}.leavies{iLeaf}.rootInfo;
    traceIsNode=bioTree{rootFrame}.root{iRoot}.is2Node;
    if traceIsNode==0
        nextLeaf=bioTree{rootFrame}.root{iRoot}.leafInfo;
        bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=0;
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=nextLeaf;
        bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[];
        for i=iframe-preRoot(1)+1:rootFrame-preRoot(1)
            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList{i}=leafPixel;
            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{i}=leafMeasurment;
        end
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList(rootFrame-preRoot(1)+1:nextLeaf(1)-preRoot(1)+1)=bioTree{rootFrame}.root{iRoot}.traceInfo.pixelIdxList;
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment(rootFrame-preRoot(1)+1:nextLeaf(1)-preRoot(1)+1)=bioTree{rootFrame}.root{iRoot}.traceInfo.measurment;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=preRoot;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
    end
    if traceIsNode==1
        nextNode=bioTree{rootFrame}.root{iRoot}.nodeInfo;
        bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
        bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=nextNode;
        for i=iframe-preRoot(1)+1:rootFrame-preRoot(1)
            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList{i}=leafPixel;
            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{i}=leafMeasurment;
        end
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList(rootFrame-preRoot(1)+1:nextNode(1)-preRoot(1))=bioTree{rootFrame}.root{iRoot}.traceInfo.pixelIdxList;
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment(rootFrame-preRoot(1)+1:nextNode(1)-preRoot(1))=bioTree{rootFrame}.root{iRoot}.traceInfo.measurment;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=preRoot;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
    end
end
if preIsNode==1
    preNode=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
    traceIsNode=bioTree{rootFrame}.root{iRoot}.is2Node;
    if traceIsNode==0
        nextLeaf=bioTree{rootFrame}.root{iRoot}.leafInfo;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=0;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=nextLeaf;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[];
        for i=iframe-preNode(1)+1:rootFrame-preNode(1)
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList{i}=leafPixel;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{i}=leafMeasurment;
        end
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList(rootFrame-preNode(1)+1:nextLeaf(1)-preNode(1)+1)=bioTree{rootFrame}.root{iRoot}.traceInfo.pixelIdxList;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment(rootFrame-preNode(1)+1:nextLeaf(1)-preNode(1)+1)=bioTree{rootFrame}.root{iRoot}.traceInfo.measurment;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=preNode;
    end
    if traceIsNode==1
        nextNode=bioTree{rootFrame}.root{iRoot}.nodeInfo;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=nextNode;
        for i=iframe-preNode(1)+1:rootFrame-preNode(1)
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList{i}=leafPixel;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{i}=leafMeasurment;
        end
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList(rootFrame-preNode(1)+1:nextNode(1)-preNode(1))=bioTree{rootFrame}.root{iRoot}.traceInfo.pixelIdxList;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment(rootFrame-preNode(1)+1:nextNode(1)-preNode(1))=bioTree{rootFrame}.root{iRoot}.traceInfo.measurment;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=preNode;
    end
end
end
function bioTree=bioTreeResize(bioTree,iframe,iLeaf,rootFrame,iRoot)
bioTree{iframe}.leavies(iLeaf)=[];
for iLeaf=1:size(bioTree{iframe}.leavies,2)
    is2Node=bioTree{iframe}.leavies{iLeaf}.is2Node;
    if is2Node==1
        preNode=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[iframe,iLeaf];
    end
    if is2Node==0
        preRoot=bioTree{iframe}.leavies{iLeaf}.rootInfo;
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[iframe,iLeaf];
    end
end
bioTree{rootFrame}.root(iRoot)=[];
for iRoot=1:size(bioTree{rootFrame}.root,2)
    is2Node=bioTree{rootFrame}.root{iRoot}.is2Node;
    if is2Node==1
        nextNode=bioTree{rootFrame}.root{iRoot}.nodeInfo;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[rootFrame,iRoot];
    end
    if is2Node==0
        nextLeaf=bioTree{rootFrame}.root{iRoot}.leafInfo;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[rootFrame,iRoot];
    end
end
end
%% especially for cluster 
function bioTree=bioTreeReductionAccordingForbiddenImage(bioTree,forbiddenImage)
cc=bwconncomp(forbiddenImage);
forbiddenRegionIdx=[];
for i=1:cc.NumObjects
    forbiddenRegionIdx=[forbiddenRegionIdx;cc.PixelIdxList{i}];
end
for iframe=1:size(bioTree,2)-1
    iLeaf=0;
    while iLeaf~=size(bioTree{iframe}.leavies,2)
        iLeaf=iLeaf+1;
        pixelIdxList=bioTree{iframe}.leavies{iLeaf}.leafPixelDetail;
        if all(ismember(pixelIdxList,forbiddenRegionIdx)==1)
            [bioTree,iLeaf]=leafForbiddenImage(bioTree,iframe,iLeaf,forbiddenRegionIdx);
        end
    end
end    
end
function [bioTree,iLeaf]=leafForbiddenImage(bioTree,iframe,iLeaf,forbiddenRegionIdx)
preIsNode=bioTree{iframe}.leavies{iLeaf}.is2Node;
if preIsNode==0
    preRoot=bioTree{iframe}.leavies{iLeaf}.rootInfo;
    pixelIdxList=bioTree{preRoot(1)}.root{preRoot(2)}.rootPixelDetail;
    if all(ismember(pixelIdxList,forbiddenRegionIdx)==1) && numel(bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList)<=100
        bioTree=bioTreeResize(bioTree,iframe,iLeaf,preRoot(1),preRoot(2));
        iLeaf=iLeaf-1;
    end
end
if preIsNode==1
    preNode=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
    if size(bioTree{preNode(1)}.node{preNode(2)}.Out,2)>1 && numel(bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList)<=100
        bioTree=bioTreeNodeLeafResize(bioTree,preNode,iframe,iLeaf);
        iLeaf=iLeaf-1;
    end
end
end
function bioTree=bioTreeNodeLeafResize(bioTree,preNode,iframe,iLeaf)
bioTree{iframe}.leavies(iLeaf)=[];
for iLeaf=1:size(bioTree{iframe}.leavies,2)
    is2Node=bioTree{iframe}.leavies{iLeaf}.is2Node;
    if is2Node==1
        leafNode=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
        bioTree{leafNode(1)}.node{leafNode(2)}.Out{leafNode(3)}.leafInfo=[iframe,iLeaf];
    end
    if is2Node==0
        preRoot=bioTree{iframe}.leavies{iLeaf}.rootInfo;
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[iframe,iLeaf];
    end
end
bioTree{preNode(1)}.node{preNode(2)}.Out(preNode(3))=[];
for iOut=1:size(bioTree{preNode(1)}.node{preNode(2)}.Out,2)
    is2Node=bioTree{preNode(1)}.node{preNode(2)}.Out{iOut}.is2Node;
    if is2Node==1
        nextNode=bioTree{preNode(1)}.node{preNode(2)}.Out{iOut}.nodeInfo;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[preNode(1),preNode(2),iOut];
    end
    if is2Node==0
        nextLeaf=bioTree{preNode(1)}.node{preNode(2)}.Out{iOut}.leafInfo;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[preNode(1),preNode(2),iOut];
    end
end
end

%% always border is a confused problem, delete some nonsense bacteria, which would made the cluster ambiguous.
function bioTree=bioTreeClearBorder(bioTree)
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
for iframe=1:size(bioTree,2)
    iRoot=0;
    while iRoot~=size(bioTree{iframe}.root,2)
        iRoot=iRoot+1;
        is2Node=bioTree{iframe}.root{iRoot}.is2Node;
        rootPixelDetail=bioTree{iframe}.root{iRoot}.rootPixelDetail;
%         a=false(2160,2560);
%         a(rootPixelDetail)=1;
%         imshow(a)
%         close all
        if is2Node==0 && size(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList,2)<=100
            if any(ismember(rootPixelDetail,limitRegionIdx))==1
                finalPixelDetail=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{end};
                if any(ismember(finalPixelDetail,limitRegionIdx))==1
                    leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
                    bioTree=bioTreeResize(bioTree,leafInfo(1),leafInfo(2),iframe,iRoot);
                    iRoot=iRoot-1;
                end
            end
        end
    end
end
end