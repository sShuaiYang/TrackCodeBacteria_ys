function bioTree=bioTreeOptimize(bioTree)
bioTree=autoTidyBioTree(bioTree);
bioTree=bioTreeReduction(bioTree);
% bioTree=detectWrongNode(bioTree);
end
%% delete one frame bacteria(flow from the surface)
function bioTree=bioTreeReduction(bioTree)
% delete the bacteria which only exist for one frame
for iframe=1:size(bioTree,2)
    rootNum=[];
    leafNum=[];
    for iRoot=1:size(bioTree{iframe}.root,2)
        if bioTree{iframe}.root{iRoot}.is2Node==0
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            if leafInfo(1)==iframe
                rootNum=[rootNum;iRoot];
                leafNum=[leafNum;leafInfo(2)];
            end
        end
    end
    bioTree{iframe}.root(rootNum)=[];
    for iRoot=1:size(bioTree{iframe}.root,2)
        if bioTree{iframe}.root{iRoot}.is2Node==0
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iRoot];
        end
        if bioTree{iframe}.root{iRoot}.is2Node==1
            nodeInfo=bioTree{iframe}.root{iRoot}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[iframe,iRoot];
        end
    end
    bioTree{iframe}.leavies(leafNum)=[];
    for iLeaf=1:size(bioTree{iframe}.leavies,2)
        if bioTree{iframe}.leavies{iLeaf}.is2Node==0
            rootInfo=bioTree{iframe}.leavies{iLeaf}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=[iframe,iLeaf];
        end
        if bioTree{iframe}.leavies{iLeaf}.is2Node==1
            nodeInfo=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[iframe,iLeaf];
        end
    end
end
end

%% find node which was produced by wrong image recoginize£¨1 to 2 £©
function bioTree=detectWrongNode(bioTree)
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        if size(bioTree{iframe}.node{iNode}.In,2)~=1 || size(bioTree{iframe}.node{iNode}.Out,2)~=2
            continue
        end
        nodeInfo=[iframe,iNode];
        pixelIdxListIn=getInputMask(bioTree,nodeInfo);
        for i=1:2
            pixelIdxListOut{i}=bioTree{iframe}.node{iNode}.Out{i}.traceInfo.pixelIdxList{1};
        end
        iswrong=iswrongNode(pixelIdxListIn,pixelIdxListOut,bioTree{1}.imageSize);
        if iswrong==1
        end
    end
end
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
        end
    end
end
end
function [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
xyMin=[xMin,yMin];
end
function pixelIdxList=xy2Idx(xyMin,BWImage,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function iswrong=iswrongNode(pixelIdxListIn,pixelIdxListOut,imageSize)
iswrong=0;
[~,inImage]=idx2Xy(pixelIdxListIn{1},imageSize);
inImage=imfill(inImage,'holes');
inImage=bwmorph(inImage,'thin','inf');
inImage=findEndPoint(inImage);
inCC=bwconncomp(inImage);
if inCC.NumObjects~=4
    return
end
for i=1:2
    [~,outImage]=idx2Xy(pixelIdxListOut{i},imageSize);
    outImage=imfill(outImage,'holes');
    outImage=bwmorph(outImage,'thin','inf');
    stats=regionprops(outImage,'Area');
    area(i)=stats.Area;
    outImage=findEndPoint(outImage);
    outCC=bwconncomp(outImage);
    outNum(i)=outCC.NumObjects;
end
if (area(2)>=2*area(1) && outNum(2)==3) ||  (area(1)>=2*area(2) && outNum(1)==3)
    iswrong=1;
end
end
function out=findEndPoint(im)
% endpoints for skel image (bwhitmiss)
se(:,:,1) = [   -1  1 -1   ;...
    -1  1 -1   ;...
    -1 -1 -1   ];

se(:,:,2) = [   -1 -1  1   ;...
    -1  1 -1   ;...
    -1 -1 -1   ];

se(:,:,3) = [   -1 -1 -1   ;...
    -1  1  1   ;...
    -1 -1 -1   ];

se(:,:,4) = [   -1 -1 -1   ;...
    -1  1 -1   ;...
    -1 -1  1   ];

se(:,:,5) = [   -1 -1 -1   ;...
    -1  1 -1   ;...
    -1  1 -1   ];

se(:,:,6) = [   -1 -1 -1   ;...
    -1  1 -1   ;...
    1 -1 -1   ];

se(:,:,7) = [   -1 -1 -1   ;...
    1  1 -1   ;...
    -1 -1 -1   ];

se(:,:,8) = [    1 -1 -1   ;...
    -1  1 -1   ;...
    -1 -1 -1   ];
out = false(size(im));
for i=1:size(se,3)
    hom = bwhitmiss(im, se(:,:,i));
    out = max(out, hom);
end
end