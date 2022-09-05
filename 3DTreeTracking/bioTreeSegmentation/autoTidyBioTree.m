function bioTree=autoTidyBioTree(bioTree)
for iframe=1:size(bioTree,2)
    if iframe==135
        p=1;
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                nodeList=[iframe,iNode,iOut];
                bioTree=tidyBioTreePos(bioTree,nodeList);
            end
        end
    end
end
for iframe=size(bioTree,2):-1:1
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            try
                for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
                    nodeList=[iframe,iNode,iIn];
                    bioTree=tidyBioTreeNeg(bioTree,nodeList);
                end
            catch err
            end
        end
    end
end
end
function bioTree=tidyBioTreePos(bioTree,nodeList)
traceSize=size(bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList,2);
pictureSize=bioTree{1}.imageSize;
if traceSize==0
    return
end
pixelIdxList=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList{1};
[regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize);
if regionNum==1
    return
else
    for i=1:regionNum
        newPixelIdxList{i}=xy2Idx(xyMin,regionImage{i},pictureSize);
    end
end
if traceSize==1
    is2Node=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.is2Node;
    if is2Node==0
        leafInfo=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.leafInfo;
        for i=1:regionNum
            if i==1
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.is2Node=0;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.leafInfo=leafInfo;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.nodeInfo=[];
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList=newPixelIdxList(i);
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=1;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeList;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=newPixelIdxList{i};
            end
            if i>1
                iLeaf=size(bioTree{leafInfo(1)}.leavies,2)+1;
                iOut=size(bioTree{nodeList(1)}.node{nodeList(2)}.Out,2)+1;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.is2Node=0;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.leafInfo=[leafInfo(1),size(bioTree{leafInfo(1)}.leavies,2)+1];
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.nodeInfo=[];
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.traceInfo.pixelIdxList=newPixelIdxList(i);
                bioTree{leafInfo(1)}.leavies{iLeaf}.is2Node=1;
                bioTree{leafInfo(1)}.leavies{iLeaf}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{iLeaf}.nodeInfo=[nodeList(1),nodeList(2),iOut];
                bioTree{leafInfo(1)}.leavies{iLeaf}.leaviesPixelDetail=newPixelIdxList{i};
            end
        end
    end
    if is2Node==1
        nodeInfo=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.nodeInfo;
        for i=1:regionNum
            if i==1
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.is2Node=1;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.nodeInfo=nodeInfo;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.leafInfo=[];
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList=newPixelIdxList(i);
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=nodeList;
            end
            if i>1
                iOut=size(bioTree{nodeList(1)}.node{nodeList(2)}.Out,2)+1;
                iIn=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)+1;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.is2Node=1;
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.leafInfo=[];
                bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.traceInfo.pixelIdxList=newPixelIdxList(i);
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode=1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo=[nodeList(1),nodeList(2),iOut];
            end
        end
    end
end
if traceSize>1
    is2Node=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.is2Node;
    if is2Node==0
        newNode=[nodeList(1)+1,size(bioTree{nodeList(1)+1}.node,2)+1];
        leafInfo=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.leafInfo;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.is2Node=0;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.leafInfo=leafInfo;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.nodeInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.traceInfo.pixelIdxList=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList(2:end);
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=1;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[newNode,1];
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
    end
    if is2Node==1
        newNode=[nodeList(1)+1,size(bioTree{nodeList(1)+1}.node,2)+1];
        nodeInfo=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.nodeInfo;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.is2Node=1;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.leafInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.nodeInfo=nodeInfo;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.traceInfo.pixelIdxList=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList(2:end);
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=1;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[newNode,1];
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[];
    end
    for i=1:regionNum
        if i==1
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.is2Node=1;
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.nodeInfo=[newNode,i];
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.leafInfo=[];
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList=newPixelIdxList(i);
            bioTree{newNode(1)}.node{newNode(2)}.In{i}.isNode=1;
            bioTree{newNode(1)}.node{newNode(2)}.In{i}.rootInfo=[];
            bioTree{newNode(1)}.node{newNode(2)}.In{i}.nodeInfo=nodeList;
        end
        if i>1
            iOut=size(bioTree{nodeList(1)}.node{nodeList(2)}.Out,2)+1;
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.is2Node=1;
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.nodeInfo=[newNode,i];
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.leafInfo=[];
            bioTree{nodeList(1)}.node{nodeList(2)}.Out{iOut}.traceInfo.pixelIdxList=newPixelIdxList(i);
            bioTree{newNode(1)}.node{newNode(2)}.In{i}.isNode=1;
            bioTree{newNode(1)}.node{newNode(2)}.In{i}.rootInfo=[];
            bioTree{newNode(1)}.node{newNode(2)}.In{i}.nodeInfo=[nodeList(1),nodeList(2),iOut];
        end
    end
end
end
function bioTree=tidyBioTreeNeg(bioTree,nodeList)
pictureSize=bioTree{1}.imageSize;
isNode=bioTree{nodeList(1)}.node{nodeList(2)}.In{nodeList(3)}.isNode;
if isNode==0
    preRoot=bioTree{nodeList(1)}.node{nodeList(2)}.In{nodeList(3)}.rootInfo;
    traceSize=size(bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,2);
    if traceSize==0 || traceSize==1
        return
    end
    if traceSize>1
        pixelIdxList=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList{end};
        [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize);
        if regionNum==1
            return
        else
            for i=1:regionNum
                newPixelIdxList{i}=xy2Idx(xyMin,regionImage{i},pictureSize);
            end
        end
        newNode=[nodeList(1)-1,size(bioTree{nodeList(1)-1}.node,2)+1];
        bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
        bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[newNode,1];
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList(1:end-1);
        bioTree{newNode(1)}.node{newNode(2)}.In{1}.isNode=0;
        bioTree{newNode(1)}.node{newNode(2)}.In{1}.rootInfo=preRoot;
        bioTree{newNode(1)}.node{newNode(2)}.In{1}.nodeInfo=[];
    end
end
if isNode==1
    preNode=bioTree{nodeList(1)}.node{nodeList(2)}.In{nodeList(3)}.nodeInfo;
    traceSize=size(bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,2);
    if traceSize==0 || traceSize==1
        return
    end
    if traceSize>1
        pixelIdxList=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList{end};
        [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize);
        if regionNum==1
            return
        else
            for i=1:regionNum
                newPixelIdxList{i}=xy2Idx(xyMin,regionImage{i},pictureSize);
            end
        end
        newNode=[nodeList(1)-1,size(bioTree{nodeList(1)-1}.node,2)+1];
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[newNode,1];
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList(1:end-1);
        bioTree{newNode(1)}.node{newNode(2)}.In{1}.isNode=1;
        bioTree{newNode(1)}.node{newNode(2)}.In{1}.rootInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.In{1}.nodeInfo=preNode;
    end
end
for i=1:regionNum
    if i==1
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.is2Node=1;
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.nodeInfo=nodeList;
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.leafInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.traceInfo.pixelIdxList=newPixelIdxList(i);
        bioTree{nodeList(1)}.node{nodeList(2)}.In{nodeList(3)}.isNode=1;
        bioTree{nodeList(1)}.node{nodeList(2)}.In{nodeList(3)}.nodeInfo=[newNode,1];
        bioTree{nodeList(1)}.node{nodeList(2)}.In{nodeList(3)}.rootInfo=[];
    end
    if i>1
        newInNum=size(bioTree{nodeList(1)}.node{nodeList(2)}.In,2)+1;
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.is2Node=1;
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.nodeInfo=[nodeList(1),nodeList(2),newInNum];
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.leafInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.traceInfo.pixelIdxList=newPixelIdxList(i);
        bioTree{nodeList(1)}.node{nodeList(2)}.In{newInNum}.isNode=1;
        bioTree{nodeList(1)}.node{nodeList(2)}.In{newInNum}.nodeInfo=[newNode,i];
        bioTree{nodeList(1)}.node{nodeList(2)}.In{newInNum}.rootInfo=[];
    end  
end
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
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
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize)
% find how many regions there are in a input/output
if isempty(pixelIdxList)
    regionNum=0;
    xyMin=zeros(1,2);
    regionImage{1}=[];
else
    [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize);
    CC=bwconncomp(BWImage,4);
    regionNum=CC.NumObjects;
    for i=1:regionNum
        pixelIdxList2=CC.PixelIdxList{i};
        BW=false(CC.ImageSize);
        BW(pixelIdxList2)=1;
        regionImage{i}=BW;
    end
end
end