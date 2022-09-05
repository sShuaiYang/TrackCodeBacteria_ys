function [newTrace,canDivideorNot,bioTree,bacteriaList]=fullNodeType_finalTracking(bioTree,nodeInfo)%further track the N input and one Output node, return N new complete bacterial trace
imageSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
for iIn=1:size(pixelIdxListIn,2)
    newTrace{iIn}.traceInfo.pixelIdxList=[];
    newTrace{iIn}.is2Node=[];
end
bacteriaList=zeros(size(newTrace,2),1);
for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    [pixelIdxListNew,canDivideorNot]=basicDivideSolution_deep(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1},imageSize);
    if canDivideorNot==1
        for ibranch=1:size(pixelIdxListNew,2)
            if ~isempty(pixelIdxListNew{ibranch})
                bacteriaList(ibranch)=iOut;
            end
        end   
    else
        bacteriaList=[];
        return
    end
end

newRootInfo=0;
newRoot=0;
for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    if max(bacteriaList==iOut)==0
        newRoot=newRoot+1;
        bacteriaList(end+1)=iOut;
        newTrace{1,end+1}.traceInfo.pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList;
        newTrace{1,end}.is2Node=[];
        if newRoot==1
            newRootInfo=1+nodeInfo(1);
        else
            newRootInfo=[newRootInfo;1+nodeInfo(1)];
        end
    end
end

for iIn=1:size(bacteriaList,1)
    if bacteriaList(iIn)==0
        newTrace{iIn}.is2Node=false;
        newTrace{iIn}.leafInfo=[nodeInfo(1),size(bioTree{nodeInfo(1)}.leavies,2)+1];
    end
end

for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    testis2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
    if testis2Node==true
        testnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
        inNum=size(bioTree{testnodeInfo(1)}.node{testnodeInfo(2)}.In,2);
    end
    realIn=0;
    for iIn=1:size(newTrace,2)
        if isempty(newTrace{iIn}.is2Node)
            if bacteriaList(iIn)==iOut
                newTrace{iIn}.is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
                if newTrace{iIn}.is2Node==true
                    newTrace{iIn}.nodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                    realIn=realIn+1;
                    if realIn~=1
                        newTrace{iIn}.nodeInfo(1,3)=inNum+realIn-1;
                        newNodeInfo=newTrace{iIn}.nodeInfo;
                    end
                end
                if newTrace{iIn}.is2Node==false
                    newTrace{iIn}.leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                end
            end
        else
            if isempty(newTrace{iIn}.traceInfo.pixelIdxList)
                pixelIdxListIn=getInputMask(bioTree,nodeInfo);
                newTrace{iIn}.traceInfo.pixelIdxList{1}=pixelIdxListIn{1,iIn};
            end
        end
    end
end
for iRoot=size(newRootInfo,1)
    if size(newRootInfo,1)==1 && newRootInfo==0
        return
    end
    oriRootNum=size(bioTree{newRootInfo(iRoot)}.root,2);
    bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.is2Node=newTrace{end-size(newRootInfo,1)+iRoot}.is2Node;
    if bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.is2Node==1
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.leafInfo=[];
        nextNodeInfo=newTrace{end-size(newRootInfo,1)+iRoot}.nodeInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.nodeInfo=nextNodeInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.rootPixelDetail=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList{1};
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.traceInfo.pixelIdxList=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList;
        bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{nextNodeInfo(3)}.isNode=0;
        bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{nextNodeInfo(3)}.nodeInfo=[];
        bioTree{nextNodeInfo(1)}.node{nextNodeInfo(2)}.In{nextNodeInfo(3)}.rootInfo=[newRootInfo(iRoot),oriRootNum+1];
    end
    if bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.is2Node==0
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.rootInfo=[];
        nextLeafInfo=newTrace{end-size(newRootInfo,1)+iRoot}.leafInfo;
        if nextLeafInfo(2)==0
            nextLeafInfo(2)=size(bioTree{nextLeafInfo(1)}.root,2)+1;
        end
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.leafInfo=nextLeafInfo;
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.rootPixelDetail=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList{1};
        bioTree{1,newRootInfo(iRoot)}.root{1,oriRootNum+1}.traceInfo.pixelIdxList=newTrace{end-size(newRootInfo,1)+iRoot}.traceInfo.pixelIdxList;
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.isNode=0;
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.nodeInfo=[];
        bioTree{nextLeafInfo(1)}.leavies{nextLeafInfo(2)}.rootInfo=[newRootInfo(iRoot),oriRootNum+1];
    end
end
newTrace(end-size(newRootInfo,1)+1:end)=[];
bacteriaList(end-size(newRootInfo,1)+1:end)=[];
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

%% basicDivideSolution_deep
function [afterDivide,canDivideorNot]=basicDivideSolution_deep(eachInfo,maskImage,pictureSize)
inputNum=numel(eachInfo);
realIn=0;
newRoot=0;
for i=1:inputNum
[~,xyMin,regionImage]=findRegionNum(eachInfo{i},pictureSize);
in{i}.xyMin=xyMin;
in{i}.BWImage=regionImage{1};
if ~isempty(in{i}.BWImage)
    realIn=realIn+1;
end
end
[outputNum,xyMin,regionImage]=findRegionNum(maskImage,pictureSize);
for i=1:outputNum
    out{i}.xyMin=xyMin;
    out{i}.BWImage=regionImage{i};
end
if realIn>outputNum
    for i=1:size(out,2)
        out{1}.BWImage=out{1}.BWImage|out{i}.BWImage;
    end
    [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize);
    memberNum=0;
    if canDivideorNot==1
        for i=1:size(afterDivide,2)
            if ~isempty(afterDivide{i})
                memberNum=memberNum+1;
            end
        end
        if memberNum==0
            canDivideorNot=0;
        end
    end
else
    canDivideorNot=0;
    afterDivide=[];
end
end
%% this  four function is used to divide N to 1 Node,by rotating the model and erode the out image
function [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize)
sita=-30:2:30;
outInfo=out{1};
outImage=outInfo.BWImage;
outImage=imfill(outImage,'holes');
xyMin=outInfo.xyMin;
t=0;
parfor i=1:numel(in)
    if ~isempty(in{i}.BWImage)
    imageStack=rotateStacks(in{i}.BWImage,sita);
    newInImage=gainSameSize(in{i},outInfo,pictureSize);
    imageStack=gainPossibleModel(imageStack,outImage,newInImage);
    inxyMin=in{i}.xyMin;
    afterDivide{i}=findTheBest(in{i}.BWImage,inxyMin,imageStack,pictureSize,xyMin,sita);
    else
        afterDivide{i}=[];
    end
    if ~isempty(afterDivide{i})
        t=t+1;
    end
end
if numel(in)==t;
    canDivideorNot=1;
else
    canDivideorNot=0;
end
end
function imageStack=rotateStacks(BWImage,sita)
% fill in the image and rotate them to stacks
BWImage=imfill(BWImage,'holes');
parfor j=1:numel(sita)
    imageStack{j}=imrotate(BWImage,sita(j));
end
end
function imageStack=gainPossibleModel(imageStack,outImage,inImage)
% get all possible model
parfor j=1:numel(imageStack)
    model=bwmorph(imageStack{j},'erode');
    I=imerode(outImage,model);
    cc=bwconncomp(I);
    if cc.NumObjects==0
        model2=bwmorph(model,'erode');
        I=imerode(outImage,model2);
        cc=bwconncomp(I);
        if cc.NumObjects==1
            imageStack{j}=imdilate(I,model);
            imageStack{j}=imageStack{j}&outImage;
            if max(max(imageStack{j}))==1
                imageStack{j}=bwmorph(imageStack{j},'remove');
            else
                imageStack{j}=[];
            end
        else if cc.NumObjects>=2
                result=0;
                for i=1:cc.NumObjects
                    eachone=zeros(size(I));
                    eachone(cc.PixelIdxList{i})=1;
                    eachone=eachone&inImage;
                    if max(max(eachone))==1
                        result=1;
                        imageStack{j}=imdilate(eachone,model);
                        imageStack{j}=imageStack{j}&outImage;
                        if max(max(imageStack{j}))==1
                            imageStack{j}=bwareaopen(imageStack{j},10);
                            imageStack{j}=bwmorph(imageStack{j},'remove');
                            if max(max(imageStack{j}))==0
                                imageStack{j}=[];
                            end
                        else
                            imageStack{j}=[];
                        end
                        break
                    end
                end
                if result==0
                    imageStack{j}=[];
                end
            end
        end
    end   
    if cc.NumObjects==1
        imageStack{j}=imdilate(I,imageStack{j});
        imageStack{j}=imageStack{j}&outImage;
        if max(max(imageStack{j}))==1
            imageStack{j}=bwmorph(imageStack{j},'remove');
        else
            imageStack{j}=[];
        end
    else if cc.NumObjects>=2
            result=0;
            for i=1:cc.NumObjects
                eachone=zeros(size(I));
                eachone(cc.PixelIdxList{i})=1;
                eachone=eachone&inImage;
                if max(max(eachone))==1
                    result=1;
                    imageStack{j}=imdilate(eachone,imageStack{j});
                    imageStack{j}=imageStack{j}&outImage;
                    if max(max(imageStack{j}))==1
                        imageStack{j}=bwareaopen(imageStack{j},10);
                        imageStack{j}=bwmorph(imageStack{j},'remove');
                        if max(max(imageStack{j}))==0
                            imageStack{j}=[];
                        end
                    else
                        imageStack{j}=[];
                    end
                    break
                end
            end
            if result==0
                imageStack{j}=[];
            end
        end
    end
end
end
function result=findTheBest(inImage,inxyMin,imageStack,pictureSize,xyMin,sita)
% the best choice to the input bacteria
threShold=50;
inInfo=regionprops(inImage,'Centroid','MajorAxisLength','MinorAxisLength','PixelIdx');
inImageInfo=[inInfo.Centroid(2)+inxyMin(1),inInfo.Centroid(1)+inxyMin(2),inInfo.MajorAxisLength,inInfo.MinorAxisLength,0];
outImageInfo=zeros(numel(imageStack),5);
usefulData=zeros(numel(imageStack),1);
for i=1:numel(imageStack)
    if isempty(imageStack{i})
     usefulData(i)=0;
     continue
    else
        usefulData(i)=1;
        outInfo=regionprops(imageStack{i},'Centroid','MajorAxisLength','MinorAxisLength','PixelIdx');
        outImageInfo(i,:)=[outInfo(1,1).Centroid(2)+xyMin(1),outInfo(1,1).Centroid(1)+xyMin(2),outInfo(1,1).MajorAxisLength,outInfo(1,1).MinorAxisLength,outInfo(1,1).MajorAxisLength*sita(i)/180*pi];
    end
end
caculateDist=pdist2(outImageInfo,inImageInfo);
caculateDist(usefulData==0)=max(caculateDist);
[minDist,Num]=min(caculateDist);
if minDist>threShold
    result=[];
else
    correctImage=imageStack{Num};
    result=xy2Idx(xyMin,correctImage,pictureSize);
end
end
%% here are some basic functions
function [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize)
% find how many regions there are in a input/output
if isempty(pixelIdxList)
    regionNum=0;
    xyMin=zeros(1,2);
    regionImage{1}=[];
else
    [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize);
    CC=bwconncomp(BWImage);
    regionNum=CC.NumObjects;
    for i=1:regionNum
        pixelIdxList2=CC.PixelIdxList{i};
        BW=false(CC.ImageSize);
        BW(pixelIdxList2)=1;
        regionImage{i}=BW;
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
function result=AminusRow(A,row)
%used to make a differ between matrixA and a row vector
rowNum=size(A,1);
B=repmat(row,rowNum,1);
result=A-B;
end
function [sumArea,minArea]=getSumArea(info,num)
areaInfo=zeros(num,1);
for i=1:num
    filledImage=imfill(info{i}.BWImage,'holes');
    areaInfo(i)=numel(filledImage(filledImage==1));
end
sumArea=sum(areaInfo);
minArea=min(areaInfo);
end
function newA=gainSameSize(A,B,pictureSize)
A.BWImage=imfill(A.BWImage,'holes');
% change picture in A to the same size that in B
pixelIdxList=xy2Idx(A.xyMin,A.BWImage,pictureSize);
new=false(pictureSize);
new(pixelIdxList)=true;
bSize=size(B.BWImage);
newA=new(B.xyMin(1):B.xyMin(1)+bSize(1)-1,B.xyMin(2):B.xyMin(2)+bSize(2)-1);
end