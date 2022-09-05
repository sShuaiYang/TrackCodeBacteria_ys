function [newTrace,inNum,outNum,canDivideorNot]=fullNodeType3Tracking(bioTree,nodeInfo,threShold)%further track the N input and one Output node, return N new complete bacterial trace
imageSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
[newTrace,inNum,outNum,canDivideorNot]=n22Solution(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,imageSize,threShold);
end
function [newTrace,inNum,outNum,canDivideorNot]=n22Solution(pixelIdxListIn,bioTreeOut,pictureSize,threShold)
inputNum=numel(pixelIdxListIn);
outputNum=size(bioTreeOut,2);
for i=1:inputNum
[regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListIn{i},pictureSize);
if regionNum>1
    if ~((inputNum==3 && outputNum==3) || (inputNum==3 && outputNum==4) || (inputNum==4 && outputNum==5)|| inputNum==2)
        newTrace=[];
        inNum=[];
        outNum=[];
        canDivideorNot=0;
        return
    end
   in{i}.xyMin=[];
   in{i}.BWImage=[];
else
in{i}.xyMin=xyMin;
in{i}.BWImage=regionImage{1};
end
end
realOut=0;
for i=1:outputNum
[regionNum,xyMin,regionImage]=findRegionNum(bioTreeOut{1,i}.traceInfo.pixelIdxList{1},pictureSize);
if regionNum>1
    out{i}.xyMin=[];
    out{i}.BWImage=[];
else
    realOut=realOut+1;
    out{i}.xyMin=xyMin;
    out{i}.BWImage=regionImage{1};
end
end
if realOut==0;
    newTrace=[];
    inNum=[];
    outNum=[];
    canDivideorNot=0;
    return
end
if inputNum==4 && outputNum==3 || inputNum==3 && outputNum==4 || inputNum==3 && outputNum==3
    [inNum,outNum,isbreak]=four23NodeSolution(in,out,threShold);
    if isbreak==1
        newTrace=[];
        inNum=[];
        outNum=[];
        canDivideorNot=0;
        return
    end
else
    if inputNum==4 && outputNum==2
        [inNum,outNum,isbreak]=four22NodeSolution(in,out,threShold);
        if isbreak==1
            newTrace=[];
            inNum=[];
            outNum=[];
            canDivideorNot=0;
            return
        end
    else
        if inputNum==4 && outputNum==5 || inputNum==5 && outputNum==4
            [inNum,outNum,isbreak]=five24NodeSolution(in,out,threShold);
            if isbreak==1
                newTrace=[];
                inNum=[];
                outNum=[];
                canDivideorNot=0;
                return
            end
        else
            [inNum,outNum,isBreak]=n22NodeSolution(in,out,threShold);
            if numel(inNum)>1 ...
                    || inputNum==2 && realOut==2 && outputNum==2 && inNum==0 ...
                    || isBreak==1 ...
                    || (inputNum==2 && outputNum>=3 && outputNum~=realOut)
                newTrace=[];
                inNum=[];
                outNum=[];
                canDivideorNot=0;
                return
            end
        end
    end
end
for i=1:numel(outNum)
    newTrace{i}.is2Node=bioTreeOut{outNum(i)}.is2Node;
    newTrace{i}.traceInfo.pixelIdxList=bioTreeOut{outNum(i)}.traceInfo.pixelIdxList;
    if newTrace{i}.is2Node==true
        newTrace{i}.nodeInfo=bioTreeOut{outNum(i)}.nodeInfo;
    end
    if newTrace{i}.is2Node==false
        newTrace{i}.leafInfo=bioTreeOut{outNum(i)}.leafInfo;
    end
    canDivideorNot=1;
end
end

function [inputNum,outputNum,isBreak]=n22NodeSolution(in,out,threShold)
% used to deal with N-to-M node and divide them,here N>=M
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
    Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,1:4)=0;
    end
end
outNum=numel(out);
for i=1:outNum;
    if ~isempty(out{i}.BWImage)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    outArea(i,:)=Info.FilledArea;
    else
    outInfo(i,1:4)=inf;
    end
end
distanceMat=zeros(inNum,outNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
%     inputNum=0;
%     [~,outputNum]=min(outArea);
isBreak=1;
inputNum=[];
outputNum=[];
else
    [inputNum,outputNum]=find(distanceMat==minDist);
    isBreak=0;
end
end

function [inputNum,outputNum,isbreak]=four23NodeSolution(in,out,threShold)
% used to deal with N-to-M node and divide them,here N>=M
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
    Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,1:4)=0;
    end
end
outNum=numel(out);
for i=1:outNum;
    if ~isempty(out{i}.BWImage)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    outArea(i,:)=Info.FilledArea;
    else
    outInfo(i,1:4)=inf;
    end
end
distanceMat=zeros(inNum,outNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
else
    [iIn,iOut]=find(distanceMat==minDist);
    if numel(iIn)==1
        inputNum(1)=iIn;
        outputNum(1)=iOut;
    else
        isbreak=1;
        inputNum=[];
        outputNum=[];
        return
    end
end
outInfo(outputNum(1),:)=inf;
inInfo(inputNum(1),:)=0;
distanceMat=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
else
    [iIn,iOut]=find(distanceMat==minDist);
    if numel(iIn)==1
        inputNum(2)=iIn;
        outputNum(2)=iOut;
        isbreak=0;
    else
        isbreak=1;
        inputNum=[];
        outputNum=[];
        return
    end
end
end
function [inputNum,outputNum,isbreak]=four22NodeSolution(in,out,threShold)
% used to deal with N-to-M node and divide them,here N>=M
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
    Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    inArea(i,:)=Info.FilledArea;
    else
        inInfo(i,1:4)=0;
    end
end
outNum=numel(out);
for i=1:outNum;
    if ~isempty(out{i}.BWImage)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    outArea(i,:)=Info.FilledArea;
    else
    outInfo(i,1:4)=inf;
    end
end
distanceMat=zeros(inNum,outNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
else
    [iIn,iOut]=find(distanceMat==minDist);
    if numel(iIn)==1
        if abs(sum(inArea)-sum(outArea))<150 && abs(inArea(iIn)-outArea(iOut))<50
        inputNum(1)=iIn;
        outputNum(1)=iOut;
        isbreak=0;
        return
        end
    end
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
end
end
function [inputNum,outputNum,isbreak]=five24NodeSolution(in,out,threShold)
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
    Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,1:4)=0;
    end
end
outNum=numel(out);
for i=1:outNum;
    if ~isempty(out{i}.BWImage)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','FilledArea');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    outArea(i,:)=Info.FilledArea;
    else
    outInfo(i,1:4)=inf;
    end
end
distanceMat=zeros(inNum,outNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
else
    [iIn,iOut]=find(distanceMat==minDist);
    if numel(iIn)==1
        inputNum(1)=iIn;
        outputNum(1)=iOut;
    else
        isbreak=1;
        inputNum=[];
        outputNum=[];
        return
    end
end
outInfo(outputNum(1),:)=inf;
inInfo(inputNum(1),:)=0;
distanceMat=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
else
    [iIn,iOut]=find(distanceMat==minDist);
    if numel(iIn)==1
        inputNum(2)=iIn;
        outputNum(2)=iOut;
    else
        isbreak=1;
        inputNum=[];
        outputNum=[];
        return
    end
end
outInfo(outputNum(2),:)=inf;
inInfo(inputNum(2),:)=0;
distanceMat=pdist2(inInfo,outInfo);
minDist=min(min(distanceMat));
if minDist>threShold
    isbreak=1;
    inputNum=[];
    outputNum=[];
    return
else
    [iIn,iOut]=find(distanceMat==minDist);
    if numel(iIn)==1
        inputNum(3)=iIn;
        outputNum(3)=iOut;
        isbreak=0;
    else
        isbreak=1;
        inputNum=[];
        outputNum=[];
        return
    end
end
end
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
% BWImageGain=imfill(BWImageGain,'holes');
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