function [newTrace,isN2nNode]=fullNodeType2Tracking(bioTree,nodeInfo,strength)%further track the N input and one Output node, return N new complete bacterial trace
imageSize=bioTree{1}.imageSize;
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
for iIn=1:size(pixelIdxListIn,2)
    newTrace{iIn}.traceInfo.pixelIdxList=[];
    newTrace{iIn}.is2Node=[];
end
[newTrace,isN2nNode]=n2nSolution(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,imageSize,strength);
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
function [newTrace,isN2nNode]=n2nSolution(pixelIdxListIn,bioTreeOut,pictureSize,strength)
inputNum=numel(pixelIdxListIn);
for i=1:inputNum
[regionNum,xyMin,regionImage]=findRegionNum(pixelIdxListIn{i},pictureSize);
if regionNum>1
    isN2nNode=0;
    newTrace=[];
    return
else
in{i}.xyMin=xyMin;
in{i}.BWImage=regionImage{1};
end
end
outputNum=size(bioTreeOut,2);
for i=1:outputNum
[regionNum,xyMin,regionImage]=findRegionNum(bioTreeOut{1,i}.traceInfo.pixelIdxList{1},pictureSize);
if regionNum>1
    isN2nNode=0;
    newTrace=[];
    return
else
out{i}.xyMin=xyMin;
out{i}.BWImage=regionImage{1};
end
end
[afterDivide,isN2nNode,largeDist]=n2nNodeSolution(in,out,strength);
if isN2nNode==0
    newTrace=[];
else
    needChange=0;
    if largeDist==1 && strength==2
        for iOut=1:size(bioTreeOut,2)
            if size(bioTreeOut{i}.traceInfo.pixelIdxList,2)>1
                break
            else if iOut==size(bioTreeOut,2)
                    needChange=1;
                end
            end
        end
    end
    for i=1:numel(afterDivide)
        newTrace{i}.is2Node=bioTreeOut{afterDivide(i)}.is2Node;
        newTrace{i}.traceInfo.pixelIdxList=bioTreeOut{afterDivide(i)}.traceInfo.pixelIdxList;
        if needChange==1
            newTrace{i}.traceInfo.pixelIdxList(1)=pixelIdxListIn(i);
        end
        if newTrace{i}.is2Node==true
            newTrace{i}.nodeInfo=bioTreeOut{afterDivide(i)}.nodeInfo;
        end
        if newTrace{i}.is2Node==false
            newTrace{i}.leafInfo=bioTreeOut{afterDivide(i)}.leafInfo;
        end
    end
end
end
function [afterDivide,isN2nNode,largeDist]=n2nNodeSolution(in,out,strength)
% used to deal with N-to-M node and divide them,here N>=M
inNum=numel(in);
if strength==1
    threShold=15;
else
    threShold=40;
end
for i=1:inNum
    if ~isempty(in{i}.BWImage)
        Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength,Info.Orientation/360*2*pi];
    else
        inInfo(i,:)=zeros(1,4);
    end
end
outNum=numel(out);
for i=1:numel(out)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength,Info.Orientation/360*2*pi];
end
distanceMat=zeros(inNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo(:,1:4),outInfo(:,1:4));
morphMatrix=zeros(inNum);
morphMatrix(1:inNum,1:outNum)=pdist2(inInfo(:,3:4),outInfo(:,3:4));
angleMatrix=zeros(inNum);
angleMatrix(1:inNum,1:outNum)=pdist2(inInfo(:,5),outInfo(:,5));
angleMatrix=abs(asin(sin(angleMatrix)))*50;
morphMatrix=(morphMatrix.^2+angleMatrix.^2).^0.5;
if inNum>=9 && inNum==outNum
    for iLine=1:inNum
        order(iLine)=find(distanceMat(iLine,:)==min(distanceMat(iLine,:)));
    end
    orderNew=unique(order);
    if numel(orderNew)==inNum
        afterDivide=order;
        largeDist=0;
        isN2nNode=1;
    else
        afterDivide=[];
        isN2nNode=0;
        largeDist=0;
    end
    return
end
possibleSol=perms(1:inNum);
caculateDist=zeros(size(possibleSol,1),1);
morphDist=zeros(size(possibleSol,1),1);
for i=1:size(possibleSol,1)
    caculateDist(i)=trace(distanceMat(:,possibleSol(i,:)));
    morphDist(i)=trace(morphMatrix(:,possibleSol(i,:)));
end
[minDist,Num]=min(caculateDist);
if inNum==outNum
    [minMorphDist]=min(morphDist);
    if minMorphDist<=5
        [minDist,Num]=min(caculateDist);
        [minNums]=find(morphDist<=5);
        if ismember(Num,minNums)
            threShold=80;
        end
    end
    if minMorphDist<=10 && minMorphDist>5
        [minDist,Num]=min(caculateDist);
        [minNums]=find(morphDist<=10);
        if ismember(Num,minNums)
            threShold=40;
        end
    end
end
if inNum==2 && outNum==2
    if min(min(distanceMat))<=4
        [line,row]=find(distanceMat==min(min(distanceMat)));
        largeDist=0;
        if line==row
            isN2nNode=true;
            afterDivide=possibleSol(2,:);
        else
             isN2nNode=true;
            afterDivide=possibleSol(1,:);
        end
        return
    end
end
% if size(in,2)==2 && (outInfo(1,4)>20 && any(distanceMat(:,2)<20) || outInfo(2,4)>20 && any(distanceMat(:,1)<20))
%     threShold=35;
% end
if minDist>threShold*outNum^0.5
    isN2nNode=false;
    afterDivide=[];
else
    isN2nNode=true;
    afterDivide=possibleSol(Num,:);
end
if minDist>15*outNum^0.5
    largeDist=1;
else
    largeDist=0;
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