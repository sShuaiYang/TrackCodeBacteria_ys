function [afterDivide,canDivideorNot]=basicDivideSolution(eachInfo,maskImage,pictureSize,strength,finalOneToOneMatchThreShold)
if strength==1
    areaDiffThr=0.3;
end
if strength==2
    areaDiffThr=0.6;
end
if strength==3
    areaDiffThr=+inf;
end
inputNum=numel(eachInfo);
realIn=0;
inImage=[];
for i=1:inputNum
    realNum(i)=0;
    [eachNum,xyMin,regionImage]=findRegionNum(eachInfo{i},pictureSize);
    if eachNum>=2
        canDivideorNot=0;
        afterDivide=[];
        return
    end
    in{i}.xyMin=xyMin;
    in{i}.BWImage=regionImage{1};
    if ~isempty(in{i}.BWImage)
        realIn=realIn+1;
        realNum(i)=1;
        stats=regionprops(in{i}.BWImage,'Area','MajorAxisLength','MinorAxisLength');
        if abs(stats.MajorAxisLength*stats.MinorAxisLength*pi/4-stats.Area)/stats.Area>=areaDiffThr
            canDivideorNot=0;
            afterDivide=[];
            return
        end
    end
    inImage=[inImage;eachInfo{i}];
end
if all(realNum==0)
    canDivideorNot=0;
    afterDivide=[];
    return
end
[outputNum,xyMin,regionImage]=findRegionNum(maskImage,pictureSize);
for i=1:outputNum
    out{i}.xyMin=xyMin;
    out{i}.BWImage=regionImage{i};
end
% close all
% for i=1:inputNum
%     figure;imshow(in{i}.BWImage)
% end
% for i=1:outputNum
%     figure;imshow(out{i}.BWImage)
% end
if realIn<outputNum
    canDivideorNot=false;
    afterDivide=[];
    return
end
if realIn==outputNum && outputNum<=9
    [afterDivide,canDivideorNot]=n2mNodeSolution(in,out,pictureSize);
    return
end
if realIn>outputNum || realIn==outputNum && outputNum>=10
    [inArea,minArea]=getSumArea(inImage,in,inputNum,pictureSize);
    [outArea,~]=getSumArea(maskImage,[],outputNum,pictureSize);
    if (inArea-outArea)/(realIn-outputNum)/minArea>=0.8
        [afterDivide,canDivideorNot]=n2mNodeSolution(in,out,pictureSize);
    else
        for i=1:size(out,2)
            out{1}.BWImage=out{1}.BWImage|out{i}.BWImage;
        end
        [afterDivide(realNum==1),canDivideorNot]=n21NodeDivide(in(realNum==1),out,pictureSize,finalOneToOneMatchThreShold);
        if numel(afterDivide)<numel(in)
            afterDivide{numel(in)}=[];
        end
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
    end
end
% figure;imshow(out{1}.BWImage)
% p=[];
% for i=1:numel(afterDivide)
%    p=[p;afterDivide{i}];
% end
% [~,image]=idx2Xy(p,pictureSize);
% figure;imshow(image)
end
%% this  four function is used to divide N to 1 Node,by rotating the model and erode the out image
function [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize,finalOneToOneMatchThreShold)
outInfo=out{1};
[afterDivide,canDivideorNot]=ASolution(in,outInfo,pictureSize,finalOneToOneMatchThreShold);
end
%% this function is used to divide N to M Node directly,choose the nearest bacteria
function [afterDivide,isN2mNode]=n2mNodeSolution(in,out,pictureSize)
% used to deal with N-to-M node and divide them,here N>=M
threShold=26;
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
        Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
        xyMin=in{i}.xyMin;
        inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,:)=zeros(1,4);
    end
end
outNum=numel(out);
for i=1:numel(out)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
end
distanceMat=zeros(inNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
possibleSol=perms(1:inNum);
caculateDist=zeros(size(possibleSol,1),1);
for i=1:size(possibleSol,1)
    caculateDist(i)=trace(distanceMat(:,possibleSol(i,:)));
end
[minDist,Num]=min(caculateDist);
if minDist>threShold*outNum
    isN2mNode=false;
    afterDivide=[];
else
    isN2mNode=true;
    correctVec=possibleSol(Num,:);
    for i=1:inNum
        if correctVec(i)>outNum
            afterDivide{i}=[];
        else
            afterDivide{i}=xy2Idx(out{correctVec(i)}.xyMin,out{correctVec(i)}.BWImage,pictureSize);
        end
    end
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
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
% BWImage=bwmorph(BWImage,'remove');
pixelIdxListOri=find(BWImage==1);
if size(pixelIdxListOri,2)~=1
    pixelIdxListOri=pixelIdxListOri';
end
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function [sumArea,minArea]=getSumArea(idxInfo,info,num,pictureSize)
[~,BWImage]=idx2Xy(idxInfo,pictureSize);
sumArea=numel(BWImage(BWImage==1));
minArea=[];
if ~isempty(info)
    areaInfo=zeros(num,1);
    for i=1:num
        filledImage=info{i}.BWImage;
        areaInfo(i)=numel(filledImage(filledImage==1));
    end
    minArea=min(areaInfo(areaInfo~=0));
end
end