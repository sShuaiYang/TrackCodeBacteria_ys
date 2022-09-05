function [afterDivide,canDivideorNot]=ASolution(in,outInfo,pictureSize,finalOneToOneMatchThreShold)
% this function use the idea(find junction points) to divide the BW image
% into small pieces. Then get the orientation of each short line. And
% rotate the model by different angles.
if nargin==3
    finalOneToOneMatchThreShold=25;
end
[regionInfo,outInfo.BWImageNew,isTriadius,openNum]=findJunctionPieces(outInfo.BWImage,'Orientation','Centroid');
if openNum==numel(in)
    if openNum<=9
        [afterDivide,canDivideorNot]=n2nNodeSolution(in,outInfo,pictureSize);
        if canDivideorNot==1
            return
        end
    end
end
[needSeg,afterDivide]=oriInSegmentation(in,outInfo,pictureSize);
if needSeg==0
    canDivideorNot=1;
    return
end
clear afterDivide
aveLen=40;
if numel(regionInfo)~=0
    if isTriadius==0
        threShold=10;
        for i=1:size(in,2)
            in{i}.regionInfo=regionprops(in{i}.BWImage,'Orientation','Centroid','MajorAxisLength','MinorAxisLength');
            in{i}.realCentroid=[in{i}.regionInfo.Centroid(1)+in{i}.xyMin(2),in{i}.regionInfo.Centroid(2)+in{i}.xyMin(1)];
            inInfo(i,:)=in{i}.realCentroid;
            angleIn(i,1)=in{i}.regionInfo.Orientation;
        end
        for i=1:numel(regionInfo)
            regionInfo(i).realCentroid=[regionInfo(i).Centroid(1)+outInfo.xyMin(2),regionInfo(i).Centroid(2)+outInfo.xyMin(1)];
            possibleInfo(i,:)=regionInfo(i).realCentroid;
            anglePos(i,1)=regionInfo(i).Orientation;
        end
        angleDif=pdist2(angleIn,anglePos);
        angleDif(angleDif>90)=180-angleDif(angleDif>90);
        angleDif=angleDif*aveLen/180*pi;
        distMatrix=pdist2(inInfo,possibleInfo);
        distMatrix=(distMatrix.^2+angleDif.^2).^0.5;
        caculateNum=0;
        for i=1:numel(in)
            if min(min(distMatrix))==100
                continue
            end
            caculateNum=caculateNum+1;
            [iIn,iRegion]=find(distMatrix==min(min(distMatrix)));
            [afterDivide{iIn(1)},minDistance(iIn(1)),leftImage{iIn(1)}]=bestMatchingOneNoTri(in{iIn(1)},outInfo,pictureSize,regionInfo(iRegion(1)));
            if ~isempty(afterDivide{iIn(1)})
                distMatrix(iIn(1),:)=100;
                distMatrix(:,iRegion(1))=100;
            end
        end
        if numel(afterDivide)<i
            afterDivide{i}=[];
        end
        if numel(minDistance)<i
            minDistance(i)=0;
        end
        if numel(leftImage)<i
            leftImage{i}=[];
        end
        if numel(in)==numel(regionInfo)+1
            outImage=outInfo.BWImage;
            for i=1:numel(in);
                if ~isempty(leftImage{i})
                    outImage=outImage-leftImage{i};
                end
            end
            outImage=logical(outImage);
            outImage=bwToOneArea(outImage);
            outImage=imerode(outImage,true(3));
            outImage=imdilate(outImage,true(3));
            outImage=bwToOneArea(outImage);
            if numel(outImage(outImage==1))<30 && max(minDistance(minDistance~=0))<=threShold
                canDivideorNot=1;
            else
                canDivideorNot=0;
                isTriadius=1;
            end
        else
            if max(minDistance(minDistance~=0))<=threShold
                canDivideorNot=1;
            else
                canDivideorNot=0;
                isTriadius=1;
            end
        end
    end
    if isTriadius==1
        clear afterDivide angleDif distMatrix leftImage outImage
        for i=1:size(in,2)
            in{i}.regionInfo=regionprops(in{i}.BWImage,'Orientation','Centroid','MajorAxisLength','MinorAxisLength');
            in{i}.realCentroid=[in{i}.regionInfo.Centroid(1)+in{i}.xyMin(2),in{i}.regionInfo.Centroid(2)+in{i}.xyMin(1)];
            angleIn(i,1)=in{i}.regionInfo.Orientation;
        end
        for i=1:numel(regionInfo)
            anglePos(i,1)=regionInfo(i).Orientation;
        end
        angleDif=pdist2(angleIn,anglePos);
        angleDif(angleDif>90)=180-angleDif(angleDif>90);
        angleDif=angleDif*aveLen/180*pi;
        distMatrix=angleDif;
        caculateNum=0;
        for i=1:numel(in)
            if min(min(distMatrix))==100
                continue
            end
            caculateNum=caculateNum+1;
            [iIn,iRegion]=find(distMatrix==min(min(distMatrix)));
            [afterDivide{iIn(1)},canDivideorNotEach(iIn(1)),minDist(iIn(1)),filledImage{iIn(1)}]=bestMatchingOne2Tri(in{iIn(1)},outInfo,pictureSize,regionInfo(iRegion(1)));
            distMatrix(iIn(1),:)=100;
            distMatrix(:,iRegion(1))=100;
        end
        if numel(afterDivide)<i
            afterDivide{i}=[];
            minDist(i)=100;
            canDivideorNotEach(i)=0;
        end
        canDivideorNot=all(canDivideorNotEach);
        if canDivideorNot==0 && min(minDist)<10 && numel(canDivideorNotEach(canDivideorNotEach==0))==1
            rightOne=canDivideorNotEach==1;
            rightOrder=find(canDivideorNotEach==1);
            for i=1:numel(rightOrder)
                if i==1
                    leftImage=outInfo.BWImage-filledImage{rightOrder(1)};
                else
                    leftImage=leftImage-filledImage{rightOrder(i)};
                end
            end
            leftImage=im2bw(leftImage);
            leftImage=imopen(leftImage,ones(4));
            leftImage=bwmorph(leftImage,'bridge');
            leftImage=bwToOneArea(leftImage);
            if max(max(leftImage))==0
                canDivideorNot=0;
            else
                stats=regionprops(leftImage,'Centroid','MajorAxisLength','MinorAxisLength');
                canDivideorNot=0;
                if numel(stats)==1
                    realCentroid=[stats.Centroid(1)+outInfo.xyMin(2),stats.Centroid(2)+outInfo.xyMin(1)];
                    distanceNew=pdist2([realCentroid(1),realCentroid(2),stats.MajorAxisLength,stats.MinorAxisLength],[in{~rightOne}.realCentroid(1),in{~rightOne}.realCentroid(2),in{~rightOne}.regionInfo.MajorAxisLength,in{~rightOne}.regionInfo.MinorAxisLength]);
                    if (distanceNew<=finalOneToOneMatchThreShold && canDivideorNot==0)
                        canDivideorNot=1;
                        afterDivide{~rightOne}=xy2Idx(outInfo.xyMin,leftImage,pictureSize);
                    end
                end
            end
        end
    end
else
    canDivideorNot=0;
end
if canDivideorNot==0;
    for i=1:numel(in)
        afterDivide{i}=xy2Idx(in{i}.xyMin,in{i}.BWImage,pictureSize);
    end
    canDivideorNot=1;
end
end
function [afterDivide,minDistance,outImage]=bestMatchingOneNoTri(in,outInfo,pictureSize,regionInfo)
outImage=false(size(outInfo.BWImage));
regionInfo.position=ceil([regionInfo.Centroid(2),regionInfo.Centroid(1)]);
modelImage.image=imrotate(in.BWImage,regionInfo.Orientation-in.regionInfo.Orientation);
modelImage.image=bwmorph(modelImage.image,'bridge');
modelImage.image=bwToOneArea(modelImage.image);
modelImage.rotateInfo=regionprops(modelImage.image,'Centroid');
modelImage.newImage=false(size(outInfo.BWImage));
if ~(round(regionInfo.position(1)-modelImage.rotateInfo(1).Centroid(2)+(size(modelImage.image,1)+1)/2)<=0 || round(regionInfo.position(1)-modelImage.rotateInfo(1).Centroid(2)+(size(modelImage.image,1)+1)/2)>pictureSize(1) ||...
        round(regionInfo.position(2)-modelImage.rotateInfo(1).Centroid(1)+(size(modelImage.image,2)+1)/2)<=0 || round(regionInfo.position(2)-modelImage.rotateInfo(1).Centroid(1)+(size(modelImage.image,2)+1)/2)>pictureSize(2))
    modelImage.newImage(round(regionInfo.position(1)-modelImage.rotateInfo(1).Centroid(2)+(size(modelImage.image,1)+1)/2),round(regionInfo.position(2)-modelImage.rotateInfo(1).Centroid(1)+(size(modelImage.image,2)+1)/2))=1;
end
modelImage.newImage=imdilate(modelImage.newImage,modelImage.image);
if ~isequal(size(modelImage.newImage),size(outInfo.BWImage))
    afterDivide=[];
    minDistance=100;
    return
end
modelImage.newImage=modelImage.newImage & outInfo.BWImage;
modelImage.newImage=bwmorph(modelImage.newImage,'open');
if max(max(modelImage.newImage))==0
    modelImage.realCentroid=[];
    minDistance=100;
else
    modelImage.newImage=bwToOneArea(modelImage.newImage);
    modelImage.regionInfo=regionprops(modelImage.newImage,'Centroid','MajorAxisLength','MinorAxisLength');
    if numel(modelImage.regionInfo)>=2
        modelImage.realCentroid=[];
        minDistance=100;
    else
        modelImage.realCentroid=[modelImage.regionInfo.Centroid(1)+outInfo.xyMin(2),modelImage.regionInfo.Centroid(2)+outInfo.xyMin(1)];
        minDistance=pdist2([modelImage.realCentroid(1),modelImage.realCentroid(2),modelImage.regionInfo.MajorAxisLength,modelImage.regionInfo.MinorAxisLength],[in.realCentroid(1),in.realCentroid(2),in.regionInfo.MajorAxisLength,in.regionInfo.MinorAxisLength]);
    end
end
if minDistance==100
    afterDivide=[];
    return
end
afterDivide=xy2Idx(outInfo.xyMin,modelImage.newImage,pictureSize);
outImage=modelImage.newImage;
end
function image=bwToOneArea(image)
cc=regionprops(image,'Area');
eachArea=[];
if size(cc,1)>1
    for iArea=1:size(cc,1)
        eachArea(:,iArea)=cc(iArea).Area;
    end
    image=bwareaopen(image,max(eachArea)-1);
end
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
function [afterDivide,canDivideorNot,minDistance,filledImage]=bestMatchingOne2Tri(in,outInfo,pictureSize,regionInfo)
threShold=10;
modelImage.image=imrotate(in.BWImage,regionInfo.Orientation-in.regionInfo.Orientation);
modelImage.image=bwmorph(modelImage.image,'bridge');
modelImage.image=bwToOneArea(modelImage.image);
modelImage.image1=bwmorph(modelImage.image,'spur');
modelImage.image2=bwmorph(modelImage.image,'erode');
modelImage.newImage=imerode(outInfo.BWImage,modelImage.image1);
if max(max(modelImage.newImage))==0
    modelImage.newImage=imerode(outInfo.BWImage,modelImage.image2);
    if max(max(modelImage.newImage))==0
        minDistance=100;
    else
        [minDistance,modelImage]=erodeDilateModel(modelImage,outInfo,in);
    end
else
    [minDistance,modelImage]=erodeDilateModel(modelImage,outInfo,in);
end
if minDistance>=10
    clear minDistance;
    modelImage.image3=imdilate(outInfo.BWImage,ones(4));
    modelImage.image3(1,:)=false;
    modelImage.image3(end,:)=false;
    modelImage.image3(:,1)=false;
    modelImage.image3(:,end)=false;
    modelImage.newImage=imerode(modelImage.image3,modelImage.image2);
    modelImage.newImage=bwmorph(modelImage.newImage,'shrink');
    if max(max(modelImage.newImage))==0
        minDistance=100;
    else
        [minDistance,modelImage]=erodeDilateModel(modelImage,outInfo,in);
    end
    if min(minDistance)==100
        afterDivide=[];
        canDivideorNot=0;
        minDistance=100;
        filledImage=[];
        return
    end
end
if min(minDistance)>threShold
    canDivideorNot=0;
    afterDivide=[];
    filledImage=[];
else
    canDivideorNot=1;
    afterDivide=xy2Idx(outInfo.xyMin,modelImage.newImage,pictureSize);
    filledImage=modelImage.newImage;
end
end
function [minDistance,modelImage]=erodeDilateModel(modelImage,outInfo,in)
modelImage.newImage=bwmorph(modelImage.newImage,'thin',inf);
modelImage.newImage=bwmorph(modelImage.newImage,'endpoints');
modelImage.stats=regionprops(modelImage.newImage,'PixelList');
centroidDis=[];
for iSmallPiece=1:numel(modelImage.stats)
    centroidDis(iSmallPiece)=pdist2([modelImage.stats(iSmallPiece).PixelList(1,1)+outInfo.xyMin(2),modelImage.stats(iSmallPiece).PixelList(1,2)+outInfo.xyMin(1)],[in.realCentroid(1),in.realCentroid(2)]);
end
modelImage.newImage=false(size(outInfo.BWImage));
smallProperNum=find(centroidDis==min(centroidDis));
if isempty(smallProperNum)
    minDistance=100;
    modelImage=[];
    return
end
smallProperNum=smallProperNum(1);
modelImage.newImage(modelImage.stats(smallProperNum).PixelList(1,2),modelImage.stats(smallProperNum).PixelList(1,1))=1;
modelImage.newImage=imdilate(modelImage.newImage,modelImage.image);
modelImage.newImage=imdilate(modelImage.newImage,ones(2));
modelImage.newImage=modelImage.newImage & outInfo.BWImage;
modelImage.newImage=bwmorph(modelImage.newImage,'close');
modelImage.newImage=bwmorph(modelImage.newImage,'bridge');
minDistance=100;
if max(max(modelImage.newImage))==1
    modelImage.newImage=bwToOneArea(modelImage.newImage);
    modelImage.regionInfo=regionprops(modelImage.newImage,'Centroid','MajorAxisLength','MinorAxisLength');
    if numel(modelImage.regionInfo)==1
        modelImage.realCentroid=[modelImage.regionInfo.Centroid(1)+outInfo.xyMin(2),modelImage.regionInfo.Centroid(2)+outInfo.xyMin(1)];
        minDistance=pdist2([modelImage.realCentroid(1),modelImage.realCentroid(2),modelImage.regionInfo.MajorAxisLength,modelImage.regionInfo.MinorAxisLength],[in.realCentroid(1),in.realCentroid(2),in.regionInfo.MajorAxisLength,in.regionInfo.MinorAxisLength]);
    end
end
end
function [afterDivide,isN2nNode]=n2nNodeSolution(in,outOri,pictureSize)
% used to deal with N-to-M node and divide them,here N>=M
threShold=15;
outOri.BWImage=outOri.BWImageNew;
pixelIdxList=xy2Idx(outOri.xyMin,outOri.BWImage,pictureSize);
[regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize);
for i=1:regionNum
    out{i}.xyMin=xyMin;
    out{i}.BWImage=regionImage{i};
end
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
        Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
        xyMin=in{i}.xyMin;
        inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,:)=zeros(1,5);
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
if minDist>threShold*numel(in)
    isN2nNode=false;
    for i=1:outNum
        afterDivide{i}=[];
    end
else
    isN2nNode=true;
    rightOrder=possibleSol(Num,:);
    for i=1:outNum
        afterDivide{i}=xy2Idx(out{rightOrder(i)}.xyMin,out{rightOrder(i)}.BWImage,pictureSize);
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
function [needSeg,afterDivide]=oriInSegmentation(in,out,pictureSize)
threShold=10;
for i=1:size(in,2)
    inNew{i}.BWImage=gainSameSize(in{i},out,pictureSize);
    inNew{i}.xyMin=out.xyMin;
    inNew{i}.BWImage=inNew{i}.BWImage & out.BWImage;
    inNew{i}.BWImage=bwmorph(inNew{i}.BWImage,'bridge');
    inNew{i}.BWImage=bwToOneArea(inNew{i}.BWImage);
    inNew{i}.regionInfo=regionprops(inNew{i}.BWImage,'Orientation','Centroid','MajorAxisLength','MinorAxisLength');
    if numel(inNew{i}.regionInfo)>=2 || max(max(inNew{i}.BWImage))==0
        minDistance(i)=100;
        continue
    end
    inNew{i}.realCentroid=[inNew{i}.regionInfo.Centroid(1)+inNew{i}.xyMin(2),inNew{i}.regionInfo.Centroid(2)+inNew{i}.xyMin(1)];
    in{i}.regionInfo=regionprops(in{i}.BWImage,'Orientation','Centroid','MajorAxisLength','MinorAxisLength');
    in{i}.realCentroid=[in{i}.regionInfo.Centroid(1)+in{i}.xyMin(2),in{i}.regionInfo.Centroid(2)+in{i}.xyMin(1)];
    minDistance(i)=pdist2([inNew{i}.realCentroid(1),inNew{i}.realCentroid(2),inNew{i}.regionInfo.MajorAxisLength,inNew{i}.regionInfo.MinorAxisLength],[in{i}.realCentroid(1),in{i}.realCentroid(2),in{i}.regionInfo.MajorAxisLength,in{i}.regionInfo.MinorAxisLength]);
end
if max(minDistance)<=threShold
    needSeg=0;
    for i=1:size(in,2)
        afterDivide{i}=xy2Idx(in{i}.xyMin,in{i}.BWImage,pictureSize);
    end
else
    needSeg=1;
    afterDivide=[];
end
end
function newA=gainSameSize(A,B,pictureSize)
% change picture in A to the same size that in B
pixelIdxList=xy2Idx(A.xyMin,A.BWImage,pictureSize);
new=false(pictureSize);
new(pixelIdxList)=true;
bSize=size(B.BWImage);
newA=false(bSize);
newA(2:end-1,2:end-1)=new(B.xyMin(1):B.xyMin(1)+bSize(1)-3,B.xyMin(2):B.xyMin(2)+bSize(2)-3);
end