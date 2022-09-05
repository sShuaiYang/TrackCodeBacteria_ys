function afterDivision=morphologySegmentation_new(maskImage,eachInfo,pictureSize)
sigma=0.9;
sita=-30:2:30;
alpha=0.5;
[OriginalXY,beforeDivision]=Idx2xy(maskImage,pictureSize);
n=0;
for i=1:numel(eachInfo)
    if ~isempty(eachInfo{i})
        n=n+1;
    end
end
if n==1
    AllInfo=regionprops(beforeDivision,'PixelIdx');
    for i=1:numel(AllInfo)
        afterDivision{i}=AllInfo(i).PixelIdxList;
    end
    if numel(eachInfo)>numel(AllInfo)
        for i=numel(AllInfo)+1:numel(eachInfo)
            afterDivision{i}=[];
        end
    end
else
    rudeCase=bwconncomp(beforeDivision);
    if rudeCase.NumObjects==numel(eachInfo)
        AllInfo=regionprops(rudeCase,'Centroid','MajorAxisLength','MinorAxisLength','PixelIdx');
        for i=1:numel(AllInfo)
            beforeDivisionInfo(i,:)=[AllInfo(i,1).Centroid(1),AllInfo(i,1).Centroid(2),AllInfo(i,1).MajorAxisLength,AllInfo(i,1).MajorAxisLength];
        end
        for i=1:numel(eachInfo)
            beforeImageCC=makeCC(eachInfo{i},pictureSize);
            Info=regionprops(beforeImageCC,'Centroid','MajorAxisLength','MinorAxisLength','PixelIdx');
            beforeInfo=[Info.Centroid(1),Info.Centroid(2),Info.MajorAxisLength,Info.MinorAxisLength];
            caculateDistance=sqrt(sum((AminusB(beforeDivisionInfo,beforeInfo).^2),2));
            [~,Num]=min(caculateDistance);
            afterDivision{i}=AllInfo(Num,1).PixelIdxList;
            %                 [~,beforeImage]=revertIdx2xy(eachInfo{i},pictureSize);
            %                 figure;imshow(beforeImage);
            %                 [~,u]=revertIdx2xy(afterDivision{i},pictureSize);
            %                 figure;imshow(u)
        end
        %             figure;imshow(beforeDivision)
    else
        beforeDivision=imfill(beforeDivision,'holes');
        [isonemiss,canDivideorNot]=isOneMissing(maskImage,eachInfo,pictureSize,alpha,sigma);
        if canDivideorNot==0
            afterDivision{1}=maskImage;
            if numel(eachInfo)>1
                for i=2:numel(eachInfo)
                    afterDivision{i}=[];
                end
            end
        else
            for i=1:numel(eachInfo)
                [~,afterDivision{i},centroidInfo(i,1)]=getBacteria(pictureSize,eachInfo{i},OriginalXY,beforeDivision,sita);
            end
            if isonemiss==1
                [~,num]=find(centroidInfo==max(centroidInfo));
                afterDivision{num}=[];
            end
        end
        %         figure;imshow(beforeDivision1)
    end
end
end
function [xyMin,BWImage]=Idx2xy(pixelIdxList,pictureSize)
xSize=pictureSize(1);
yReal=fix(pixelIdxList./xSize)+1;
xReal=pixelIdxList-(yReal-1).*xSize;
xMin=min(xReal);
xMax=max(xReal);
yMin=min(yReal);
yMax=max(yReal);
xReal=xReal-xMin+1;
yReal=yReal-yMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList=xReal+(yReal-1).*(xMax-xMin+1);
BWImage(pixelIdxList)=1;
xyMin=[xMin,yMin];
end
function result=AminusB(A,B)
for i=1:size(A,1)
    result(i,:)=A(i,:)-B;
end
end
function [beforeDivision1,afterDivision,centroidInfo]=getBacteria(pictureSize,eachInfo,OriginalXY,beforeDivision,sita)
if isempty(eachInfo)
    afterDivision=[];
    beforeDivision1=[];
    centroidInfo=10000;
else
[xyResult,beforeImage]=revertIdx2xy(eachInfo,pictureSize);
xSize=[min(min(OriginalXY(:,1)),min(xyResult(:,1))),max(max(OriginalXY(:,1)),max(xyResult(:,1)))];
ySize=[min(min(OriginalXY(:,2)),min(xyResult(:,2))),max(max(OriginalXY(:,2)),max(xyResult(:,2)))];
beforeDivision1=beforeDivision(xSize(1):xSize(2),ySize(1):ySize(2));
beforeImageCac=beforeImage(xSize(1):xSize(2),ySize(1):ySize(2));
beforeImageCac=imfill(beforeImageCac,'holes');
ccbefore=regionprops(beforeImageCac,'Centroid','MajorAxisLength','MinorAxisLength','Orientation');
ccbeforeInfo=[ccbefore.Centroid(1),ccbefore.Centroid(2),ccbefore.MajorAxisLength,ccbefore.MinorAxisLength];
eachXSize=[min(xyResult(:,1)),max(xyResult(:,1))];
eachYSize=[min(xyResult(:,2)),max(xyResult(:,2))];
beforeImage=beforeImageCac(eachXSize(1)-xSize(1)+1:eachXSize(2)-xSize(1)+1,eachYSize(1)-ySize(1)+1:eachYSize(2)-ySize(1)+1);
t=0;
ellipseInfo=[];
eachBWInfo=[];
orientationInfo=[];
for j=1:numel(sita)
    BW=imrotate(beforeImage,sita(j));
    BW2=bwmorph(BW,'erode');
    I=imerode(beforeDivision1,BW2);
    if max(max(I))~=0
        cc=bwconncomp(I);
        for p=1:cc.NumObjects
            erodePaper=false(xSize(2)-xSize(1)+1,ySize(2)-ySize(1)+1);
            erodePaper(cc.PixelIdxList{p})=1;
            erodePaper=imdilate(erodePaper,BW);
            erodePaper=bwmorph(erodePaper,'thin');
            erodePaper=bwmorph(erodePaper,'close');
            %                 erodePaper=bwmorph(erodePaper,'dilate');
            STATS=regionprops(erodePaper,'Centroid','MajorAxisLength','MinorAxisLength','Orientation','PixelIdxList');
            if numel(STATS)==1
                t=t+1;
                ellipseInfo(t,:)=[STATS.Centroid(1),STATS.Centroid(2),STATS.MajorAxisLength,STATS.MinorAxisLength];
                orientationInfo(t,:)=abs(STATS.Orientation-ccbefore.Orientation);
                orientationInfo(orientationInfo>=90)=orientationInfo(orientationInfo>=90)-90;
                eachBWInfo(t).pixelList=STATS.PixelIdxList;
            end
        end
    end
end
if t==0
    afterDivision=[];
    centroidInfo=10000;
else
caculateDistance=sqrt(sum((AminusB(ellipseInfo,ccbeforeInfo).^2),2)+(orientationInfo/180*pi*ccbefore.MajorAxisLength).^2);
[~,modelNum]=min(caculateDistance);
centroidInfo=sqrt((sum(AminusB(ellipseInfo(modelNum,1:2),ccbeforeInfo(1,1:2))).^2));
erodePaper=false(pictureSize(1),pictureSize(2));
erodePaper1=false(xSize(2)-xSize(1)+1,ySize(2)-ySize(1)+1);
erodePaper1(eachBWInfo(modelNum).pixelList)=1;
properBacteria=erodePaper1;
properBacteria=bwmorph(properBacteria,'remove');
%     properBacteria=bwmorph(properBacteria,'dilate');
erodePaper(xSize(1):xSize(2),ySize(1):ySize(2))=properBacteria;
properBacteriaInfo=regionprops(erodePaper,'PixelIdxList');
afterDivision=properBacteriaInfo.PixelIdxList;
% figure;imshow(beforeImageCac)
% figure;imshow(properBacteria)
end
end
end
function [isonemiss,canDivideorNot]=isOneMissing(maskImage,eachInfo,pictureSize,alpha,sigma)
for i=1:numel(eachInfo)
    if ~isempty(eachInfo{i})
    [~,beforeImage]=revertIdx2xy(eachInfo{i},pictureSize);
    Info=regionprops(beforeImage,'FilledArea');
    Area(i,1)=Info.FilledArea;
    end
end
AllAreaBefore=sum(Area);
AreaMin=min(Area);
[~,beforeDivision]=revertIdx2xy(maskImage,pictureSize);
Info=regionprops(beforeDivision,'FilledArea');
maskImageArea=0;
for i=1:numel(Info)
    maskImageArea=maskImageArea+Info(i).FilledArea;
end
if abs(double((maskImageArea-AllAreaBefore)/maskImageArea))>alpha
    canDivideorNot=0;
else
    canDivideorNot=1;
end
if double((AllAreaBefore-maskImageArea)/AreaMin)>sigma
    isonemiss=1;
else
    isonemiss=0;
end
end
function ImageCC=makeCC(eachInfo,pictureSize)
ImageCC.Connectivity=8;
ImageCC.ImageSize=pictureSize;
ImageCC.NumObjects=1;
ImageCC.PixelIdxList{1,1}=eachInfo;
end