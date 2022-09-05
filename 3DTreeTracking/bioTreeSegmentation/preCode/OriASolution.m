function [afterDivide,canDivideorNot]=OriASolution(in,outInfo,pictureSize,outEndInfo)
outImage=outInfo.BWImage;
outImage=imfill(outImage,'holes');
xyMin=outInfo.xyMin;
parfor i=1:numel(in)
    emptyDivide(i)=0;
    if ~isempty(in{i}.BWImage)
        sita=-30:3:30;
        imageStack=rotateStacks(in{i}.BWImage,sita);
        newInImage=gainSameSize(in{i},outInfo,pictureSize);
        [imageStack,skelImage]=gainPossibleModel(imageStack,outImage,newInImage);
        inxyMin=in{i}.xyMin;
        remainImage=outImage-newInImage;
        afterDividePre{i}=findTheBest(in{i}.BWImage,inxyMin,imageStack,xyMin,sita,inf,remainImage,skelImage);
        
        sita=afterDividePre{i}.sitaNum-3:1:afterDividePre{i}.sitaNum+3;
        imageStack=rotateStacks(in{i}.BWImage,sita);
        [imageStack,skelImage]=gainPossibleModel(imageStack,outImage,newInImage);
        afterDividePre{i}=findTheBest(in{i}.BWImage,inxyMin,imageStack,xyMin,sita,30,remainImage,skelImage);
    else
        afterDividePre{i}=[];
    end
    if ~isempty(afterDividePre{i})
        emptyDivide(i)=1;
    end
end
if all(emptyDivide==0)
    afterDivide=[];
    canDivideorNot=0;
    return
end
% isN2n=isN2nMatching(in,outImage);
% if isN2n==1 && any(emptyDivide==0)
%     afterDivide=[];
%     canDivideorNot=0;
%     return
% end
[afterDivide,canDivideorNot,properNum,bacteriaImage,afterDividePre]=getProperbacteria(afterDividePre,xyMin,pictureSize,in,outImage,outEndInfo);
if size(in,2)==2
    if (any(properNum==2) || any(properNum==1)) && any(properNum==3) ||...
            any(properNum==1) && any(properNum==4)
     % here num in properNum means: 1 means minDist<=2,
     % 2---minDist<=threShold ,3--minDist>threShold ,4--no proper erode model
     % 5--wrong segementation
    [afterDivide,canDivideorNot]=tryToFindBacteria(afterDivide,bacteriaImage,afterDividePre,canDivideorNot,properNum,in,outImage,xyMin,pictureSize);
    end
end
isN2n=isN2nMatching(in,outImage);
if isN2n==1 && any(cellfun(@isempty,afterDivide)==1)
    afterDivide=[];
    canDivideorNot=0;
    return
end
end
function [afterDivide,canDivideorNot]=tryToFindBacteria(afterDivide,bacteriaImage,afterDividePre,canDivideorNot,properNum,in,outImage,xyMin,pictureSize)
rightNum=find(properNum==1);
if isempty(rightNum)
    rightNum=find(properNum==2);
end
if isempty(in{3-rightNum}.BWImage)
    return
end
leftImage=outImage-bacteriaImage{rightNum};
leftImage=bwmorph(leftImage,'open');
leftImage=bwToOneArea(leftImage);
if max(max(leftImage))==0
    return
end
if numel(leftImage(leftImage==1))<=120
    return
end
inxyMin=in{3-rightNum}.xyMin;
inImage=imfill(in{3-rightNum}.BWImage,'holes');
inImage=bwToOneArea(inImage);
inInfo=regionprops(inImage,'Centroid','PixelIdx','MajorAxisLength','MinorAxisLength');
inImageInfo=[inInfo.Centroid(2)+inxyMin(1),inInfo.Centroid(1)+inxyMin(2),inInfo.MajorAxisLength,inInfo.MinorAxisLength];
outInfo=regionprops(leftImage,'Centroid','PixelIdx','FilledArea','MajorAxisLength','MinorAxisLength');
if numel(outInfo)>1
    return
end
outImageInfo=[outInfo.Centroid(2)+xyMin(1),outInfo.Centroid(1)+xyMin(2),outInfo.MajorAxisLength,outInfo.MinorAxisLength];
dist=min(pdist2(inImageInfo(1:2),outImageInfo(1:2)),pdist2(inImageInfo(3:4),outImageInfo(3:4)));
if dist>25
    return
end
if ~isempty(afterDivide{3-rightNum})
    if dist>=afterDividePre{3-rightNum}.minDist
        return
    end
end

canDivideorNot=1;
leftImage=bwmorph(leftImage,'remove');
afterDivide{3-rightNum}=xy2Idx(xyMin,leftImage,pictureSize);
end

function isN2n=isN2nMatching(in,outImage)
areaIn=0;
for i=1:numel(in)
    if ~isempty(in{i}.BWImage)
    cc=regionprops(in{i}.BWImage,'FilledArea');
    areaIn=areaIn+cc.FilledArea;
    end
end
image=imfill(outImage,'holes');
areaOut=numel(image(image==1));
if abs(areaIn-areaOut)<=50
    isN2n=1;
else
    isN2n=0;
end
end
function imageStack=rotateStacks(BWImage,sita)
% fill in the image and rotate them to stacks
BWImage=imfill(BWImage,'holes');
parfor j=1:numel(sita)
    imageStack{j}=imrotate(BWImage,sita(j));
end
end
function [imageStack,beforeOverlap]=gainPossibleModel2(imageStack,outImage,inImage)
% get all possible model
parfor j=1:numel(imageStack)
    neednext=1;
    testI=imerode(outImage,imageStack{j});
    cc=bwconncomp(testI);
    if cc.NumObjects==1
        image=imdilate(testI,imageStack{j});
        image=bwareaopen(image,10);
        overlapArea=image&inImage;
        if double(numel(overlapArea(overlapArea==1))/numel(inImage(inImage==1)))>0.5
            beforeOverlap{j}=image;
            image=image&outImage;
            if max(max(image))~=0
                imageStack{j}=image;
                neednext=0;
            end
        end
    else if cc.NumObjects>=2
            for i=1:cc.NumObjects
                eachone=zeros(cc.ImageSize);
                eachone(cc.PixelIdxList{i})=1;
                imageGain=imdilate(eachone,imageStack{j});
                overlapArea=imageGain&inImage;
                if double(numel(overlapArea(overlapArea==1))/numel(inImage(inImage==1)))>0.5
                    imageGain=bwareaopen(imageGain,10);
                    beforeOverlap{j}=imageGain;
                    imageGain=imageGain&outImage;
                    if max(max(imageGain))~=0
                        imageStack{j}=imageGain;
                        neednext=0;
                        break
                    end
                end
            end
        end
    end
    if neednext==1
        model=bwmorph(imageStack{j},'thin');
        I=imerode(outImage,model);
        I=bwmorph(I,'shrink');
        I=bwmorph(I,'thin',1);
        I=imclearborder(I);
        cc=bwconncomp(I);
        if cc.NumObjects==0
            model2=bwmorph(imageStack{j},'shrink');
            I=imerode(outImage,model2);
            I=bwmorph(I,'shrink');
            I=bwmorph(I,'thin',2);
            I=imclearborder(I);
            cc=bwconncomp(I);
            if cc.NumObjects==1
                imageStack{j}=imdilate(I,imageStack{j});
                imageStack{j}=bwareaopen(imageStack{j},10);
                beforeOverlap{j}=imageStack{j};
                imageStack{j}=imageStack{j}&outImage;
                if max(max(imageStack{j}))==0
                    imageStack{j}=[];
                end
            else if cc.NumObjects==0
                    imageStack{j}=[];
                else if cc.NumObjects>=2
                        result=0;
                        for i=1:cc.NumObjects
                            eachone=false(cc.ImageSize);
                            eachone(cc.PixelIdxList{i})=1;
                            imageGain=imdilate(eachone,imageStack{j});
                            imageGain=imageGain&inImage;
                            if double(numel(imageGain(imageGain==1))/numel(inImage(inImage==1)))>0.3
                                result=1;
                                imageStack{j}=imdilate(eachone,imageStack{j});
                                imageStack{j}=bwareaopen(imageStack{j},10);
                                beforeOverlap{j}=imageStack{j};
                                imageStack{j}=imageStack{j}&outImage;
                                if max(max(imageStack{j}))==0
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
        else if cc.NumObjects==1
                imageStack{j}=imdilate(I,imageStack{j});
                imageStack{j}=bwareaopen(imageStack{j},10);
                beforeOverlap{j}=imageStack{j};
                imageStack{j}=imageStack{j}&outImage;
                if max(max(imageStack{j}))==0
                    imageStack{j}=[];
                end
            else if cc.NumObjects>=2
                    result=0;
                    for i=1:cc.NumObjects
                        eachone=false(cc.ImageSize);
                        eachone(cc.PixelIdxList{i})=1;
                        imageGain=imdilate(eachone,imageStack{j});
                        imageGain=imageGain&inImage;
                        if double(numel(imageGain(imageGain==1))/numel(inImage(inImage==1)))>0.3
                            result=1;
                            imageStack{j}=imdilate(eachone,imageStack{j});
                            imageStack{j}=bwareaopen(imageStack{j},10);
                            beforeOverlap{j}=imageStack{j};
                            imageStack{j}=imageStack{j}&outImage;
                            if max(max(imageStack{j}))==0
                                imageStack{j}=[];
                            end
                            break
                        end
                    end
                    if result==0
                        imageStack{j}=[];
                    end
                else
                    imageStack{j}=[];
                end
            end
        end
    end
end
end
function [imageStack,skelImage]=gainPossibleModel(imageStack,outImage,inImage)
% get all possible model
parfor j=1:numel(imageStack)
    neednext=1;
    testI=imerode(outImage,imageStack{j});
    testI=testI & outImage;
    cc=bwconncomp(testI);
    if cc.NumObjects==1
        skelImage{j}=testI;
        image=imdilate(testI,imageStack{j});
        image=bwareaopen(image,10);
        overlapArea=image&inImage;
        if double(numel(overlapArea(overlapArea==1))/numel(inImage(inImage==1)))>0.7
            image=image&outImage;
            image=bwmorph(image,'open');
            image=bwToOneArea(image);
            if max(max(image))~=0
                imageStack{j}=image;
                neednext=0;
            end
        end
    else if cc.NumObjects>=2
            for i=1:cc.NumObjects
                eachone=zeros(cc.ImageSize);
                eachone(cc.PixelIdxList{i})=1;
                imageGain=imdilate(eachone,imageStack{j});
                overlapArea=imageGain&inImage;
                if double(numel(overlapArea(overlapArea==1))/numel(inImage(inImage==1)))>0.7
                    skelImage{j}=eachone;
                    imageGain=bwareaopen(imageGain,10);
                    imageGain=imageGain&outImage;
                    imageGain=bwmorph(imageGain,'open');
                    imageGain=bwToOneArea(imageGain);
                    if max(max(imageGain))~=0
                        imageStack{j}=imageGain;
                        neednext=0;
                        break
                    end
                end
            end
        end
    end
    if neednext==1
        outImage1=bwmorph(outImage,'dilate');
        outImage2=false(size(outImage1)+2);
        outImage2(2:end-1,2:end-1)=outImage1;
        I=imerode(outImage2,imageStack{j});
        I=bwmorph(I,'thin');
        I=bwmorph(I,'spur');
        I=I(2:end-1,2:end-1);
        I=imclearborder(I);
        I=I & outImage;
        cc=bwconncomp(I);
        if cc.NumObjects==0
            model2=bwmorph(imageStack{j},'shrink');
            model2=bwmorph(model2,'spur');
            I=imerode(outImage2,model2);
            I=bwmorph(I,'shrink');
            I=bwmorph(I,'spur');
            I=I(2:end-1,2:end-1);
            I=imclearborder(I);
            I=I & outImage;
            cc=bwconncomp(I);
            if cc.NumObjects==0
                imageStack{j}=[];
                skelImage{j}=[];
            else
                [imageStack{j},skelImage{j}]=getOutImageFromModel(I,imageStack{j},cc,outImage,inImage);
            end
        else
            [imageStack{j},skelImage{j}]=getOutImageFromModel(I,imageStack{j},cc,outImage,inImage);
        end
    end
end
end
function [imageStack,skelImage]=getOutImageFromModel(erodeImage,oriModel,cc,outImage,inImage)
if cc.NumObjects==1
    skelImage=erodeImage;
    imageStack=imdilate(erodeImage,oriModel);
    imageStack=bwareaopen(imageStack,10);
    imageStack=imageStack&outImage;
    imageStack=bwmorph(imageStack,'open');
    imageStack=bwToOneArea(imageStack);
    if max(max(imageStack))==0
        imageStack=[];
        skelImage=[];
    end
else if cc.NumObjects>=2
        for i=1:cc.NumObjects
            eachone=false(cc.ImageSize);
            eachone(cc.PixelIdxList{i})=1;
            imageGain=imdilate(eachone,oriModel);
            imageGain=imageGain&inImage;
            ratio(i)=double(numel(imageGain(imageGain==1))/numel(inImage(inImage==1)));
        end
        if max(ratio)>0.3
            numOb=find(ratio==max(ratio));
            numOb=numOb(fix((numel(numOb)/2)+1));
            eachone=false(cc.ImageSize);
            eachone(cc.PixelIdxList{numOb})=1;
            skelImage=eachone;
            imageStack=imdilate(eachone,oriModel);
            imageStack=bwareaopen(imageStack,10);
            imageStack=imageStack&outImage;
            imageStack=bwmorph(imageStack,'open');
            imageStack=bwToOneArea(imageStack);
            if max(max(imageStack))==0
                imageStack=[];
                skelImage=[];
            end
        else
            imageStack=[];
            skelImage=[];
        end
    end
end
end
function result=findTheBest(inImage,inxyMin,imageStack,xyMin,sita,threShold,remainImage,skelImage)
% the best choice to the input bacteria
inImage=imfill(inImage,'holes');
inInfo=regionprops(inImage,'Centroid','PixelIdx','MajorAxisLength','MinorAxisLength');
inImageInfo=[inInfo.Centroid(2)+inxyMin(1),inInfo.Centroid(1)+inxyMin(2),inInfo.MajorAxisLength,inInfo.MinorAxisLength];
outImageInfo=zeros(numel(imageStack),4);
outAreaInfo=zeros(numel(imageStack),1);
usefulData=zeros(numel(imageStack),1);
parfor i=1:numel(imageStack)
    if isempty(imageStack{i})
        usefulData(i)=0;
        continue
    else
        usefulData(i)=1;
        outInfo=regionprops(imageStack{i},'Centroid','PixelIdx','FilledArea','MajorAxisLength','MinorAxisLength');
        outImageInfo(i,:)=[outInfo(1,1).Centroid(2)+xyMin(1),outInfo(1,1).Centroid(1)+xyMin(2),outInfo(1,1).MajorAxisLength,outInfo(1,1).MinorAxisLength];
        outAreaInfo(i,:)=outInfo(1,1).FilledArea;
    end
end
caculateDist=pdist2(outImageInfo,inImageInfo);
caculateDist(usefulData==0)=max(caculateDist);
if min(caculateDist)<=5
    minCacDist=5;
else
    minCacDist=min(caculateDist)+1;
end
outAreaInfo(caculateDist>minCacDist)=0;
maxArea=max(outAreaInfo(caculateDist<=minCacDist));

% NumMaxArea=find(outAreaInfo>=maxArea-2);
% for i=1:numel(NumMaxArea)
%     if ~isempty(imageStack{NumMaxArea(i)})
%     overlapImage=imageStack{NumMaxArea(i)} & remainImage;
%     overlapArea(i)=numel(overlapImage(overlapImage==1));
%     else
%         overlapArea(i)=200;
%     end
% end
% Num=find(overlapArea==min(overlapArea));

Num=find(outAreaInfo==maxArea);
Num=Num(ceil((numel(Num)+1)/2));
% Num=NumMaxArea(Num);
minDist=caculateDist(Num);
if minDist>threShold
    result=[];
else
    correctImage=imageStack{Num};
    result.Image=correctImage;
    result.sitaNum=sita(Num);
    result.skelImage=skelImage{Num};
    result.inInfo=inImageInfo;
end
end
function [afterDivide,canDivideorNot,properNum,bacteriaImage,afterDividePre]=getProperbacteria(afterDividePre,xyMin,pictureSize,in,outImage,outEndInfo)
canDivideorNot=1;
minDist=8;
for i=1:numel(in)
    bacteriaImage{i}=[];
    afterDivide{i}=[];
    if ~isempty(in{i}.BWImage)
        model=imfill(in{i}.BWImage,'holes');
        if ~isempty(afterDividePre{i})
            properNum(i,1)=5;
            inEndInfo=[];
            model=imrotate(model,afterDividePre{i}.sitaNum);
            model=bwToOneArea(model);
            image=afterDividePre{i}.skelImage;
            image=bwToOneArea(image);
            image=imfill(image,'holes');
            image=bwmorph(image,'thin',inf);
%             image1=image-bwmorph(image,'endPoints');
%             if max(max(image1))~=0
%                 image=bwmorph(image1,'endPoints');
%             else
                image=bwmorph(image,'endPoints');
%             end
            [xInfo,yInfo]=find(image==1);
            inEndInfo(:,1)=xInfo;
            inEndInfo(:,2)=yInfo;
            pointNum=findpointInImage(outEndInfo,afterDividePre{i}.Image);
            image=false(size(outImage));
            if isempty(pointNum)
                distanceMat=pdist2(inEndInfo,outEndInfo);
                [inNum,outNum]=find(distanceMat==min(min(distanceMat)));
                image(xInfo(inNum(1)),yInfo(inNum(1)))=1;
                way=1;
            else
                distanceMat=pdist2(inEndInfo,outEndInfo(pointNum,:));
                [inNum,outNum]=find(distanceMat==min(min(distanceMat)));
                image(xInfo(inNum(1)),yInfo(inNum(1)))=1;
                way=2;
            end
            inBWImage=imfill(in{i}.BWImage,'holes');
            if numel(model(model==1))<500
                cc=regionprops(model,'MajorAxisLength','MinorAxisLength','Orientation');
                cc=cc(1);
                model=getBacteria(cc.MajorAxisLength,cc.MinorAxisLength,cc.Orientation);
                image=imdilate(image,model);
                image=image&outImage;
                if -numel(image(image==1))+numel(inBWImage(inBWImage==1))>numel(inBWImage(inBWImage==1))*0.4
                    properNum(i,1)=5;
                    canDivideorNot=0;
                    return
                end
                image=bwmorph(image,'open');
            else
                image=imdilate(image,model);
                image=image&outImage;
                if -numel(image(image==1))+numel(inBWImage(inBWImage==1))>numel(inBWImage(inBWImage==1))*0.4
                    properNum(i,1)=5;
                    canDivideorNot=0;
                    return
                end
            end
            image=bwToOneArea(image);
            if max(max(image))==1
                outInfo=regionprops(image,'Centroid','PixelIdx','MajorAxisLength','MinorAxisLength');
                outInfo=outInfo(1);
%                 if numel(outInfo)==2
%                     p=1;
%                 end
                outImageInfo=[outInfo.Centroid(2)+xyMin(1),outInfo.Centroid(1)+xyMin(2),outInfo.MajorAxisLength,outInfo.MinorAxisLength];
                afterDividePre{i}.minDist=min(pdist2(afterDividePre{i}.inInfo(1:2),outImageInfo(1:2)),pdist2(afterDividePre{i}.inInfo(3:4),outImageInfo(3:4)));
                if afterDividePre{i}.inInfo(3)<20 && afterDividePre{i}.inInfo(4)<13
                    bestFitThr=3;
                else
                    bestFitThr=2;
                end
                if afterDividePre{i}.minDist<=bestFitThr;
                    properNum(i,1)=1;
                else
                    if afterDividePre{i}.minDist<=minDist
                        properNum(i,1)=2;
                    else
                        properNum(i,1)=3;
                    end
                end
                bacteriaImage{i}=image;
                image=bwmorph(image,'remove');
                afterDivide{i}=xy2Idx(xyMin,image,pictureSize);
                [afterDivide,properNum]=overlapJudge(bacteriaImage,afterDividePre,afterDivide,properNum);
                if way==1
                    outEndInfo(outNum(1),:)=[];
                else if way==2
                        outEndInfo(pointNum(outNum(1)),:)=[];
                    end
                end
                continue
            else
                canDivideorNot=0;
                return
            end
        else
            properNum(i,1)=4;
        end
    end
end
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
function pointNum=findpointInImage(outEndInfo,BWImage)
pointNum=[];
BWImage=bwmorph(BWImage,'dilate');
for i=1:size(outEndInfo,1)
   image=false(size(BWImage));
   image(outEndInfo(i,1),outEndInfo(i,2))=1;
   if max(max(image & BWImage))==1
       pointNum=[pointNum;i];
   end
end
end
function newA=gainSameSize(A,B,pictureSize)
A.BWImage=imfill(A.BWImage,'holes');
% change picture in A to the same size that in B
pixelIdxList=xy2Idx(A.xyMin,A.BWImage,pictureSize);
new=false(pictureSize);
new(pixelIdxList)=true;
bSize=size(B.BWImage);
newA=false(bSize);
newA(2:end-1,2:end-1)=new(B.xyMin(1):B.xyMin(1)+bSize(1)-3,B.xyMin(2):B.xyMin(2)+bSize(2)-3);
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
function model=getBacteria(length,width,oritation)
% 100Neo
% a=length/2;
% b=width/2;
% if b<=5
%     b=5;
% end
% length=sqrt(1-5/(b.^2))*a*2;
% length=ceil(length);
% if length>30
% width=11;
% else
%     width=max(fix(width),11);
% end
%
a=length/2;
b=width/2;
if b<=5
    b=5;
end
length=sqrt(1-1/(b.^2))*a*2;
length=ceil(length);
if length>23
width=11;
else
    width=max(fix(width),11);
end
model=false(width,length);
model(3:width-2,1)=true;
model(3:width-2,length)=true;
model(2:width-1,2)=true;
model(2:width-1,length-1)=true;
model(:,3:length-2)=true;
model=imrotate(model,oritation);
% model=bwmorph(model,'remove');
end
function [afterDivide,properNum]=overlapJudge(bacteriaImage,afterDividePre,afterDivide,properNum)
rightNum=size(afterDivide,2);
for i=1:rightNum-1
    if ~isempty(bacteriaImage{i}) 
        overlapArea=bacteriaImage{i} & bacteriaImage{rightNum};
        if numel(overlapArea(overlapArea==1))/min(numel(bacteriaImage{i}(bacteriaImage{i}==1)),numel(bacteriaImage{i}(bacteriaImage{rightNum}==1)))>=0.8
            if afterDividePre{i}.minDist>=afterDividePre{rightNum}.minDist
                afterDivide{i}=[];
                properNum(i)=4;
                bacteriaImage{i}=[];
            else
                afterDivide{rightNum}=[];
                properNum(rightNum)=4;
                bacteriaImage{rightNum}=[];
                return
            end
        end
    end
end
end
                