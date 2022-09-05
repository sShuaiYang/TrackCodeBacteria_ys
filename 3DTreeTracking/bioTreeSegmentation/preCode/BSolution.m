function [afterDivide,canDivideorNot]=BSolution(in,sita,outInfo,outImage,pictureSize,xyMin)
canDivideorNot=1;
for i=1:numel(in)
    emptyDivide(i)=0;
    if ~isempty(in{i}.BWImage)
        sita=-30:2:30;
        imageStack=rotateStacks(in{i}.BWImage,sita);
        newInImage=gainSameSize(in{i},outInfo,pictureSize);
        imageStack=gainPossibleModel(imageStack,outImage,newInImage);
        inxyMin=in{i}.xyMin;
        [afterDivide{i},sitaBest]=findTheBest(in{i}.BWImage,inxyMin,imageStack,pictureSize,xyMin,sita,inf);
        
        sita=sitaBest-2:1:sitaBest+2;
        imageStack=rotateStacks(in{i}.BWImage,sita);
        imageStack=gainPossibleModel(imageStack,outImage,newInImage);
        [afterDivide{i},~]=findTheBest(in{i}.BWImage,inxyMin,imageStack,pictureSize,xyMin,sita,25);
    else
        afterDivide{i}=[];
    end   
    if ~isempty(afterDivide{i})
        emptyDivide(i)=1;
    end
end
if all(emptyDivide==0)
    afterDivide=[];
    canDivideorNot=0;
end
isN2n=isN2nMatching(in,outImage);
if isN2n==1 && any(cellfun(@isempty,afterDivide)==1)
    afterDivide=[];
    canDivideorNot=0;
    return
end
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
function imageStack=gainPossibleModel(imageStack,outImage,inImage)
% get all possible model
parfor j=1:numel(imageStack)
    neednext=1;
    testI=imerode(outImage,imageStack{j});
    cc=bwconncomp(testI);
    if cc.NumObjects==1
        image=imdilate(testI,imageStack{j});
%         image=bwmorph(image,'dilate');
        image=image&outImage;
        if max(max(image))==1
            image=bwareaopen(image,10);
            if max(max(image))~=0
                imageStack{j}=bwmorph(image,'remove');
%                 imageStack{j}=image;
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
                    imageGain=imageGain&outImage;
                    if max(max(imageGain))==1
                        imageGain=bwareaopen(imageGain,10);
                        imageGain=bwmorph(imageGain,'remove');
                        if max(max(imageGain))~=0
                            imageStack{j}=imageGain;
                            neednext=0;
                            break
                        end
                    end
                end
            end
        end
    end
    if neednext==1
        model=bwmorph(imageStack{j},'thin');
        I=imerode(outImage,model);
        I=bwmorph(I,'thin',inf);
        I=imclearborder(I);
        cc=bwconncomp(I);
        if cc.NumObjects==0
            model2=bwmorph(model,'shrink');
            model2=bwmorph(model2,'spur');
            I=imerode(outImage,model2);
            I=bwmorph(I,'thin',inf);
            I=imclearborder(I);
            cc=bwconncomp(I);
            if cc.NumObjects==1
                I=getCentroid(cc.ImageSize,cc.PixelIdxList{1});
                imageStack{j}=imdilate(I,model);
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
            else if cc.NumObjects==0
                    imageStack{j}=[];
                else if cc.NumObjects>=2
                        result=0;
                        for i=1:cc.NumObjects
                            eachone=getCentroid(cc.ImageSize,cc.PixelIdxList{i});
                            imageGain=imdilate(eachone,model);
                            imageGain=imageGain&inImage;
                            if double(numel(imageGain(imageGain==1))/numel(inImage(inImage==1)))>0.3
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
        else if cc.NumObjects==1
                I=getCentroid(cc.ImageSize,cc.PixelIdxList{1});
                imageStack{j}=imdilate(I,imageStack{j});
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
            else if cc.NumObjects>=2
                    result=0;
                    for i=1:cc.NumObjects
                        eachone=getCentroid(cc.ImageSize,cc.PixelIdxList{i});
                        imageGain=imdilate(eachone,model);
                        imageGain=imageGain&inImage;
                        if double(numel(imageGain(imageGain==1))/numel(inImage(inImage==1)))>0.3
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
                else
                    imageStack{j}=[];
                end
            end
        end
    end
end
end
function [result,sita]=findTheBest(inImage,inxyMin,imageStack,pictureSize,xyMin,sita,threShold)
% the best choice to the input bacteria
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
    sita=sita(Num);
end
end
function newA=gainSameSize(A,B,pictureSize)
A.BWImage=imfill(A.BWImage,'holes');
% change picture in A to the same size that in B
pixelIdxList=xy2Idx(A.xyMin,A.BWImage,pictureSize);
new=false(pictureSize);
new(pixelIdxList)=true;
bSize=size(B.BWImage);
newA=new(B.xyMin(1)-1:B.xyMin(1)+bSize(1)-2,B.xyMin(2)-1:B.xyMin(2)+bSize(2)-2);
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
function I=getCentroid(pictureSize,pixelIdxList)
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=pixelIdxList-(yresult-1)*xSize;
yCentroid=round(mean(yresult));
xCentroid=round(mean(xresult));
centroidLinearIndex=xCentroid+(yCentroid-1)*xSize;
I=false(pictureSize);
I(centroidLinearIndex)=1;
end