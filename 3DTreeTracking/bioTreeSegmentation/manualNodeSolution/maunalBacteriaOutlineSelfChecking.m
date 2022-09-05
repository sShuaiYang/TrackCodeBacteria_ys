function maunalBacteriaOutlineSelfChecking()
dirFile=uigetdir();
dirMap=[dirFile,'\manual result'];
dirResult=[dirFile,'\checking result'];
mkdir(dirResult);
nameList=dir(dirMap);
for i=1:numel(nameList)-2
    image=import_tiff_stack([dirMap,'\',nameList(i+2).name]);
    image=uint8(image);
    imageStack(:,:,:,i)=image;
    bwProfile(:,:,i)=im2bw(image(:,:,3));
end
for i=1:numel(nameList)-2
    imwrite(im2uint8(bwProfile(:,:,i)),[dirResult,'\',num2str(i),'.tif'])
end
return
for i=1:size(bwProfile,3)
    image=bwProfile(:,:,i);
    maskImageBig(:,:,i)=logical(image);
    image1=imfill(image,'holes');
    image=image1-image;
    maskImage(:,:,i)=logical(image);
end
maskImage=imageTrans(maskImage,maskImageBig);
for i=1:size(bwProfile,3)
    iImage=imageStack(:,:,:,i);
    image1=iImage(:,:,1);
    image2=iImage(:,:,2);
    image3=iImage(:,:,3);
    image1(maskImage(:,:,i))=255;
    image2(maskImage(:,:,i))=255;
    image3(maskImage(:,:,i))=255;
    iImage=cat(3,image1,image2,image3);
    imwrite(iImage,[dirResult,'\',num2str(i),'.tif'])
end
end
function imageStack= import_tiff_stack( fname )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
warning off all
infoImage=imfinfo(fname);
frameNum=size(infoImage,1);
imageWith=infoImage(1).Width;
imageHeight=infoImage(1).Height;
imageBit=infoImage(1).BitsPerSample;
imageBit=strcat('uint',num2str(imageBit));
imageStack=zeros(imageHeight,imageWith,3,1);
imageCurrent=Tiff(fname,'r');
for iframe=1:frameNum
    imageCurrent.setDirectory(iframe);
    imageStack(:,:,:,iframe)= imageCurrent.read();
end
imageCurrent.close();
end
function image2= imageTrans(image2,imageBig)
for i=1:size(image2,3)
    image=image2(:,:,i);
    image=bwareaopen(image,5,4);
    cc=bwconncomp(image,4);
    imageNew=false(size(image));
    for iCC=1:cc.NumObjects
        [xyMin,image1]=idx2Xy(cc.PixelIdxList{iCC},cc.ImageSize);
        image1=imerode(image1,[0,1,0;1,1,1;0,1,0]);
        image1=imdilate(image1,[0,1,0;1,1,1;0,1,0]);
        image1=bwmorph(image1,'thin','inf');
        image1=bwareaopen(image1,1);
        image1=imdilate(image1,ones(3));
        pixelIdxList=xy2Idx(xyMin,image1,cc.ImageSize);
        imageNew(pixelIdxList)=true;
    end
    imageNew(imageBig(:,:,i))=false;
    image2(:,:,i)=imageNew;
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
