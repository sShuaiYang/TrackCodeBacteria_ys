function plotOneBranchAndMarkerWithAgeColor(bioTree,bacteriaFrameInfo)
dirFile=uigetdir();
aimBranchList=1;
color=colormap(jet(10));
for iframe=1:numel(bioTree)
    image=zeros([bioTree{1}.imageSize,3]);
    image1=image(:,:,1);
    image2=image(:,:,2);
    image3=image(:,:,3);
    for i=1:numel(aimBranchList)
        aimBranch=aimBranchList(i);
        bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
        bacteriaInfo=bacteriaInfo(bacteriaInfo(:,5)==aimBranch,:);
        
        for ibac=1:size(bacteriaInfo,1)
            bacInfo=bacteriaInfo(ibac,:);
            if bacInfo(3)==0
                pixelIdxList=bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)};
                colorOrder=max(bioTree{bacInfo(1)}.root{bacInfo(2)}.isOld1,bioTree{bacInfo(1)}.root{bacInfo(2)}.isOld2);
                [xyMin,BWImageGain]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
                pixelIdxList=xy2Idx(xyMin,BWImageGain,bioTree{1}.imageSize);
                image1(pixelIdxList)=255*color(colorOrder,1);
                image2(pixelIdxList)=255*color(colorOrder,2);
                image3(pixelIdxList)=255*color(colorOrder,3);
            end
            if bacInfo(3)~=0
                pixelIdxList=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)};
                colorOrder=max(bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.isOld1,bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.isOld2);
                [xyMin,BWImageGain]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
                pixelIdxList=xy2Idx(xyMin,BWImageGain,bioTree{1}.imageSize);
                image(pixelIdxList)=0;
                image1(pixelIdxList)=255*color(colorOrder,1);
                image2(pixelIdxList)=255*color(colorOrder,2);
                image3(pixelIdxList)=255*color(colorOrder,3);
            end
        end
    end
    image=cat(3,image1,image2,image3);
    imshow(uint8(image(559:1025,819:1371,:)));
    saveas(gcf,strcat(dirFile,'\',num2str(iframe),'.tif'))
    close all
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
pixelIdxList2=round(xresult2+(yresult2-1)*(xMax-xMin+1));
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
xyMin=[xMin,yMin];
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
