function [aveIntensityT,bacteriaFrameInfo]=fluoIntensityVsTime(bioTree,bacteriaFrameInfo,gfpImage,rfpImage)
%  get gfpImage and rfpImage data
aveIntensityT(:,1)=0:100:size(bioTree,2);
imageSize=bioTree{1}.imageSize;
for iframe=100:100:size(bioTree,2)
    p=im2uint16(zeros(2160,2560));
    focusImageNum=fix(iframe/100)+1;
    gfpBackGround=getBackGround(gfpImage(:,:,focusImageNum));
    rfpBackGround=getBackGround(rfpImage(:,:,focusImageNum));
    focusGFP=gfpImage(:,:,focusImageNum)-gfpBackGround;
    focusRFP=rfpImage(:,:,focusImageNum)-rfpBackGround;
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    for ibac=1:size(bacteriaInfo,1)
        bacInfo=bacteriaInfo(ibac,:);
        if bacInfo(3)==0
            pixelIdxList=bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)};
        end
        if bacInfo(3)~=0;
            pixelIdxList=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)};
        end
        [xyMin,BWImage]=idx2Xy(pixelIdxList,imageSize);
        pixelIdxList=xy2Idx(xyMin,BWImage,imageSize);
        p(pixelIdxList)=focusGFP(pixelIdxList);
        bacteriaFrameInfo{iframe}.gfpImageInfo{ibac}=focusGFP(pixelIdxList);
        bacteriaFrameInfo{iframe}.t_DimerImageInfo{ibac}=focusRFP(pixelIdxList);
    end
end
%
for iframe=100:100:size(bioTree,2)
    focusImageNum=fix(iframe/100)+1;
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    sumIntensity=[];
    for ibac=1:size(bacteriaInfo,1)
        sumEach=mean(bacteriaFrameInfo{iframe}.t_DimerImageInfo{ibac});
        sumIntensity=[sumIntensity;sumEach];
    end
    aveIntensityT(focusImageNum,2)=mean(sumIntensity);
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
BWImageGain=bwmorph(BWImageGain,'dilate');
xyMin=[xMin,yMin];
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
% BWImage=bwmorph(BWImage,'remove');
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function backGround=getBackGround(image)
image=image(:);
image=sort(image);
imageSize=numel(image);
backGround=mean(image(1:imageSize/5));
end
