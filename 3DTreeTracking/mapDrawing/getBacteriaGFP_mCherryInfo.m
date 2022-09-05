function [markPicture,correctRate]=getBacteriaGFP_mCherryInfo(bioTree,gfpImage,rfpImage,bacteriaFrameInfo,branchList,stackNum)
% 同时拍摄GFP和RFP两个channel，用来验证程序准确性的一个程序
% stackNum=95;
% gfpImage=gfpImage(:,:,1:stackNum);
% rfpImage=rfpImage(:,:,1:stackNum);
imageSize=bioTree{1}.imageSize;
% [bioTree,branchList,~,~]=myBiograph_new2(bioTree);

%  initialize bioTree
bioTree=bioTreeInitialize(bioTree,branchList);

% processing gfpImage
greenBack=100;
redBack=100;
gfpImage=gfpImage-greenBack;
rfpImage=rfpImage-redBack;
% gfpImage=processingFluoImage(gfpImage);
% rfpImage=processingFluoImage(rfpImage);

% mark bioTree branch with R & G
greenThresh=250;
redThresh=50;
smallStackNum=100;
branchOrder=[];
orderNum=[1,100:smallStackNum:(stackNum-1)*smallStackNum];
for iframe=orderNum
    focusGFPImage=gfpImage(:,:,floor((iframe+1)/smallStackNum)+1);
    focusGFPImage(874,:)=0;
    focusGFPImage(:,865)=0;
    focusRFPImage=rfpImage(:,:,floor((iframe+1)/smallStackNum)+1);
    focusRFPImage(874,:)=0;
    focusRFPImage(:,865)=0;
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    bacteriaInfo=bacteriaInfo(~ismember(bacteriaInfo(:,5),branchOrder),:);
    branchOrder=[branchOrder;bacteriaInfo(:,5)];
    branchOrder=unique(branchOrder);
    for ibac=1:size(bacteriaInfo,1)
        bacInfo=bacteriaInfo(ibac,:);
        if bacInfo(3)==0
            pixelIdxInfo=bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)};
        else
            pixelIdxInfo=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)};
        end
        [xyMin,BWImage]=idx2Xy(pixelIdxInfo,imageSize);
        pixelIdxInfo=xy2Idx(xyMin,BWImage,imageSize);
        bacBranch=bacInfo(5);
        branchInfo=branchList(bacBranch,:);
        if branchInfo(3)==0
            if bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo==0
                if mean(focusGFPImage(pixelIdxInfo))>=greenThresh
                    bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo=1;
                end
                if mean(focusRFPImage(pixelIdxInfo))>=redThresh
                    bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo=2;
                end
            end
        end
        if branchInfo(3)==1
            if bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo==0
                if mean(focusGFPImage(pixelIdxInfo))>=greenThresh
                    bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo=1;
                end
                if mean(focusRFPImage(pixelIdxInfo))>=redThresh
                    bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo=2;
                end
            end
        end
    end
end

% get the final color picture
% finalPage=orderNum(end);
finalPage=100*stackNum;
focusGFPImage=gfpImage(:,:,stackNum+1);
focusGFPImage(874,:)=0;
focusGFPImage(:,865)=0;
focusRFPImage=rfpImage(:,:,stackNum+1);
focusRFPImage(874,:)=0;
focusRFPImage(:,865)=0;
markPicture1=uint8(zeros(imageSize));
markPicture2=uint8(zeros(imageSize));
markPicture3=uint8(zeros(imageSize));
bacteriaInfo=bacteriaFrameInfo{finalPage}.bacteriaInfo;
for ibac=1:size(bacteriaInfo,1)
    bacInfo=bacteriaInfo(ibac,:);
    if bacInfo(3)==0
        pixelIdxInfo=bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)};
    else
        pixelIdxInfo=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)};
    end
    bacBranch=bacInfo(5);
    branchInfo=branchList(bacBranch,:);
    if branchInfo(3)==0
        if bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo==0
            [xyMin,BWImage]=idx2Xy(pixelIdxInfo,imageSize);
            pixelIdxInfo=xy2Idx(xyMin,BWImage,imageSize);
            if mean(focusGFPImage(pixelIdxInfo))<greenThresh
                    markPicture1(pixelIdxInfo)=255;
                else
                    markPicture2(pixelIdxInfo)=255;
            end
            %             markPicture1(pixelIdxInfo)=255;
            %             markPicture2(pixelIdxInfo)=255;
        else
            [xyMin,BWImage]=idx2Xy(pixelIdxInfo,imageSize);
            pixelIdxInfo=xy2Idx(xyMin,BWImage,imageSize);
            if bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo==1
%                 markPicture2(pixelIdxInfo)=255;  % add 5-22
                if mean(focusGFPImage(pixelIdxInfo))<greenThresh
                    markPicture1(pixelIdxInfo)=255;
                else
                    markPicture2(pixelIdxInfo)=255;
                end
%                 close all;imshow(markPicture2);
            end
            if bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo==2
                if mean(focusGFPImage(pixelIdxInfo))>greenThresh
                    markPicture2(pixelIdxInfo)=255;
                else
                    markPicture1(pixelIdxInfo)=255;
                end
%                 markPicture1(pixelIdxInfo)=255;
%                 close all;imshow(markPicture1)
            end
        end
    end
    if branchInfo(3)==1
        if bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo==0
            [xyMin,BWImage]=idx2Xy(pixelIdxInfo,imageSize);
            pixelIdxInfo=xy2Idx(xyMin,BWImage,imageSize);
            if mean(focusGFPImage(pixelIdxInfo))<greenThresh
                    markPicture1(pixelIdxInfo)=255;
                else
                    markPicture2(pixelIdxInfo)=255;
            end
                
%             markPicture1(pixelIdxInfo)=255;
%             markPicture2(pixelIdxInfo)=255;
        else
            [xyMin,BWImage]=idx2Xy(pixelIdxInfo,imageSize);
            pixelIdxInfo=xy2Idx(xyMin,BWImage,imageSize);
            if bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo==1
%                 markPicture2(pixelIdxInfo)=255;
                if mean(focusGFPImage(pixelIdxInfo))<greenThresh
                    markPicture1(pixelIdxInfo)=255;
                else
                    markPicture2(pixelIdxInfo)=255;
                end
%                 close all;imshow(markPicture2);
            end
            if bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo==2
%                 markPicture1(pixelIdxInfo)=255;
                if mean(focusGFPImage(pixelIdxInfo))>greenThresh
                    markPicture2(pixelIdxInfo)=255;
                else
                    markPicture1(pixelIdxInfo)=255;
                end
%                 close all;imshow(markPicture1);
            end
        end
    end
end
% markPicture1(focusRFPImage>=redThresh)=255;
% markPicture2(focusGFPImage>=greenThresh & focusRFPImage<redThresh)=255;
markPicture=cat(3,markPicture1,markPicture2,markPicture3);
close all
imshow(markPicture)

% newPicture1=uint8(zeros(imageSize));
% newPicture2=uint8(zeros(imageSize));
% newPicture3=uint8(zeros(imageSize));
% newPicture1(focusRFPImage>=redThresh)=255;
% newPicture2(focusGFPImage>=greenThresh & focusRFPImage<redThresh)=255;
% newPic=cat(3,newPicture1,newPicture2,newPicture3);
% figure;imshow(newPic)
% caculate correctRate
correctRate=bioTreeCorrectRate(bioTree,markPicture,bacteriaFrameInfo,finalPage);
end
function bioTree=bioTreeInitialize(bioTree,branchList)
for iBra=1:size(branchList,1)
    branchInfo=branchList(iBra,:);
    if branchInfo(3)==0
        bioTree{branchInfo(1)}.root{branchInfo(2)}.colorInfo=0;
    else
        bioTree{branchInfo(1)}.node{branchInfo(2)}.colorInfo=0;
    end
end
end
function correctRate=bioTreeCorrectRate(bioTree,markPicture,bacteriaFrameInfo,frameInfo)
bacteriaInfo=bacteriaFrameInfo{frameInfo}.bacteriaInfo;
bacteriaNum=size(bacteriaInfo,1);
markPicture1=markPicture(:,:,1);
markPicture2=markPicture(:,:,2);
markPicture3=markPicture(:,:,3);
rightNum=0;
for ibac=1:bacteriaNum
    bacInfo=bacteriaInfo(ibac,:);
    if bacInfo(3)==0
        pixelIdxInfo=bioTree{bacInfo(1)}.root{bacInfo(2)}.traceInfo.pixelIdxList{bacInfo(4)};
    else
        pixelIdxInfo=bioTree{bacInfo(1)}.node{bacInfo(2)}.Out{bacInfo(3)}.traceInfo.pixelIdxList{bacInfo(4)};
    end
    [xyMin,BWImage]=idx2Xy(pixelIdxInfo,bioTree{1}.imageSize);
    pixelIdxInfo=xy2Idx(xyMin,BWImage,bioTree{1}.imageSize);
    if all(markPicture2(pixelIdxInfo)==255)
        focusGFPInfo=markPicture3(pixelIdxInfo);
        if double(numel(focusGFPInfo(focusGFPInfo==255)))/numel(focusGFPInfo)>0.2
            rightNum=rightNum+1;
        end
    end
    if all(markPicture2(pixelIdxInfo)==0)
        focusGFPInfo=markPicture3(pixelIdxInfo);
        if double(numel(focusGFPInfo(focusGFPInfo==255)))/numel(focusGFPInfo)<0.2
            rightNum=rightNum+1;
        end
    end
end
correctRate=rightNum/bacteriaNum;
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
% BWImageGain=bwmorph(BWImageGain,'dilate');
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
function gfpImage=processingFluoImage(gfpImage)
gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
for i=1:size(gfpImage,3)
     gfpImage(:,:,i)=imfilter(gfpImage(:,:,i),gaussianFilter); % use Guassian blur filter process
    gfpImage(:,:,i)=imfilter(gfpImage(:,:,i),edgeFilter); %use edgeFilter process
    gfpImage(:,:,i)=imfilter(gfpImage(:,:,i),gaussianFilter); % use Guassi
    gfpImage(:,:,i)=imclearborder(gfpImage(:,:,i));
end
gfpImage=uint8(gfpImage);
for i=1:size(gfpImage,3)
    gfpImage(:,:,i)=imadjust(gfpImage(:,:,i),[double(min(min(gfpImage(:,:,i))))/255,double(max(max(gfpImage(:,:,i))))/255],[0,1]);
end
end