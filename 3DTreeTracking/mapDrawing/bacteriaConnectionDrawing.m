function bacteriaConnectionDrawing(bioTree,clusterTree)
% 画出一个branch如何与其他细菌相连的示意图
dirFile=uigetdir();
dirImage=strcat(dirFile,'\tiff2matlab');
nameList=dir(dirImage);
coreBranch=[3,10,6,10];
colorB(1,:)=[1,1,0];
colorB(2,:)=[1,0,1];
colorB(3,:)=[0,1,1];
colorB(4,:)=[0,1,0];
aimNum=[52,87,37,12];
for i=201:200:numel(bioTree)
    stackNum=fix((i-1)/189)+1;
    smallNum=i-(stackNum-1)*189;
    imageStack=load(strcat(dirImage,'\',nameList(stackNum+2).name));
    imageStack=imageStack.imageStack;
    imageStack=imageStack(:,:,smallNum);
    imageStack1=imageStack;
    imageStack2=imageStack;
    imageStack3=imageStack;
    bacteriaInfo=clusterTree{i}.bacteriaList;
    distMatrix=full(clusterTree{i}.distMatrix);
    centroidInfo=[];
    % get all centroid
    for iBac=1:size(bacteriaInfo,1)
        if bacteriaInfo(iBac,3)==0
            rootInfo=bacteriaInfo(iBac,:);
            centroidInfo=[centroidInfo;bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{rootInfo(4)}.Centroid];
        end
        if bacteriaInfo(iBac,3)~=0
            nodeInfo=bacteriaInfo(iBac,:);
            centroidInfo=[centroidInfo;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{nodeInfo(4)}.Centroid];
        end
    end
    for iBranch=1:numel(coreBranch)
        coreB=coreBranch(iBranch);
        for iBac=aimNum(iBranch);
            if iBac<=size(bacteriaInfo,1)
                %         for iBac=1:size(bacteriaInfo,1)
                if bacteriaInfo(iBac,3)==0
                    rootInfo=bacteriaInfo(iBac,:);
                    pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{rootInfo(4)};
                end
                if bacteriaInfo(iBac,3)~=0
                    nodeInfo=bacteriaInfo(iBac,:);
                    pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{nodeInfo(4)};
                end
                [xyMin,BWImage]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
                pixelIdxList=xy2Idx(xyMin,BWImage,bioTree{1}.imageSize);
                if bacteriaInfo(iBac,1)+bacteriaInfo(iBac,4)-1==i
                    if bacteriaInfo(iBac,5)==coreB
                        imageStack1(pixelIdxList)=255;
                        linkBac=distMatrix(iBac,:);
                        for iLink=1:numel(linkBac)
                            if linkBac(iLink)~=0
                                if bacteriaInfo(iLink,5)~=coreB
                                    if bacteriaInfo(iLink,3)==0
                                        rootInfo=bacteriaInfo(iLink,:);
                                        pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{rootInfo(4)};
                                    end
                                    if bacteriaInfo(iLink,3)~=0
                                        nodeInfo=bacteriaInfo(iLink,:);
                                        pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{nodeInfo(4)};
                                    end
                                    [xyMin,BWImage]=idx2Xy(pixelIdxList,bioTree{1}.imageSize);
                                    pixelIdxList=xy2Idx(xyMin,BWImage,bioTree{1}.imageSize);
                                    imageStack3(pixelIdxList)=255;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    imageStack=cat(3,imageStack1,imageStack2,imageStack3);
    imshow(imageStack)
    hold on
    for iBranch=1:numel(coreBranch)
        coreB=coreBranch(iBranch);
        for iBac=aimNum(iBranch);
            if iBac<=size(bacteriaInfo,1)
                %         for iBac=1:size(bacteriaInfo,1)
                if bacteriaInfo(iBac,5)==coreB;
                    linkBac=distMatrix(iBac,:);
                    for iLink=1:numel(linkBac)
                        if linkBac(iLink)~=0
                            if bacteriaInfo(iLink,5)~=coreB
                                line([centroidInfo(iBac,1),centroidInfo(iLink,1)],[centroidInfo(iBac,2),centroidInfo(iLink,2)],'color',colorB(iBranch,:));
                            end
                        end
                    end
                end
            end
        end
    end
    saveas(gcf,strcat('C:\Users\jzy\Desktop\to Fan Jin\F3\',num2str(i),'.tif'))
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