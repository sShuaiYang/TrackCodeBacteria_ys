function [glueMask,glueMaskReal,glueMaskSearching]=getGlueMask(bioTree,startFrame,stepFrame,endFrame)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
allSize=xSize*ySize;
timeMap=startFrame:stepFrame:endFrame;
iImage=1:size(timeMap,2);
glueMask=false(xSize,ySize,size(timeMap,2));
glueMaskTemp=false(xSize,ySize,size(timeMap,2));
glueMaskSearching=false(xSize,ySize,size(timeMap,2));
coreBranch=[3,4,6];
for iframe=1:size(bioTree,2)
    dispFrame(iframe);
    if iframe>endFrame
        break;
    end
    if ~isempty(bioTree{iframe}.root)
        for  iroot=1:size(bioTree{iframe}.root,2)
            rootInfo=[iframe,iroot];
            if bioTree{iframe}.root{iroot}.is2Node==true
                if ismember(bioTree{iframe}.root{iroot}.branchIndex,coreBranch)
                    for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2);
                        if iframe+iTrace-1<=endFrame
                            pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{iTrace};
                            if ~isempty(iImage(timeMap==iframe+iTrace-1))
                                glueMaskTemp(pixelIdxList+(iImage(timeMap==iframe+iTrace-1)-1)*allSize)=1;
                            end
                            if iframe~=1
                                pixelIdxList1=transferFrame1(pixelIdxList,iframe,iTrace,timeMap,allSize,size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2));
                                glueMaskSearching(pixelIdxList1)=1;
                            end
                            pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize);
                            glueMask(pixelIdxList)=1;
                            
                        else
                            continue;
                        end
                    end
                end
            end
            if bioTree{iframe}.root{iroot}.is2Node==false
                if ismember(bioTree{iframe}.root{iroot}.branchIndex,coreBranch)
                    for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2);
                        if iframe+iTrace-1<=endFrame
                            pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{iTrace};
                            if ~isempty(iImage(timeMap==iframe+iTrace-1))
                                glueMaskTemp(pixelIdxList+(iImage(timeMap==iframe+iTrace-1)-1)*allSize)=1;
                            end
                            if iframe~=1
                                pixelIdxList1=transferFrame1(pixelIdxList,iframe,iTrace,timeMap,allSize,size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2));
                                glueMaskSearching(pixelIdxList1)=1;
                            end
                            pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize);
                            glueMask(pixelIdxList)=1;
                        else
                            continue;
                        end
                    end
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            nodeInfo=[iframe,iNode];
            if ismember(bioTree{iframe}.node{iNode}.branchIndex,coreBranch)
                for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
                    for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList,2);
                        if iframe+iTrace-1<=endFrame
                            pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{iTrace};
                            if ~isempty(iImage(timeMap==iframe+iTrace-1))
                                glueMaskTemp(pixelIdxList+(iImage(timeMap==iframe+iTrace-1)-1)*allSize)=1;
                            end
                            pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize);
                            glueMask(pixelIdxList)=1;
                        else
                            continue;
                        end
                    end
                end
            end
        end
    end
end
for iImage=1:size(timeMap,2)-1
    glueMask(:,:,iImage+1)=glueMask(:,:,iImage+1)|glueMask(:,:,iImage);
end
glueMaskReal=glueMask;
for iImage=1:size(timeMap,2)
    glueMaskTemp(:,:,iImage)=imfill(glueMaskTemp(:,:,iImage),'holes');
%     glueMask(:,:,iImage)=glueMask(:,:,iImage)&(glueMaskSearching(:,:,iImage));
    glueMaskSearching(:,:,iImage)=imfill(glueMaskSearching(:,:,iImage),'holes');
    glueMask(:,:,iImage)=glueMask(:,:,iImage)&(~glueMaskTemp(:,:,iImage));
end
for iImage=1:size(timeMap,2)-1
    glueMaskSearching(:,:,iImage+1)=glueMaskSearching(:,:,iImage+1)|glueMaskSearching(:,:,iImage);
end
for iImage=1:size(timeMap,2)
    glueMaskSearching(:,:,iImage)=bwmorph(glueMaskSearching(:,:,iImage),'remove');
end
end
function pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize)
iImage=1:size(timeMap,2);
iMax=max(iImage(timeMap<=iframe+iTrace-1));
pixelIdxList=pixelIdxList+(iMax-1)*allSize;
end
function pixelIdxListAfter=transferFrame1(pixelIdxList,iframe,iTrace,timeMap,allSize,durationTime)
pixelIdxListAfter=[];
iImage=1:size(timeMap,2);
iMax1=max(iImage(timeMap<=iframe+iTrace-1));
iMax2=max(iImage(timeMap<=iframe+durationTime-1));
allI=iImage(iImage>=iMax1 & iImage<=iMax2);
for i=1:size(allI,2)
pixelIdxListAfter=[pixelIdxListAfter;pixelIdxList+(allI(i)-1)*allSize];
end
end
function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
