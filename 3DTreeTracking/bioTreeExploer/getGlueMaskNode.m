function glueMaskNode=getGlueMaskNode(bioTree,startFrame,stepFrame,endFrame)
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
allSize=xSize*ySize;
timeMap=startFrame:stepFrame:endFrame;
iImage=1:size(timeMap,2);
glueMaskNode=false(xSize,ySize,size(timeMap,2));
glueMaskTemp=false(xSize,ySize,size(timeMap,2));
for iframe=1:size(bioTree,2)
    dispFrame(iframe);
    if iframe>endFrame
        break;
    end
    if ~isempty(bioTree{iframe}.root)
        for  iroot=1:size(bioTree{iframe}.root,2)
            rootInfo=[iframe,iroot];
            if bioTree{iframe}.root{iroot}.is2Node==true
                for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2);
                    if iframe+iTrace-1<=endFrame
                        pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{iTrace};
                        if ~isempty(iImage(timeMap==iframe+iTrace-1))
                            glueMaskTemp(pixelIdxList+(iImage(timeMap==iframe+iTrace-1)-1)*allSize)=1;
                        end
                        pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize);                       
                        glueMaskNode(pixelIdxList)=1;               
                    else
                        continue;
                    end
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            nodeInfo=[iframe,iNode];
            for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
                for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList,2);
                    if iframe+iTrace-1<=endFrame
                        pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{iTrace};
                        if ~isempty(iImage(timeMap==iframe+iTrace-1))
                            glueMaskTemp(pixelIdxList+(iImage(timeMap==iframe+iTrace-1)-1)*allSize)=1;
                        end
                        pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize);
                        glueMaskNode(pixelIdxList)=1;
                    else
                        continue;
                    end
                end
            end
        end
    end
end
for iImage=1:size(timeMap,2)
    glueMaskNode(:,:,iImage)=imfill(glueMaskNode(:,:,iImage),'holes');
    glueMaskNode(:,:,iImage)=bwmorph(glueMaskNode(:,:,iImage),'skel','Inf');
end
for iImage=1:size(timeMap,2)-1
    glueMaskNode(:,:,iImage+1)=glueMaskNode(:,:,iImage+1)|glueMaskNode(:,:,iImage);
end
% for iImage=1:size(timeMap,2)
%     glueMaskNode(:,:,iImage)=bwmorph(glueMaskNode(:,:,iImage),'skel','Inf');
% end
end
function pixelIdxList=transferFrame(pixelIdxList,iframe,iTrace,timeMap,allSize)
iImage=1:size(timeMap,2);
iMax=max(iImage(timeMap<=iframe+iTrace-1));
pixelIdxList=pixelIdxList+(iMax-1)*allSize;
end
function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end