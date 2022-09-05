function bioTree=bioTreeSizeTransform(bioTree,xSize,ySize)
% 将bioTree的视野范围进行任意的缩小的程序。xSize,ySize为指定的x,y的新范围
% 这样就可以方便的得到不同大小视野的数据

% 先把bioTree里面所有的pixelIdxList转化为row column坐标，便于鉴别细菌是否在边界内
bioTree=bioTreePixelTrans(bioTree);

for iframe=1:size(bioTree,2)
    deleteRootList=[];
    deleteLeafList=[];
    for iRoot=1:size(bioTree{iframe}.root,2)
        is2Node=bioTree{iframe}.root{iRoot}.is2Node;
        tracePixel=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList;
        n1=0;
        n2=0;
        for iTrace=1:numel(tracePixel)
            rcInfo=tracePixel{iTrace};
            if any(rcInfo(:,1)<xSize(1)) || any(rcInfo(:,1)>xSize(2)) || any(rcInfo(:,2)<ySize(1)) || any(rcInfo(:,2)>ySize(2))
                n1=iTrace;
                break
            end
        end
        if n1~=0
            for iTrace=n1:numel(tracePixel)
                rcInfo=tracePixel{iTrace};
                if all(rcInfo(:,1)>=xSize(1)) && all(rcInfo(:,1)<=xSize(2)) && all(rcInfo(:,2)>=ySize(1)) && all(rcInfo(:,2)<=ySize(2))
                    n2=iTrace;
                    break
                end
            end
        end
        
        % n1=0,n2=0 说明trace的在视野中，不用处理
        % n1~=0,n2=0 说明出去之后就再也没有回来，前置的root得到leaf，后置的删除进入的节点（node 删除in，leaf整个都删除）
        % n1~=0,n2~=0 说明出去之后又回来了，生成新的root，前置的root得到leaf
        if n1==0 && n2==0
            continue
        end
        if n1~=0 && n2==0
            rootDetail=bioTree{iframe}.root{iRoot};
            if n1==1
                deleteRootList=[deleteRootList;iframe,iRoot];
            else
                leafFrame=iframe+n1-1;
                leafInfo=[leafFrame,size(bioTree{leafFrame}.leavies,2)+1];
                bioTree{iframe}.root{iRoot}.is2Node=0;
                bioTree{iframe}.root{iRoot}.nodeInfo=[];
                bioTree{iframe}.root{iRoot}.leafInfo=leafInfo;
                bioTree{iframe}.root{iRoot}.traceInfo.pixelIndexList=bioTree{iframe}.root{iRoot}.traceInfo.pixelIndexList(1:n1-1);  % 生成root的新leaf信息
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=0;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iRoot];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{iframe}.root{iRoot}.traceInfo.pixelIndexList{end};  % 生成新的leaf
            end
            if n2==0
                if rootDetail.is2Node==0  % 指向的是leaf就删除leaf
                    deleteLeafList=[deleteLeafList;rootDetail.leafInfo];
                end
                if rootDetail.is2Node==1  % 指向的是node就删除一个in
                    nodeInfo=rootDetail.nodeInfo;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}=[];
                end
            end
            if n2~=0
                if rootDetail.is2Node==0 % 指向的是leaf就生成一个root
                    
                    
            
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
xyMin=[xMin,yMin];
end
function bioTree=bioTreePixelTrans(bioTree)
imageSize=bioTree{1}.imageSize;
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        bioTree{iRoot}.rootPixelDetail=changeIndex_to_rcPixel(bioTree{iRoot}.rootPixelDetail,imageSize);
        for iTrace=1:size(bioTree{iRoot}.traceInfo.pixelIdxList,2)
            bioTree{iRoot}.traceInfo.pixelIdxList{iTrace}=changeIndex_to_rcPixel(bioTree{iRoot}.traceInfo.pixelIdxList{iTrace},imageSize);
        end
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out)
            for iTrace=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList,2)
                bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}=changeIndex_to_rcPixel(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace},imageSize);
            end
        end
    end
    for iLeaf=1:size(bioTree{iframe}.leavies,2)
        bioTree{iframe}.leavies{iLeaf}.leaviesPixelDetail=changeIndex_to_rcPixel(bioTree{iframe}.leavies{iLeaf}.leaviesPixelDetail);
    end
end
end
function rcPixel=changeIndex_to_rcPixel(pixelIdxList,imageSize)
% change the pixelIdxInfo to row column coordinate
rSize=imageSize(1);
cresult=ceil(pixelIdxList/rSize);
rresult=round(pixelIdxList-(cresult-1)*xSize);
rcPixel=cat(2,rresult,cresult);
end
