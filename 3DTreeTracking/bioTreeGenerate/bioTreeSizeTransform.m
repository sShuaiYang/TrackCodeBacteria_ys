function bioTree=bioTreeSizeTransform(bioTree,xSize,ySize)
% ��bioTree����Ұ��Χ�����������С�ĳ���xSize,ySizeΪָ����x,y���·�Χ
% �����Ϳ��Է���ĵõ���ͬ��С��Ұ������

% �Ȱ�bioTree�������е�pixelIdxListת��Ϊrow column���꣬���ڼ���ϸ���Ƿ��ڱ߽���
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
        
        % n1=0,n2=0 ˵��trace������Ұ�У����ô���
        % n1~=0,n2=0 ˵����ȥ֮�����Ҳû�л�����ǰ�õ�root�õ�leaf�����õ�ɾ������Ľڵ㣨node ɾ��in��leaf������ɾ����
        % n1~=0,n2~=0 ˵����ȥ֮���ֻ����ˣ������µ�root��ǰ�õ�root�õ�leaf
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
                bioTree{iframe}.root{iRoot}.traceInfo.pixelIndexList=bioTree{iframe}.root{iRoot}.traceInfo.pixelIndexList(1:n1-1);  % ����root����leaf��Ϣ
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=0;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iRoot];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=bioTree{iframe}.root{iRoot}.traceInfo.pixelIndexList{end};  % �����µ�leaf
            end
            if n2==0
                if rootDetail.is2Node==0  % ָ�����leaf��ɾ��leaf
                    deleteLeafList=[deleteLeafList;rootDetail.leafInfo];
                end
                if rootDetail.is2Node==1  % ָ�����node��ɾ��һ��in
                    nodeInfo=rootDetail.nodeInfo;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}=[];
                end
            end
            if n2~=0
                if rootDetail.is2Node==0 % ָ�����leaf������һ��root
                    
                    
            
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
