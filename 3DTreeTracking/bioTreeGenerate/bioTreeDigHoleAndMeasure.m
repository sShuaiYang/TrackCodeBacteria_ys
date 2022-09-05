function bioTree=bioTreeDigHoleAndMeasure(bioTree,beginFrame,xSize,ySize)
parfor iframe=beginFrame:size(bioTree,2)
   if ~isempty(bioTree{iframe}.root)
       bioTree{iframe}.root=DigRoot(bioTree{iframe}.root,xSize,ySize);
   end
   if ~isempty(bioTree{iframe}.node)
       bioTree{iframe}.node=DigNode(bioTree{iframe}.node,xSize,ySize);
   end
   if ~isempty(bioTree{iframe}.leavies)
       bioTree{iframe}.leavies=DigLeaf(bioTree{iframe}.leavies,xSize,ySize);
   end
end
end
function bioRoot=DigRoot(bioRoot,xSize,ySize)
imageSize=[xSize,ySize];
for iRoot=1:size(bioRoot,2)
    pixelIdxList=bioRoot{iRoot}.rootPixelDetail;
    [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
    bioRoot{iRoot}.rootPixelDetail=pixelIdxList;
    bioRoot{iRoot}.rootMeasurment=pos;
    bioRoot{iRoot}.traceInfo.measurment=[];
    for iTrace=1:size(bioRoot{iRoot}.traceInfo.pixelIdxList,2)
        pixelIdxList=bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace};
        [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
        bioRoot{iRoot}.traceInfo.measurment{iTrace}=pos;
        bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace}=pixelIdxList;
    end
end
end
function bioLeaf=DigLeaf(bioLeaf,xSize,ySize)
imageSize=[xSize,ySize];
for iLeaf=1:size(bioLeaf,2)
    pixelIdxList=bioLeaf{iLeaf}.leaviesPixelDetail;
    [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
    bioLeaf{iLeaf}.leaviesPixelDetail=pixelIdxList;
    bioLeaf{iLeaf}.leafMeasurment=pos;
    
end
end
function bioNode=DigNode(bioNode,xSize,ySize)
imageSize=[xSize,ySize];
for iNode=1:size(bioNode,2)
    for iNodeOut=1:size(bioNode{iNode}.Out,2)
        bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment=[];
        for iNodeTrace=1:size(bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
            pixelIdxList=bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
            [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
            bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}=pos;
            bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}=pixelIdxList;
        end
    end
end
end
function [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize)
[xyMin,BWImage]=idx2Xy(pixelIdxList,imageSize);
pos=regionprops(BWImage,'FilledArea','Centroid','Eccentricity','MajorAxisLength','Orientation','MinorAxisLength');
% BWImage=imfill(BWImage,'holes');  % 原来是屏蔽的
for i=1:size(pos,1)
pos(i).Centroid(1)=pos(i).Centroid(1)+xyMin(2)-1;
pos(i).Centroid(2)=pos(i).Centroid(2)+xyMin(1)-1;
end
pixelIdxList=xy2Idx(xyMin,BWImage,imageSize);
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
BWImageGain=bwmorph(BWImageGain,'close');
% BWImageGain=bwmorph(BWImageGain,'remove');
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

% %% 对于模拟细菌的分析
% function bioTree=bioTreeDigHoleAndMeasure(bioTree,beginFrame,xSize,ySize)
% parfor iframe=beginFrame:size(bioTree,2)
%    if ~isempty(bioTree{iframe}.root)
%        bioTree{iframe}.root=DigRoot(bioTree{iframe}.root,xSize,ySize);
%    end
%    if ~isempty(bioTree{iframe}.node)
%        bioTree{iframe}.node=DigNode(bioTree{iframe}.node,xSize,ySize);
%    end
%    if ~isempty(bioTree{iframe}.leavies)
%        bioTree{iframe}.leavies=DigLeaf(bioTree{iframe}.leavies,xSize,ySize);
%    end
% end
% end
% function bioRoot=DigRoot(bioRoot,xSize,ySize)
% imageSize=[xSize,ySize];
% for iRoot=1:size(bioRoot,2)
%     pixelIdxList=bioRoot{iRoot}.rootPixelDetail;
%     [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
%     bioRoot{iRoot}.rootPixelDetail=pixelIdxList;
%     bioRoot{iRoot}.rootMeasurment=pos;
%     bioRoot{iRoot}.traceInfo.measurment=[];
%     for iTrace=1:size(bioRoot{iRoot}.traceInfo.pixelIdxList,2)
%         pixelIdxList=bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace};
% %         [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
% %         bioRoot{iRoot}.traceInfo.measurment{iTrace}=pos;
%         bioRoot{iRoot}.traceInfo.measurment{iTrace}.Centroid=bioRoot{iRoot}.traceInfo.traceCentroid(iTrace,:);
%         bioRoot{iRoot}.traceInfo.pixelIdxList{iTrace}=pixelIdxList;
%     end
% end
% end
% function bioLeaf=DigLeaf(bioLeaf,xSize,ySize)
% imageSize=[xSize,ySize];
% for iLeaf=1:size(bioLeaf,2)
%     pixelIdxList=bioLeaf{iLeaf}.leaviesPixelDetail;
%     [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
%     bioLeaf{iLeaf}.leafPixelDetail=pixelIdxList;
%     bioLeaf{iLeaf}.leafMeasurment=pos;
%     
% end
% end
% function bioNode=DigNode(bioNode,xSize,ySize)
% imageSize=[xSize,ySize];
% for iNode=1:size(bioNode,2)
%     for iNodeOut=1:size(bioNode{iNode}.Out,2)
%         bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment=[];
%         for iNodeTrace=1:size(bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
%             pixelIdxList=bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
% %             [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize);
% %             bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}=pos;
%             bioNode{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}.Centroid=bioNode{iNode}.Out{iNodeOut}.traceInfo.traceCentroid(iNodeTrace,:);
%             bioNode{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}=pixelIdxList;
%         end
%     end
% end
% end
% function [pos,pixelIdxList]=fastGetProps(pixelIdxList,imageSize)
% [xyMin,BWImage]=idx2Xy(pixelIdxList,imageSize);
% pos=regionprops(BWImage,'Centroid');
% % BWImage=imfill(BWImage,'holes');  % 原来是屏蔽的
% for i=1:size(pos,2)
% pos(i).Centroid(1)=pos(i).Centroid(1)+xyMin(2)-1;
% pos(i).Centroid(2)=pos(i).Centroid(2)+xyMin(1)-1;
% end
% pixelIdxList=xy2Idx(xyMin,BWImage,imageSize);
% end
% function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% % this function is used to convert linearIdx to a small BWImage
% xSize=pictureSize(1);
% yresult=ceil(pixelIdxList/xSize);
% xresult=round(pixelIdxList-(yresult-1)*xSize);
% xMin=min(xresult);
% xMax=max(xresult);
% yMin=min(yresult);
% yMax=max(yresult);
% yresult2=yresult-yMin+1;
% xresult2=xresult-xMin+1;
% BWImage=false(xMax-xMin+1,yMax-yMin+1);
% pixelIdxList2=round(xresult2+(yresult2-1)*(xMax-xMin+1));
% BWImage(pixelIdxList2)=true;
% BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
% BWImageGain(2:end-1,2:end-1)=BWImage;
% xyMin=[xMin,yMin];
% end
% function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
% %this function is used to convert a small BWImage to its original pixelIdx
% BWImage=BWImageGain(2:end-1,2:end-1);
% % BWImage=bwmorph(BWImage,'remove');
% pixelIdxListOri=find(BWImage==1);
% smallxSize=size(BWImage,1);
% yresult=ceil(pixelIdxListOri/smallxSize);
% xresult=pixelIdxListOri-(yresult-1)*smallxSize;
% xresult2=xresult+xyMin(1)-1;
% yresult2=yresult+xyMin(2)-1;
% xSize=pictureSize(1);
% pixelIdxList=xresult2+(yresult2-1)*xSize;
% end