function bioTree=bioTreeMeasure(bioTree,frameShift,xSize,ySize)
parfor iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        bioTree{iframe}=measureRoot(bioTree{iframe},xSize,ySize);
    end
    if~isempty(bioTree{iframe}.leavies)
        bioTree{iframe}=measureLeaf(bioTree{iframe},xSize,ySize);
    end
    if ~isempty(bioTree{iframe}.node)
        bioTree{iframe}=measureNode(bioTree{iframe},xSize,ySize);
    end
end
end
function bioTree=measureRoot(bioTree,xSize,ySize)
% bwTemp=false(xSize,ySize);
imageSize=[xSize,ySize];
for iroot=1:size(bioTree.root,2)
    pixelIdxList=bioTree.root{iroot}.rootPixelDetail;
    pos=fastGetProps(pixelIdxList,imageSize);
    bioTree.root{iroot}.rootMeasurment=pos;
    %     bwTemp(pixelIdxList)=true;
    %     pos=myMeasurment(bwTemp);
    %     bioTree.root{iroot}.rootMeasurment=pos;
    %     bwTemp(pixelIdxList)=false;
    bioTree.root{iroot}.traceInfo.measurment=[];
    for iTrace=1:size(bioTree.root{iroot}.traceInfo.pixelIdxList,2)
        pixelIdxList=bioTree.root{iroot}.traceInfo.pixelIdxList{iTrace};
        %         bwTemp(pixelIdxList)=true;
        %         pos=myMeasurment(bwTemp);
        pos=fastGetProps(pixelIdxList,imageSize);
        bioTree.root{iroot}.traceInfo.measurment{iTrace}=pos;
%         bwTemp(pixelIdxList)=false;
    end
end
end
function bioTree=measureLeaf(bioTree,xSize,ySize)
imageSize=[xSize,ySize];
% bwTemp=false(xSize,ySize);
for iLeaf=1:size(bioTree.leavies,2)
    pixelIdxList=bioTree.leavies{iLeaf}.leaviesPixelDetail;
    %     bwTemp(pixelIdxList)=true;
    %     pos=myMeasurment(bwTemp);
    pos=fastGetProps(pixelIdxList,imageSize);
    bioTree.leavies{iLeaf}.leafMeasurment=pos;
    %     bwTemp(pixelIdxList)=false;
end
end
function bioTree=measureNode(bioTree,xSize,ySize)
imageSize=[xSize,ySize];
% bwTemp=false(xSize,ySize);
for iNode=1:size(bioTree.node,2)
    for iNodeOut=1:size(bioTree.node{iNode}.Out,2)
        bioTree.node{iNode}.Out{iNodeOut}.traceInfo.measurment=[];
        for iNodeTrace=1:size(bioTree.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
            pixelIdxList=bioTree.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
            %             bwTemp(pixelIdxList)=true;
            %             pos=myMeasurment(bwTemp);
            pos=fastGetProps(pixelIdxList,imageSize);
            bioTree.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}=pos;
            %             bwTemp(pixelIdxList)=false;
        end
    end
end
end
function pos=myMeasurment(bwTemp)
pos=regionprops(bwTemp,'FilledArea','Centroid','Eccentricity','MajorAxisLength','Orientation','MinorAxisLength');
end % Here you can set your measurement