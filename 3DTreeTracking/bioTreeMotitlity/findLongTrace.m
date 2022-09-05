function [bioTree,allData]=findLongTrace(bioTree)
frameThreshold=200;
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
iData=1;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        nodeInfo=[];
        isNode=false;
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==false;
                leafInfo=bioTree{iframe}.root{iroot}.leafInfo;
                rootInfo=[iframe,iroot];
                if (leafInfo(1)-rootInfo(1))>=frameThreshold
                    [bioTree, dataList]=getAllData(bioTree,rootInfo,nodeInfo,isNode,frameThreshold);
                    if ~isempty(dataList.velocityData) && size(dataList.velocityData,1)>=200
                        allData{iData}=dataList;
                        allData{iData}.branchIndex=[];
                        iData=iData+1;
                    end;
                end
            end
        end
    end
end
for iList=1:size(branchList,1)
    nodeInfo=branchList(iList,1:3);
    if nodeInfo(3)==0
        continue
    end
    allRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allRoot;
    allNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allNode;
    branchIndex=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex;
    for iroot=1:size(allRoot,1)
        isNode=false;
        rootInfo=allRoot(iroot,1:2);
        nodeInfo1=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
        if (nodeInfo1(1)-rootInfo(1))>=frameThreshold
            [bioTree, dataList]=getAllData(bioTree,rootInfo,nodeInfo,isNode,frameThreshold);
            if ~isempty(dataList.velocityData) && size(dataList.velocityData,1)>=200
                allData{iData}=dataList;
                allData{iData}.branchIndex=branchIndex;
                iData=iData+1;
            end
        end
    end
    for iNode=1:size(allNode,1)
        isNode=true;
        nodeInfo2=allNode(iNode,1:2);
        for iOut=1:size(bioTree{nodeInfo2(1)}.node{nodeInfo2(2)}.Out,2)
            if bioTree{nodeInfo2(1)}.node{nodeInfo2(2)}.Out{iOut}.is2Node==false
                leafInfo=bioTree{nodeInfo2(1)}.node{nodeInfo2(2)}.Out{iOut}.leafInfo;
                if (leafInfo(1)-nodeInfo2(1))>=frameThreshold
                    [bioTree, dataList]=getAllData(bioTree,rootInfo,[nodeInfo2,iOut],isNode,frameThreshold);
                    if ~isempty(dataList.velocityData) && size(dataList.velocityData,1)>=200
                        allData{iData}=dataList;
                        allData{iData}.branchIndex=branchIndex;
                        iData=iData+1;
                    end
                end
            end
             if bioTree{nodeInfo2(1)}.node{nodeInfo2(2)}.Out{iOut}.is2Node==true
                nodeInfo_next=bioTree{nodeInfo2(1)}.node{nodeInfo2(2)}.Out{iOut}.nodeInfo;
                if (nodeInfo_next(1)-nodeInfo2(1))>=frameThreshold
                    [bioTree, dataList]=getAllData(bioTree,rootInfo,[nodeInfo2,iOut],isNode,frameThreshold);
                    if ~isempty(dataList.velocityData) && size(dataList.velocityData,1)>=200
                        allData{iData}=dataList;
                        allData{iData}.branchIndex=branchIndex;
                        iData=iData+1;
                    end
                end
            end
        end
    end
end

end
function [bioTree, dataList]=getAllData(bioTree,rootInfo,nodeInfo,isNode,~)
dataList=[];
p1xPosition=[];p1yPosition=[];
p2xPosition=[];p2yPosition=[];
centroidx=[];centroidy=[];
orientation=[];
majorAxisLength=[];
vectorOrientation=[];
denosieData=[];
velocityData=[];
frameInfo=[];
if isNode==false
    if ~isempty(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment)
        traceSize=size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment,2);
        if bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node==1 && traceSize>=200
            traceSize=traceSize-199;
        end
        for iTrace=1:traceSize-1
            if numel(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1))>=2
                break
            end
            frameInfo=[frameInfo,rootInfo(1)+iTrace-1];
            p1_x=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p1Position(1);
            p1_y=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p1Position(2);
            p2_x=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p2Position(1);
            p2_y=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).p2Position(2);
            p1xPosition=[ p1xPosition,p1_x];
            p1yPosition=[ p1yPosition,p1_y];
            p2xPosition=[ p2xPosition,p2_x];
            p2yPosition=[ p2yPosition,p2_y];
            centroidx=[ centroidx,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).Centroid(1)];
            centroidy=[ centroidy,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).Centroid(2)];
            orientation=[orientation,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).Orientation];
            majorAxisLength=[majorAxisLength,bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
            denosieData=[denosieData;bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1).denosieData];
            velocityData=[ velocityData;bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}(1). velocityData];
        end
        dataList.rootInfo=rootInfo;
        dataList.nodeInfo=[];
        dataList.isNode=false;
        dataList.originalData=[ frameInfo',p1xPosition',p1yPosition',p2xPosition',p2yPosition',centroidx',centroidy',orientation',majorAxisLength'];
        dataList.denosieData=[frameInfo',denosieData];
        dataList.velocityData=[frameInfo',velocityData];
        return;
    end
end
if isNode==true
    if ~isempty(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment)
        traceSize=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment,2);
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node==1 && traceSize>=200
            traceSize=traceSize-199;
        end
        for iTrace=1:traceSize-1
            if numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace})>=2
                if ~isempty(frameInfo) && numel(frameInfo)>=1
                    break
                else
                    p1xPosition=[];p1yPosition=[];
                    p2xPosition=[];p2yPosition=[];
                    centroidx=[];centroidy=[];
                    orientation=[];
                    majorAxisLength=[];
                    vectorOrientation=[];
                    denosieData=[];
                    velocityData=[];
                    frameInfo=[];
                    continue
                end
            end
            frameInfo=[frameInfo,nodeInfo(1)+iTrace-1];
            p1_x=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p1Position(1);
            p1_y=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p1Position(2);
            p2_x=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p2Position(1);
            p2_y=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).p2Position(2);
            p1xPosition=[ p1xPosition,p1_x];
            p1yPosition=[ p1yPosition,p1_y];
            p2xPosition=[ p2xPosition,p2_x];
            p2yPosition=[ p2yPosition,p2_y];
            centroidx=[ centroidx,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).Centroid(1)];
            centroidy=[ centroidy,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).Centroid(2)];
            orientation=[orientation,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).Orientation];
            majorAxisLength=[majorAxisLength,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
            denosieData=[denosieData;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).denosieData];
            velocityData=[ velocityData;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.measurment{iTrace}(1).velocityData];
        end
        dataList.rootInfo=[];
        dataList.nodeInfo=nodeInfo;
        dataList.isNode=true;
        dataList.originalData=[frameInfo',p1xPosition',p1yPosition',p2xPosition',p2yPosition',centroidx',centroidy',orientation',majorAxisLength'];
        dataList.denosieData=[frameInfo',denosieData];
        dataList.velocityData=[frameInfo',velocityData];
        return;
    end
end
end
function angle=vectorOrientation(vector) %return the angle [-180,180] of the vector
angle=atan(vector(:,2)./abs(vector(:,1))).*(180/pi);
angle(vector(:,1)<0&vector(:,2)>=0)= 180-angle1(vector(:,1)<0&vector(:,2)>=0);
angle(vector(:,1)<0&vector(:,2)<0)= -180-angle1(vector(:,1)<0&vector(:,2)<0);
end