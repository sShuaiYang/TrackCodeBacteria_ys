function glueMap=timeGlueMap_new(bioTree,endTime,stepTime) %get the serise imgae of the Glue map with diffrent type bacteria
dirSave=uigetdir();
cd(dirSave);
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
num=fix(endTime/stepTime);
timeStop=1:stepTime:endTime;
timeType0MapNumerCount=zeros(xSize,ySize,num);
timeType1MapNumerCount=zeros(xSize,ySize,num);
timeType0f1MapNumerCount=zeros(xSize,ySize,num);
for i=1:num
    disp(i);
    type1MapNumerCount=createGlueMapType1(bioTree,1,timeStop(i+1),xSize,ySize);
    timeType1MapNumerCount(:,:,i)=type1MapNumerCount;
end
glueMap.Type1MapNumerCount=timeType1MapNumerCount;
indexT=1;
for iThreshold=20:20:120
    disp(iThreshold);
    for i=1:num
        [type0MapNumerCount,type0f1MapNumerCount]=createGlueMapType0(bioTree,1,timeStop(i+1),xSize,ySize,iThreshold);
        timeType0MapNumerCount(:,:,i)=type0MapNumerCount;
        timeType0f1MapNumerCount(:,:,i)=type0f1MapNumerCount;
    end
    glueMap.Type0MapNumerCount(indexT).Type0NumberCount=timeType0MapNumerCount;
    glueMap.Type0MapNumerCount(indexT).Type0f1NumberCount=timeType0f1MapNumerCount;
    indexT=indexT+1;
end
saveFile1=strcat(dirSave,'\glueMap');
save(saveFile1,'glueMap','-v7.3');
end
function type1MapNumerCount=createGlueMapType1(bioTree,startTime,endTime,xSize,ySize)
type1MapNumerCount=zeros(xSize,ySize);
for iframe=1:endTime
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                rootMask1=false(xSize,ySize);
                countE=0;
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace};
                    if iTrace+iframe-1<=endTime && iTrace+iframe-1>=startTime
                        rootMask1(pixelIdxList)=1;
                        countE=1;
                    else
                        break;
                    end
                end
                if countE==1
                    rootMask1=imfill(rootMask1,'holes');
                    type1MapNumerCount=type1MapNumerCount+rootMask1;
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iNodeOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                nodeMask=false(xSize,ySize);
                countE=0;
                for iNodeTrace=1:size(bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                    if iNodeTrace+iframe-1<=endTime && iNodeTrace+iframe-1>=startTime
                        pixelIdxList=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
                        nodeMask(pixelIdxList)=1;
                        countE=1;
                    else
                        break;
                    end
                end
                if countE==1
                    nodeMask=imfill(nodeMask,'holes');
                    type1MapNumerCount=type1MapNumerCount+nodeMask;
                end
            end
        end
    end
end
end
function [type0MapNumerCount,type0f1MapNumerCount]=createGlueMapType0(bioTree,startTime,endTime,xSize,ySize,lengthThreshold)
type0MapNumerCount=zeros(xSize,ySize);
type0f1MapNumerCount=zeros(xSize,ySize);
for iframe=1:endTime
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==false
                rootMask0=false(xSize,ySize);
                countE=0;
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace};
                    if iTrace+iframe-1<=endTime && iTrace+iframe-1>=startTime
                        rootMask0(pixelIdxList)=1;
                        countE=1;
                    else
                        break;
                    end
                end
                if countE==1
                    leafInfo=bioTree{iframe}.root{iroot}.leafInfo;
                    rootMask0=imfill(rootMask0,'holes');
                    x0=bioTree{iframe}.root{iroot}.rootMeasurment.Centroid(1);
                    y0=bioTree{iframe}.root{iroot}.rootMeasurment.Centroid(2);
                    x1=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafMeasurment.Centroid(1);
                    y1=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leafMeasurment.Centroid(2);
                    if (x1-x0)^2+(y1-y0)^2<=lengthThreshold^2
                        type0MapNumerCount=type0MapNumerCount+rootMask0;
                    else
                        type0f1MapNumerCount=type0f1MapNumerCount+rootMask0;
                    end
                end
            end
        end
    end
end
end