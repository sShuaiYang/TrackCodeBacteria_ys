function glueMap=timeGlueMap(bioTree,endTime,stepTime) %get the serise imgae of the Glue map with diffrent type bacteria
dirSave=uigetdir();
cd(dirSave);
xSize=bioTree{1}.imageSize(1);
ySize=bioTree{1}.imageSize(2);
num=fix(endTime/stepTime);
timeStop=stepTime:stepTime:endTime;
timeType0MapNumerCount=zeros(xSize,ySize,num);
timeType1MapNumerCount=zeros(xSize,ySize,num);
timeType0MapTimeCount=zeros(xSize,ySize,num);
timeType1MapTimeCount=zeros(xSize,ySize,num);
for i=1:num
    disp(i);
    [type0MapNumerCount,type1MapNumerCount,type0MapTimeCount,type1MapTimeCount]=createGlueMap(bioTree,timeStop(i),xSize,ySize);
    timeType0MapNumerCount(:,:,i)=type0MapNumerCount;
    timeType1MapNumerCount(:,:,i)=type1MapNumerCount;
    timeType0MapTimeCount(:,:,i)=type0MapTimeCount;
    timeType1MapTimeCount(:,:,i)=type1MapTimeCount;
end
glueMap.Type0MapNumerCount=timeType0MapNumerCount;
glueMap.Type1MapNumerCount=timeType1MapNumerCount;
glueMap.Type0MapTimeCount=timeType0MapTimeCount;
glueMap.Type1MapTimeCount=timeType1MapTimeCount;
saveFile1=strcat(dirSave,'\glueMap');
save(saveFile1,'glueMap','-v7.3');
end
function [type0MapNumerCount,type1MapNumerCount,type0MapTimeCount,type1MapTimeCount]=createGlueMap(bioTree,endTime,xSize,ySize)
type0MapNumerCount=zeros(xSize,ySize);
type1MapNumerCount=zeros(xSize,ySize);
type0MapTimeCount=zeros(xSize,ySize);
type1MapTimeCount=zeros(xSize,ySize);
for iframe=1:endTime
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==false
                rootMask0=false(xSize,ySize);
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace};
                    if iTrace+iframe-1<=endTime
                        rootMask0(pixelIdxList)=1;
                    else
                        break;
                    end
                end
                temp=zeros(xSize,ySize);
                rootMask0=imfill(rootMask0,'holes');
                temp=temp+rootMask0;
                type0MapNumerCount=type0MapNumerCount+rootMask0;
                type0MapTimeCount=type0MapTimeCount+temp.*size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2);
            end
            if bioTree{iframe}.root{iroot}.is2Node==true
                rootMask1=false(xSize,ySize);
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)                    
                    pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace};
                    if iTrace+iframe-1<=endTime
                        rootMask1(pixelIdxList)=1;
                        
                    else
                        break;
                    end
                end
                temp=zeros(xSize,ySize);
                rootMask1=imfill(rootMask1,'holes');
                temp=temp+rootMask1;
                type1MapNumerCount=type1MapNumerCount+rootMask1;
                type1MapTimeCount=type1MapTimeCount+temp.*size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2);
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iNodeOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                nodeMask=false(xSize,ySize);
                for iNodeTrace=1:size(bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                    if iNodeTrace+iframe-1<=endTime
                        pixelIdxList=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace};
                        nodeMask(pixelIdxList)=1;                        
                    else
                        break;
                    end
                end
                temp=zeros(xSize,ySize);
                nodeMask=imfill(nodeMask,'holes');
                temp=temp+nodeMask;
                type1MapNumerCount=type1MapNumerCount+nodeMask;
                type1MapTimeCount=type1MapTimeCount+temp.*size(bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2);
            end
        end
    end
end
end