function demoMovie=treeDemoMaker(imageStack,bioTree,frameShift) % here frameShift is can be calulate by (i-1)*Size of image Stack
xSize=1024;
ySize=1024;
allSize=xSize*ySize;

rootType0Color=[0,1,0]; %green
rootType1Color=[1,0,0]; %red
leafType0Color=[0,1,1];
leafType1Color=[0,1,0];
traceType0Color=[0,1,0]; %green
traceType1Color=[1,0,0]; %red
nodeType1Color=[1,0,1];
nodeType2Color=[0,0.5,0.5];
nodeType3Color=[1,1,0];

demoImage=zeros(xSize,ySize,3,size(imageStack,3),'uint8');
tic;demoImage=makeRootType0Demo(demoImage,imageStack,bioTree,frameShift,xSize,ySize,allSize,rootType0Color);toc;
tic;demoImage=makeRootType1Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,rootType1Color);toc;
tic;demoImage=makeType0TraceDemo(demoImage,bioTree,frameShift,xSize,ySize,allSize,traceType0Color);toc;
tic;demoImage=makeType1TraceDemo(demoImage,bioTree,frameShift,xSize,ySize,allSize,traceType1Color);toc;
% demoImage=makeNodeType1Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,nodeType1Color);
% demoImage=makeNodeType2Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,nodeType2Color);
% demoImage=makeNodeType3Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,nodeType3Color);
% demoImage=makeLeafType0Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,leafType0Color);
% demoImage=makeLeafType1Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,leafType1Color);
demoMovie=immovie(demoImage);
end
function demoImage=makeRootType0Demo(demoImage,imageStack,bioTree,frameShift,xSize,ySize,allSize,rootType0Color)
disp('1.Make rootType0 demo')
rootMask=false(xSize,ySize,size(imageStack,3));
for iframe=1:size(demoImage,4)+frameShift
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.rootType==0
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    if (iTrace+iframe-2-frameShift)>=0
                        pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace}+(iTrace+iframe-2-frameShift)*allSize;
                        if iTrace==1;
                            rootMask(pixelIdxList)=1;
                            if iframe~=1+frameShift
                                rootMask(pixelIdxList-allSize)=1;
                            end
                        else
                            if iTrace+iframe-1-frameShift<=size(demoImage,4)
                                rootMask(pixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(imageStack(:,:,iframe),rootMask(:,:,iframe),rootType0Color);
end
clear rootMask;
clear imageStack;
end
function demoImage=makeRootType1Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,rootType1Color)
disp('2.Make rootType1 demo')
rootMask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)+frameShift
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.rootType==1
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    if (iTrace+iframe-2-frameShift)>=0
                        pixelIdxList=bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList{iTrace}+(iTrace+iframe-2-frameShift)*allSize;
                        if iTrace==1;
                            rootMask(pixelIdxList)=1;
                            if iframe~=1+frameShift
                                rootMask(pixelIdxList-allSize)=1;
                            end
                        else
                            if iTrace+iframe-1-frameShift<=size(demoImage,4)
                                rootMask(pixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),rootMask(:,:,iframe),rootType1Color);
end
clear rootMask;
end
function demoImage=makeType0TraceDemo(demoImage,bioTree,frameShift,xSize,ySize,allSize,traceType0Color)
disp('3.Make Type0 Trace demo');
traceMask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)+frameShift
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            centroidPixelIdxList=[];
            if bioTree{iframe}.root{iroot}.rootType==0
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    if iTrace==1;
                        if (iTrace+iframe-2-frameShift)<0
                            centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                            pixelIdxListinLine=pixelInLine(centroid2,centroid2);
                            centroidPixelIdxList=[centroidPixelIdxList,pixelIdxListinLine];
                            traceMask(centroidPixelIdxList)=1;
                        end
                        if (iTrace+iframe-2-frameShift)>=0
                            centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                            pixelIdxListinLine=pixelInLine(centroid2,centroid2);
                            centroidPixelIdxList=[centroidPixelIdxList+allSize,pixelIdxListinLine+(iTrace+iframe-2-frameShift)*allSize];
                            traceMask(centroidPixelIdxList)=1;
                        end
                    else
                        if iTrace+iframe-1-frameShift<=size(demoImage,4)
                            if (iTrace+iframe-2-frameShift)<0
                                centroid1=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}.Centroid;
                                centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                                pixelIdxListinLine=pixelInLine(centroid1,centroid2);
                                centroidPixelIdxList=[centroidPixelIdxList,pixelIdxListinLine];
                                traceMask(centroidPixelIdxList)=1;
                            end
                            if (iTrace+iframe-2-frameShift)>=0
                                centroid1=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}.Centroid;
                                centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                                pixelIdxListinLine=pixelInLine(centroid1,centroid2);
                                centroidPixelIdxList=[centroidPixelIdxList+allSize,pixelIdxListinLine+(iTrace+iframe-2-frameShift)*allSize];
                                traceMask(centroidPixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
% remain the trace
for iframe=2:size(demoImage,4)
    traceMask(:,:,iframe)=traceMask(:,:,iframe)|traceMask(:,:,iframe-1);
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),traceMask(:,:,iframe),traceType0Color);
end
clear traceMask;
end
function demoImage=makeType1TraceDemo(demoImage,bioTree,frameShift,xSize,ySize,allSize,traceType1Color)
disp('4.Make Type1 Trace demo');
traceMask=false(xSize,ySize,size(demoImage,4));
tempMask=false(xSize,ySize);
for iframe=1:size(demoImage,4)+frameShift
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            centroidPixelIdxList=[];
            if bioTree{iframe}.root{iroot}.rootType==1
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    if iTrace==1;
                        if (iTrace+iframe-2-frameShift)<0
                            centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                            pixelIdxListinLine=pixelInLine(centroid2,centroid2);
                            centroidPixelIdxList=[centroidPixelIdxList,pixelIdxListinLine];
                            tempMask(centroidPixelIdxList)=1;
                        end
                    else
                        if iTrace+iframe-1-frameShift<=size(demoImage,4)
                            if (iTrace+iframe-2-frameShift)<0
                                centroid1=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}.Centroid;
                                centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                                pixelIdxListinLine=pixelInLine(centroid1,centroid2);
                                tempMask(centroidPixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iNodeOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                centroidPixelIdxList=[];
                for iNodeTrace=1:size(bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                    if iNodeTrace==1
                        if (iNodeTrace+iframe-2-frameShift)<0
                            centroid2=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}.Centroid;
                            pixelIdxListinLine=pixelInLine(centroid2,centroid2);
                            centroidPixelIdxList=[centroidPixelIdxList,pixelIdxListinLine];
                            tempMask(centroidPixelIdxList)=1;
                        end
                    else
                        if iNodeTrace+iframe-1-frameShift<=size(demoImage,4)
                            if (iNodeTrace+iframe-2-frameShift)<0
                                centroid1=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace-1}.Centroid;
                                centroid2=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}.Centroid;
                                pixelIdxListinLine=pixelInLine(centroid1,centroid2);
                                centroidPixelIdxList=[centroidPixelIdxList,pixelIdxListinLine];
                                tempMask(centroidPixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
for iframe=1:size(demoImage,4)+frameShift
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            centroidPixelIdxList=[];
            if bioTree{iframe}.root{iroot}.rootType==1
                for iTrace=1:size(bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList,2)
                    if iTrace==1;
                        if (iTrace+iframe-2-frameShift)>=0
                            centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                            pixelIdxListinLine=pixelInLine(centroid2,centroid2);
                            centroidPixelIdxList=[centroidPixelIdxList+allSize,pixelIdxListinLine+(iTrace+iframe-2-frameShift)*allSize];
                            traceMask(centroidPixelIdxList)=1;
                        end
                    else
                        if iTrace+iframe-1-frameShift<=size(demoImage,4)
                            if (iTrace+iframe-2-frameShift)>=0
                                centroid1=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace-1}.Centroid;
                                centroid2=bioTree{iframe}.root{iroot}.traceInfo.measurment{iTrace}.Centroid;
                                pixelIdxListinLine=pixelInLine(centroid1,centroid2);
                                centroidPixelIdxList=[centroidPixelIdxList+allSize,pixelIdxListinLine+(iTrace+iframe-2-frameShift)*allSize];
                                traceMask(centroidPixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iNodeOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                centroidPixelIdxList=[];
                for iNodeTrace=1:size(bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                    if iNodeTrace==1
                        if (iNodeTrace+iframe-2-frameShift)>=0
                            centroid2=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}.Centroid;
                            pixelIdxListinLine=pixelInLine(centroid2,centroid2);
                            centroidPixelIdxList=[centroidPixelIdxList+allSize,pixelIdxListinLine+(iNodeTrace+iframe-2-frameShift)*allSize];
                            traceMask(centroidPixelIdxList)=1;
                        end
                    else
                        if iNodeTrace+iframe-1-frameShift<=size(demoImage,4)
                            if (iNodeTrace+iframe-2-frameShift)>=0
                                centroid1=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace-1}.Centroid;
                                centroid2=bioTree{iframe}.node{iNode}.Out{iNodeOut}.traceInfo.measurment{iNodeTrace}.Centroid;
                                pixelIdxListinLine=pixelInLine(centroid1,centroid2);
                                centroidPixelIdxList=[centroidPixelIdxList+allSize,pixelIdxListinLine+(iNodeTrace+iframe-2-frameShift)*allSize];
                                traceMask(centroidPixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
% remain the trace
traceMask(:,:,1)=traceMask(:,:,1)|tempMask;
for iframe=2:size(demoImage,4)
    traceMask(:,:,iframe)=traceMask(:,:,iframe)+traceMask(:,:,iframe-1);
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),traceMask(:,:,iframe),traceType1Color);
end
clear traceMask;
end
function demoImage=makeNodeType1Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,nodeType1Color)
disp('3.Make nodeType1 demo');
nodeType1Mask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)
    if ~isempty(bioTree{iframe+frameShift}.node)
        for iNode=1:size(bioTree{iframe+frameShift}.node,2)
            if bioTree{iframe+frameShift}.node{iNode}.type==1 || bioTree{iframe+frameShift}.node{iNode}.type==0 %remove nodetype=0 later
                for iNodeOut=1:size(bioTree{iframe+frameShift}.node{iNode}.Out,2)
                    for iNodeTrace=1:size(bioTree{iframe+frameShift}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                        if (iNodeTrace+iframe-2)>=0
                            if iNodeTrace+iframe-1<=size(demoImage,4)
                                pixelIdxList=bioTree{iframe+frameShift}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}+(iNodeTrace+iframe-2)*allSize;
                                nodeType1Mask(pixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),nodeType1Mask(:,:,iframe),nodeType1Color);
end
clear nodeType1Mask;
end
function demoImage=makeNodeType2Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,nodeType2Color)
disp('4.Make nodeType2 demo');
nodeType2Mask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)
    if ~isempty(bioTree{iframe+frameShift}.node)
        for iNode=1:size(bioTree{iframe+frameShift}.node,2)
            if bioTree{iframe+frameShift}.node{iNode}.type==2
                for iNodeOut=1:size(bioTree{iframe+frameShift}.node{iNode}.Out,2)
                    for iNodeTrace=1:size(bioTree{iframe+frameShift}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                        if (iNodeTrace+iframe-2)>=0
                            if iNodeTrace+iframe-1<=size(demoImage,4)
                                pixelIdxList=bioTree{iframe+frameShift}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}+(iNodeTrace+iframe-2)*allSize;
                                nodeType2Mask(pixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end

for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),nodeType2Mask(:,:,iframe),nodeType2Color);
end
clear nodeType2Mask;
end
function demoImage=makeNodeType3Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,nodeType3Color)
disp('5.Make nodeType0 demo');
nodeType0Mask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)
    if ~isempty(bioTree{iframe+frameShift}.node)
        for iNode=1:size(bioTree{iframe+frameShift}.node,2)
            if bioTree{iframe+frameShift}.node{iNode}.type==3
                for iNodeOut=1:size(bioTree{iframe+frameShift}.node{iNode}.Out,2)
                    for iNodeTrace=1:size(bioTree{iframe+frameShift}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList,2)
                        if (iNodeTrace+iframe-2)>=0
                            if iNodeTrace+iframe-1<=size(demoImage,4)
                                pixelIdxList=bioTree{iframe+frameShift}.node{iNode}.Out{iNodeOut}.traceInfo.pixelIdxList{iNodeTrace}+(iNodeTrace+iframe-2)*allSize;
                                nodeType0Mask(pixelIdxList)=1;
                            end
                        end
                    end
                end
            end
        end
    end
end

for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),nodeType0Mask(:,:,iframe),nodeType3Color);
end
clear nodeType0Mask;
end
function demoImage=makeLeafType0Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,leafType0Color)
disp('6.Make leaf demo');
leafMask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)
    if~isempty(bioTree{iframe+frameShift}.leavies)
        for iLeaf=1:size(bioTree{iframe+frameShift}.leavies,2)
            if bioTree{iframe+frameShift}.leavies{iLeaf}.leafType==0
                pixelIdxList=bioTree{iframe+frameShift}.leavies{iLeaf}.leaviesPixelDetail+(iframe-1)*allSize;
                leafMask(pixelIdxList)=1;
                if iframe~=size(demoImage,4)
                    leafMask(pixelIdxList+allSize)=1;
                end
            end
        end
    end
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),leafMask(:,:,iframe),leafType0Color);
end
clear leafMask;
end
function demoImage=makeLeafType1Demo(demoImage,bioTree,frameShift,xSize,ySize,allSize,leafType1Color)
disp('6.Make leaf demo');
leafMask=false(xSize,ySize,size(demoImage,4));
for iframe=1:size(demoImage,4)
    if~isempty(bioTree{iframe+frameShift}.leavies)
        for iLeaf=1:size(bioTree{iframe+frameShift}.leavies,2)
            if bioTree{iframe+frameShift}.leavies{iLeaf}.leafType==1
                pixelIdxList=bioTree{iframe+frameShift}.leavies{iLeaf}.leaviesPixelDetail+(iframe-1)*allSize;
                leafMask(pixelIdxList)=1;
                if iframe~=size(demoImage,4)
                    leafMask(pixelIdxList+allSize)=1;
                end
            end
        end
    end
end
for iframe=1:size(demoImage,4)
    demoImage(:,:,:,iframe)=imoverlay(demoImage(:,:,:,iframe),leafMask(:,:,iframe),leafType1Color);
end
clear leafMask;
end
function pixelIdxListinLine=pixelInLine(centroid1,centroid2)
xSize=1024;

if abs((centroid2(1)-centroid1(1)))>abs((centroid2(2)-centroid1(2)))
    if fix(centroid1(1))<=fix(centroid2(1))
        xinLine=fix(centroid1(1)):1:fix(centroid2(1));
    else
        xinLine=fix(centroid1(1)):-1:fix(centroid2(1));
    end
    if fix(centroid2(1)-centroid1(1))~=0
        k=((centroid2(2)-centroid1(2))/(centroid2(1)-centroid1(1)));
        b=centroid1(2)-centroid1(1)*k;
        yinLine=fix(xinLine.*k+b);
    else
        yinLine=fix(centroid2(2));
    end
    pixelIdxListinLine=xinLine.*xSize+yinLine;
    return;
end

if abs((centroid2(1)-centroid1(1)))<=abs((centroid2(2)-centroid1(2)))
    if fix(centroid1(2))<=fix(centroid2(2))
        yinLine=fix(centroid1(2)):1:fix(centroid2(2));
    else
        yinLine=fix(centroid1(2)):-1:fix(centroid2(2));
    end
    if fix(centroid2(1)-centroid1(1))~=0
        k=((centroid2(2)-centroid1(2))/(centroid2(1)-centroid1(1)));
        b=centroid1(2)-centroid1(1)*k;
        xinLine=fix((yinLine-b)./k);
    else
        xinLine=fix(centroid1(1));
    end
    pixelIdxListinLine=xinLine.*xSize+yinLine;
    return;
end
end