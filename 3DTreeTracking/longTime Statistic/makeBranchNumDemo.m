function makeBranchNumDemo(bioTree,branchList,startFrame,stepFrame,endFrame,dirFile)
warning off all
dirImages=strcat(dirFile,'\tiff2matlab');
dirSave=strcat(dirFile,'\branchNumDemo1');
mkdir(dirSave);
cd(dirImages);
nameList=dir(dirImages);
vars = whos('-file', nameList(3).name);
branchList(branchList(:,3)==0,:)=[];
coreBranchNum=size(branchList,1);
stackSize=vars.size(3);
stackNum=fix(startFrame/(stackSize));
if stackNum~=startFrame/(stackSize)
    fileName=nameList(stackNum+3).name;
else
    fileName=nameList(stackNum+2).name;
end
imageStack=loadImageStack(fileName);
for iframe=startFrame:stepFrame:endFrame;
    disp(iframe);
    stackNumNext=fix(iframe/(stackSize));
    if stackNumNext==iframe/(stackSize)
        stackNumNext=stackNumNext-1;
    end
    if stackNumNext~=stackNum
        clear imageStack;
        fileName=nameList(stackNumNext+3).name;
        imageStack=loadImageStack(fileName);
        stackNum=stackNumNext;
    end
    demoImage=imageStack(:,:,iframe-stackNum*stackSize);
    figure('visible','off');imshow(demoImage)
    bacteriaFrameInfo=getEachBacteriaInFrame(bioTree,iframe);
    for i=1:size(bacteriaFrameInfo.centroidInfo,1)
        if bacteriaFrameInfo.bacteriaInfo(i)<=coreBranchNum
            text(bacteriaFrameInfo.centroidInfo(i,1),bacteriaFrameInfo.centroidInfo(i,2),num2str(bacteriaFrameInfo.bacteriaInfo(i)),'Color','y','FontSize',6)
        end
    end
    saveFile2=strcat(dirSave,'\',num2str(iframe),'.tif');
    saveas(gcf,saveFile2);
    close all
end
end
function bacteriaFrameInfo=getEachBacteriaInFrame(bioTree,iframe)
bacteriaFrameInfo.centroidInfo=[];
bacteriaFrameInfo.bacteriaInfo=[];
for smallFrame=1:iframe
    if ~isempty(bioTree{smallFrame}.root)
        if ~isempty(bioTree{smallFrame}.root)
            for iRoot=1:size(bioTree{smallFrame}.root,2)
                traceInfo=bioTree{smallFrame}.root{iRoot}.traceInfo.pixelIdxList;
                if isempty(traceInfo)
                    continue
                end
                if size(traceInfo,2)>=iframe+1-smallFrame
                    branchIndex=bioTree{smallFrame}.root{iRoot}.branchIndex;
                    bacteriaFrameInfo.centroidInfo=[bacteriaFrameInfo.centroidInfo;bioTree{smallFrame}.root{iRoot}.traceInfo.measurment{iframe+1-smallFrame}(1).Centroid];
                    bacteriaFrameInfo.bacteriaInfo=[bacteriaFrameInfo.bacteriaInfo;bioTree{smallFrame}.root{iRoot}.branchIndex];
                end
            end
        end
    end
    if ~isempty(bioTree{smallFrame}.node)
        for iNode=1:size(bioTree{smallFrame}.node,2)
            for iOut=1:size(bioTree{smallFrame}.node{iNode}.Out,2);
                traceInfo=bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList;
                if size(traceInfo,2)>=iframe+1-smallFrame
                    branchIndex=bioTree{smallFrame}.node{iNode}.branchIndex;
                    bacteriaFrameInfo.centroidInfo=[bacteriaFrameInfo.centroidInfo;bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.measurment{iframe+1-smallFrame}(1).Centroid];
                    bacteriaFrameInfo.bacteriaInfo=[bacteriaFrameInfo.bacteriaInfo;bioTree{smallFrame}.node{iNode}.branchIndex];
                end
            end
        end
    end
end
end
function imageStack=loadImageStack(fileName)
tempImage=load(fileName);
imageStack=tempImage.imageStack;
end