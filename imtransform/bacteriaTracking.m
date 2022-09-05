function [ CCbranch , bioTree] = bacteriaTracking( preImage , maskImage )
%BACTERIATRACKING Summary of this function goes here
%   Detailed explanation goes here
connectMask = false(size(maskImage,1),size(maskImage,2),2);
connectMask(:,:,1) = preImage;
connectMask(:,:,2) = maskImage;
CCbranch = bwconncomp(connectMask(:,:,:),26);
bioTree = twoFrameConnect(CCbranch,1);
bioTree=bioTreeMeasure(bioTree,0,CCbranch.ImageSize(1),CCbranch.ImageSize(2));

end

function bioTree=twoFrameConnect(CC,startFrame) %initial conenct i and i+1 frame as well as intial the data structure
xSize=CC.ImageSize(1);
ySize=CC.ImageSize(2);
allSize=xSize*ySize;
endFrame=startFrame+1;
bioTree{startFrame}.leavies=[];
bioTree{startFrame}.node=[];
bioTree{startFrame}.root=[];
bioTree{endFrame}.root=[];
bioTree{endFrame}.node=[];
bioTree{endFrame}.leavies=[];
count1=1;count2=1;count3=1;count4=1;
for iObject=1:CC.NumObjects
    
    pixelFrame=fix(CC.PixelIdxList{iObject}./allSize);
    startFramePixel=CC.PixelIdxList{iObject}(pixelFrame==0);
    endFramePixel=CC.PixelIdxList{iObject}(pixelFrame==1);
    
    if isempty(startFramePixel)
        bioTree{1,endFrame}.root{count1}.rootPixelDetail=endFramePixel-allSize;
        bioTree{1,endFrame}.root{count1}.rootPixelNum=numel(endFramePixel-allSize);
        bioTree{1,endFrame}.root{count1}.is2Node=false;
        bioTree{1,endFrame}.root{count1}.leafInfo=[endFrame,count4];
        bioTree{1,endFrame}.root{count1}.nodeInfo=[];
        bioTree{1,endFrame}.root{count1}.traceInfo.pixelIdxList{1}=endFramePixel-allSize;
        bioTree{1,endFrame}.leavies{count4}.is2Node=false;
        bioTree{1,endFrame}.leavies{count4}.rootInfo=[endFrame,count1];
        bioTree{1,endFrame}.leavies{count4}.nodeInfo=[];
        bioTree{1,endFrame}.leavies{count4}.leaviesPixelDetail=endFramePixel-allSize;
         bioTree{1,endFrame}.leavies{count4}.leaviesPixelNum=numel( endFramePixel-allSize);
        count1=count1+1;
        count4=count4+1;
        continue;
    end
    
    if isempty(endFramePixel)
        bioTree{1,startFrame}.root{count3}.is2Node=false;
        bioTree{1,startFrame}.root{count3}.leafInfo=[startFrame,count2];
        bioTree{1,startFrame}.root{count3}.nodeInfo=[];
        bioTree{1,startFrame}.root{count3}.rootPixelDetail= startFramePixel;
        bioTree{1,startFrame}.root{count3}.rootPixelNum= numel( startFramePixel);
        bioTree{1,startFrame}.root{count3}.traceInfo.pixelIdxList{1}=startFramePixel;
        bioTree{1,startFrame}.leavies{count2}.is2Node=false;
        bioTree{1,startFrame}.leavies{count2}.nodeInfo=[];
        bioTree{1,startFrame}.leavies{count2}.rootInfo=[startFrame,count3];
        bioTree{1,startFrame}.leavies{count2}.leaviesPixelDetail=startFramePixel;
         bioTree{1,startFrame}.leavies{count2}.leaviesPixelNum=numel( startFramePixel);
        count2=count2+1;
        count3=count3+1;
        continue;
    end
    
    if (~isempty(startFramePixel))&&(~isempty(endFramePixel))
        bioTree{1,startFrame}.root{count3}.is2Node=false;
        bioTree{1,startFrame}.root{count3}.leafInfo=[endFrame,count4];
        bioTree{1,startFrame}.root{count3}.nodeInfo=[];
        bioTree{1,startFrame}.root{count3}.rootPixelDetail=startFramePixel;
        bioTree{1,startFrame}.root{count3}.rootPixelNum=numel(startFramePixel);
        bioTree{1,startFrame}.root{count3}.traceInfo.pixelIdxList{1}=startFramePixel;
        bioTree{1,startFrame}.root{count3}.traceInfo.pixelIdxList{2}=endFramePixel-allSize;
        bioTree{1,endFrame}.leavies{count4}.is2Node=false;
        bioTree{1,endFrame}.leavies{count4}.rootInfo=[startFrame,count3];
        bioTree{1,endFrame}.leavies{count4}.nodeInfo=[];
        bioTree{1,endFrame}.leavies{count4}.leaviesPixelDetail=endFramePixel-allSize;
        bioTree{1,endFrame}.leavies{count4}.leaviesPixelNum=numel(endFramePixel-allSize);
        count3=count3+1;
        count4=count4+1;
        continue;
    end
    
end
end

