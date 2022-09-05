function bioTree=bioTreeSimulation()
rootNum=30;
frameNum=4000;
imageSize=[1666,1593];
parfor iFrame=1:frameNum
    bioTree{iFrame}.root=[];
    bioTree{iFrame}.node=[];
    bioTree{iFrame}.leavies=[];
end
bioTree{1}.imageSize=imageSize;
for iRoot=1:rootNum
    disp(iRoot)
   length=random('Poisson',25);
   orientation=-90+180*rand;
   xPos=(imageSize(1)-20)*rand;
   yPos=(imageSize(2)-20)*rand;
   bioTree{1}.root{iRoot}.modelInfo{1}.length=length;
   bioTree{1}.root{iRoot}.modelInfo{1}.orientation=orientation;
   bioTree{1}.root{iRoot}.modelInfo{1}.xyPos=[xPos,yPos];
   [pixelIdxList,~]=getIdx(bioTree{1}.root{iRoot}.modelInfo{1},imageSize);
   bioTree{1}.root{iRoot}.rootPixelDetail=pixelIdxList;
   bioTree{1}.root{iRoot}.traceInfo.pixelIdxList{1}=pixelIdxList;
   needDivide=0;
   detachingCase=0;
   i=1;
   while needDivide==0 && detachingCase==0
       i=i+1;
       bioTree{1}.root{iRoot}.modelInfo{1,i}=bacteriaMovement(bioTree{1}.root{iRoot}.modelInfo{i-1});
       [bioTree{1}.root{iRoot}.traceInfo.pixelIdxList{1,i},linkBorder]=getIdx(bioTree{1}.root{iRoot}.modelInfo{1,i},imageSize);
       [needDivide,detachingCase]=bacteriaFate(bioTree{1}.root{iRoot}.modelInfo{1}.length,bioTree{1}.root{iRoot}.modelInfo{i}.length);
       if linkBorder==1 || detachingCase==1
           break
       end
       if needDivide==1
           nodeNum=size(bioTree{i+1}.node,2);
           bioTree{i+1}.node{nodeNum+1}.In.isNode=0;
           bioTree{i+1}.node{nodeNum+1}.In.rootInfo=[1,iRoot];
           bioTree{i+1}.node{nodeNum+1}.In.nodeInfo=[];
           afterDivideCase=divisionSimulation(bioTree{1}.root{iRoot}.modelInfo{end});
           for iOut=1:2
           bioTree{i+1}.node{nodeNum+1}.Out.modelInfo{1}=afterDivideCase{iOut};
           [bioTree{i+1}.node{nodeNum+1}.Out.pixelIdxList{1},linkBorder]=getIdx(bioTree{i+1}.node{nodeNum+1}.Out.modelInfo{1},imageSize);
           end
           bioTree{1}.root{iRoot}.is2Node=1;
           bioTree{1}.root{iRoot}.nodeInfo=[i+1,nodeNum];
           bioTree{1}.root{iRoot}.leafInfo=[];  
       end
   end       
end
end
function [pixelIdxList,linkBorder]=getIdx(modelInfo,imageSize)
model=generateBacteria(modelInfo.length,modelInfo.orientation);
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=cc.PixelList(1)-cc.Centroid(2)+ceil(modelInfo.xyPos(1));
cc.PixelList(:,2)=cc.PixelList(2)-cc.Centroid(1)+ceil(modelInfo.xyPos(2));
if any(cc.PixelList(:,1)<=0 | cc.PixelList(:,1)>imageSize(1))==1
    linkBorder=1;
else
    linkBorder=0;
end
pixelIdxList=cc.PixelList(:,1)+(cc.PixelList(:,2)-1)*imageSize(1);
end
function newModel=bacteriaMovement(model)
newModel.length=model.length+double(random('Poisson',25))/900;
newModel.orientation=model.orientation-30+60*rand;
newModel.xyPos=model.xyPos+[0.1,0.1]*rand;
end
function [needDivide,detachingCase]=bacteriaFate(beginLength,recentLength)
mustDivide=2.2;
unableDivide=1.7;
if rand<(recentLength/beginLength-unableDivide)/(mustDivide-unableDivide)
    needDivide=1;
else
    needDivide=0;
end
if rand<0.00001
    detachingCase=1;
else
    detachingCase=0;
end
end
function afterDivideCase=divisionSimulation(modelInfo)
for i=1:2
    afterDivideCase{i}.orientation=modelInfo.orientation;
    afterDivideCase{i}.length=modelInfo.length/2-2;
end
afterDivideCase{1}.xyPos=[modelInfo.xyPos(1)+modelInfo.length/4*cos(modelInfo.orientation),modelInfo.xyPos(2)-modelInfo.length/4*sin(modelInfo.orientation)];
afterDivideCase{2}.xyPos=[modelInfo.xyPos(1)-modelInfo.length/4*cos(modelInfo.orientation),modelInfo.xyPos(2)+modelInfo.length/4*sin(modelInfo.orientation)];
end