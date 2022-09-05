function bioTree=stochasticMovingSimulation(rootNum,imageSize,frameNum,pesistanceLength,dirFile)
% dirFile=uigetdir();
dirImage=strcat(dirFile,'\tiff2matlab');   %put the folder of the save results
dirOriTreeSave=strcat(dirFile,'\bestBioTree');   %put the folder of the save results
mkdir(dirFile,'bioTreeResult');
mkdir(dirFile,'tiff2matlab')
mkdir(dirFile,'bestBioTree')
cd(dirImage);
clc;

stayProbability=0;
circleModel=createCircle();
imageEdgeXy=getImageEdge(imageSize);
% initialize point on the first frame
image=false(imageSize);
for iRoot=1:rootNum
    length=random('Poisson',30);
    orientation=360*rand;
    xPos=(imageSize(1)-40)*rand+20;
    yPos=(imageSize(2)-40)*rand+20;
    oriTree{iRoot}.info(1,:)=[length,orientation,xPos,yPos];
    if rand<=stayProbability
        oriTree{iRoot}.info(1,5)=0;
    else
        oriTree{iRoot}.info(1,5)=1;
    end
    oriTree{iRoot}.isRoot=1;
    oriTree{iRoot}.beginFrame=1;
    oriTree{iRoot}.is2Node=[];
    oriTree{iRoot}.hasEnd=[];
    oriTree{iRoot}.nodeInfo=[];
    oriTree{iRoot}.pixelIdxList{1,1}=(round(oriTree{iRoot}.info(1,4))-1)*imageSize(1)+round(oriTree{iRoot}.info(1,3));
    [picNew,~]=getIdx(oriTree{iRoot}.info(1,:),imageSize,circleModel);
    image=image | picNew;
end
image=imdilate(image,circleModel);
imageNum=1;
saveFile=0;
saveFile1=0;
imageStack(:,:,imageNum)=image;
% generate oriTree
for i=2:frameNum
    disp(i)
    image=false(imageSize);
    for iRoot=1:size(oriTree,2)
        if isempty(oriTree{iRoot}.hasEnd) && isempty(oriTree{iRoot}.is2Node)
            infoNum=size(oriTree{iRoot}.info,1);
            oriTree{iRoot}.info(infoNum+1,:)=bacteriaMovement(oriTree{iRoot}.info(infoNum,:),oriTree{iRoot}.info(1,5),pesistanceLength);
            oriTree{iRoot}.pixelIdxList{infoNum+1,1}=(round(oriTree{iRoot}.info(infoNum+1,4))-1)*imageSize(1)+round(oriTree{iRoot}.info(infoNum+1,3));
            [picNew,linkBorder]=getIdx(oriTree{iRoot}.info(infoNum+1,:),imageSize,circleModel);
            if linkBorder==1
                oriTree{iRoot}.hasEnd=1;
                [oriTree,image]=generateAttachingCase(oriTree,image,i,imageSize,stayProbability,circleModel,oriTree{iRoot}.info(end,1),imageEdgeXy);
            end
            if linkBorder==0
                [needDivide,detachingCase]=bacteriaFate(oriTree{iRoot}.info(end,1));
                if detachingCase==1
                    oriTree{iRoot}.hasEnd=1;
                    [picNew,~]=getIdx(oriTree{iRoot}.info(end,:),imageSize,circleModel);
                    image=image|picNew;
                    continue
                end
                if needDivide==1
                    oriTree{iRoot}.is2Node=1;
                    oriTree{iRoot}.nodeInfo=[size(oriTree,2)+1,size(oriTree,2)+2];
                    oriTree{iRoot}.pixelIdxList(end)=[];
                    afterDivideCase=divisionSimulation(oriTree{iRoot}.info(end,:));
                    oriTree{iRoot}.info(end,:)=[];
                    for divideNum=1:2
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,:)=afterDivideCase(divideNum,:);
                        if rand<stayProbability
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,5)=0;
                        else
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,5)=1;
                        end
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.beginFrame=oriTree{iRoot}.beginFrame+size(oriTree{iRoot}.info,1);
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.pixelIdxList{1,1}=(round(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,4))-1)*imageSize(1)+round(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,3));
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.nodeInfo=[];
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.isRoot=[];
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.is2Node=[];
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.hasEnd=[];
                        [picNew,~]=getIdx(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info,imageSize,circleModel);
                        image=image|picNew;
                    end
                    continue
                end
                image=image|picNew;
            end
        end
    end
    for attachingNum=1:3
        if rand<0.0002
            [oriTree,image]=generateAttachingCase(oriTree,image,i,imageSize,stayProbability,circleModel,random('Poisson',30));
        end
    end
    imageNum=imageNum+1;
    image=imdilate(image,circleModel);
    imageStack(:,:,imageNum)=image;
    if mod(imageNum,200)==0
        saveFile=saveFile+1;
        imageStack=imageStack(:,:,end-200+1:end);
        save(strcat(num2str(saveFile),'.mat'),'imageStack');
        clear imageStack
        imageStack=false(0);
    end
    if imageNum==1000 || i==frameNum
        saveFile1=saveFile1+1;
        bioTree=generateBioTree(oriTree,i,imageSize);
        save(strcat(dirOriTreeSave,'\bioTree_',num2str(saveFile1)),'bioTree');
        if i~=frameNum
            clear bioTree
        end
        %         clear largeimageStack
        imageNum=0;
    end
end
save(strcat(dirFile,'\oriTree')','oriTree');
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
bacteriaFrameInfo=getEachBacteriaInFrame(bioTree);
save(strcat(dirFile,'\bacteriaFrameInfo')','bacteriaFrameInfo');
end
function newInfo=bacteriaMovement(info,stayOrMove,pesistanceLength)
% % for bacteria situation
% if stayOrMove==0
%     newInfo(1,1)=info(1);
%     if rand<0.001
%         newInfo(1,2)=info(2)-183+6*rand;
%         newInfo(1,3:4)=info(3:4);
%         newInfo(1,5)=0;
%         return
%     end
%     vProbability=rand;
%     lamda=2;
%     velocity=-log(1-vProbability)/lamda;
%     if velocity<=0.1
%         orientationThr=-150*velocity +25;
%     else
%         if velocity<=0.5
%             orientationThr=-12*velocity+10;
%         else
%             orientationThr=3;
%         end
%     end
%     newInfo(1,2)=info(2)-orientationThr+orientationThr*2*rand;
%     newInfo(1,3:4)=[info(3),info(4)]+[-1*sin(newInfo(2)/180*pi),1*cos(newInfo(2)/180*pi)]*velocity;
%     newInfo(1,5)=0;
% else
%     newInfo=info;
%     newInfo(1,1)=newInfo(1,1)+double(random('Norm',25,3))/900;
%     newInfo(1,2)=info(2)-183+6*rand;
% end

% % for particle situation
% % for pesistanceLength
% vProbability=rand;
% lamda=0.2;
% % velocity=-log(1-vProbability)/lamda;
% velocity=3;
% % orientationRange=80;   % [0,180]
% % persistanceLength=[];
% newSita=random('norm',exp(-3/pesistanceLength),(1-exp(-3/pesistanceLength))/10);
% while newSita>1 || newSita<0
%     newSita=random('norm',exp(-3/pesistanceLength),(1-exp(-3/pesistanceLength))/10);
% end
% newSita=acos(newSita)/pi*180;
% if rand<0.5
%     newSita=-newSita;
% end
% newInfo=info;
% newInfo(1,1)=newInfo(1,1)+double(random('Norm',30,3))/900;
% newInfo(1,2)=info(2)+newSita;
% newInfo(1,3:4)=[info(3),info(4)]+[-1*sin(newInfo(2)/180*pi),1*cos(newInfo(2)/180*pi)]*velocity;
% newInfo(1,5)=0;

% for stochasitc moving 2015_1_19
% for pesistanceLength
vProbability=rand;
velocity=20*rand;
newSita=rand*360;
newInfo=info;
newInfo(1,1)=newInfo(1,1)+double(random('Norm',30,3))/900;
newInfo(1,2)=info(2)+newSita;
newInfo(1,3:4)=[info(3),info(4)]+[-1*sin(newInfo(2)/180*pi),1*cos(newInfo(2)/180*pi)]*velocity;
newInfo(1,5)=0;
end
function [picNew,linkBorder]=getIdx(info,imageSize,circleModel)
model=circleModel;  %   Elapsed time is 0.002855 seconds.
info(3:4)=round(info(3:4));
if info(4)<=1 || info(4)>=imageSize(2) || info(3)<=1 || info(3)>=imageSize(1)
    linkBorder=1;
    picNew=false(imageSize);
    return
else
    linkBorder=0;
end
picture=false(imageSize);
picture(info(3),info(4))=1;
% picture=imdilate(picture,model);
% picPix=find(picture==1);
picNew=picture;
end
function circleModel=createCircle()
se=strel('disk',10,0);
circleModel=zeros(21,21);
circleModel(ceil(end/2),ceil(end/2))=1;
circleModel=imdilate(circleModel,se);
circleModel=logical(imerode(circleModel,ones(3)));
end
function imageEdgeXy=getImageEdge(imageSize)
imageEdge=logical(ones(imageSize));
imageEdge1=bwmorph(imageEdge,'remove');
imageEdge=imageEdge-imageEdge1;
imageEdge=bwmorph(imageEdge,'remove');
[imageEdgeXy(:,1),imageEdgeXy(:,2)]=find(imageEdge==1);
end
function [oriTree,image]=generateAttachingCase(oriTree,image,i,imageSize,stayProbability,circleModel,length,imageEdgeXy)
oriTreeNum=size(oriTree,2);
orientation=360*rand;
if nargin==7
    xPos=(imageSize(1)-40)*rand+20;
    yPos=(imageSize(2)-40)*rand+20;
else
    if nargin==8
        possibleNum=size(imageEdgeXy,1);
        num=ceil(possibleNum*rand);
        xPos=imageEdgeXy(num,1);
        yPos=imageEdgeXy(num,2);
    end
end
oriTree{oriTreeNum+1}.info(1,:)=[length,orientation,xPos,yPos];
if rand<=stayProbability
    oriTree{oriTreeNum+1}.info(1,5)=0;
else
    oriTree{oriTreeNum+1}.info(1,5)=1;
end
oriTree{oriTreeNum+1}.beginFrame=i;
oriTree{oriTreeNum+1}.isRoot=1;
oriTree{oriTreeNum+1}.is2Node=[];
oriTree{oriTreeNum+1}.hasEnd=[];
oriTree{oriTreeNum+1}.nodeInfo=[];
oriTree{oriTreeNum+1}.pixelIdxList{1,1}=(round(oriTree{oriTreeNum+1}.info(1,4))-1)*imageSize(1)+round(oriTree{oriTreeNum+1}.info(1,3));
[picNew,~]=getIdx(oriTree{oriTreeNum+1}.info,imageSize,circleModel);
image=image|picNew;
end
function [needDivide,detachingCase]=bacteriaFate(recentLength)
mustDivide=2.05;
unableDivide=1.95;
if rand<(recentLength/35-unableDivide)/(mustDivide-unableDivide)
    needDivide=1;
else
    needDivide=0;
end 
if rand<0.00005
    detachingCase=1;
else
    detachingCase=0;
end
end
function afterDivideCase=divisionSimulation(info)
for i=1:2
    afterDivideCase(i,1)=random('Poisson',30);
end
afterDivideCase(1,2)=info(2);
afterDivideCase(2,2)=info(2)+180;
afterDivideCase(1,3:4)=[info(3),info(4)];
afterDivideCase(2,3:4)=[info(3),info(4)];
end
function bacteriaFrameInfo=getEachBacteriaInFrame(bioTree)
for iframe=1:size(bioTree,2)
    bacteriaFrameInfo{iframe}.centroidInfo=[];
    bacteriaFrameInfo{iframe}.lengthInfo=[];
    bacteriaFrameInfo{iframe}.bacteriaInfo=[];
end
for iframe=1:size(bioTree,2)
    disp(iframe)
    if ~isempty(bioTree{iframe}.root)
        if ~isempty(bioTree{iframe}.root)
            for iRoot=1:size(bioTree{iframe}.root,2)
                traceInfo=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList;
                if isempty(traceInfo)
                    continue
                end
                for iTrace=1:size(traceInfo,2)
%                     yCo=ceil((bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{iTrace}-1)/2160);
%                     xCo=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{iTrace}-2160*(yCo-1);
                    bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;bioTree{iframe}.root{iRoot}.traceInfo.traceCentroid(iTrace,:)];
                    bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo=[bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo;[iframe,iRoot,0,iTrace,bioTree{iframe}.root{iRoot}.branchIndex]];
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2);
                traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList;
                for iTrace=1:size(traceInfo,2)
%                     yCo=ceil((bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}-1)/2160);
%                     xCo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}-2160*(yCo-1);
                    bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.traceCentroid(iTrace,:)];
                    bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo=[bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo;[iframe,iNode,iOut,iTrace,bioTree{iframe}.node{iNode}.branchIndex]];
                end
            end
        end
    end
end
end
function bioTree=generateBioTree(oriTree,frameNum,imageSize)
for iframe=1:frameNum
    bioTree{iframe}.root=[];
    bioTree{iframe}.node=[];
    bioTree{iframe}.leavies=[];
end
bioTree{1}.imageSize=imageSize;
for iNum=1:size(oriTree,2)
    if ~isempty(oriTree{iNum}.isRoot)
        iframe=oriTree{iNum}.beginFrame;
        iRoot=size(bioTree{iframe}.root,2)+1;
        bioTree{iframe}.root{iRoot}.rootPixelDetail=oriTree{iNum}.pixelIdxList{1};
        bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList=oriTree{iNum}.pixelIdxList';
        bioTree{iframe}.root{iRoot}.traceInfo.traceCentroid=oriTree{iNum}.info(:,3:4);
        if ~isempty(oriTree{iNum}.is2Node)
            bioTree{iframe}.root{iRoot}.is2Node=1;
            bioTree{iframe}.root{iRoot}.leafInfo=[];
            nodeFrame=size(oriTree{iNum}.pixelIdxList,1)+iframe;
            bioTree{iframe}.root{iRoot}.nodeInfo=[nodeFrame,size(bioTree{nodeFrame}.node,2)+1,1];
            nodeInfo=bioTree{iframe}.root{iRoot}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode=0;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo=[iframe,iRoot];
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo=[];
            for i=1:2
                oriTree{oriTree{iNum}.nodeInfo(i)}.bioTreeNode=[nodeInfo(1),nodeInfo(2),i];
            end
        end
        if isempty(oriTree{iNum}.is2Node)
            bioTree{iframe}.root{iRoot}.is2Node=0;
            bioTree{iframe}.root{iRoot}.nodeInfo=[];
            if ~isempty(oriTree{iNum}.hasEnd)
                leafFrame=size(oriTree{iNum}.pixelIdxList,1)+iframe-1;
            else
                leafFrame=frameNum;
            end
            bioTree{iframe}.root{iRoot}.leafInfo=[leafFrame,size(bioTree{leafFrame}.leavies,2)+1];
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=0;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iRoot];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=oriTree{iNum}.pixelIdxList{end};
        end
    end
    if isempty(oriTree{iNum}.isRoot)
        bioTreeNode=oriTree{iNum}.bioTreeNode;
        bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.traceInfo.pixelIdxList=oriTree{iNum}.pixelIdxList';
        bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.traceInfo.traceCentroid=oriTree{iNum}.info(:,3:4);
        if ~isempty(oriTree{iNum}.is2Node)
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.is2Node=1;
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.leafInfo=[];
            nodeFrame=size(oriTree{iNum}.pixelIdxList,1)+bioTreeNode(1);
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.nodeInfo=[nodeFrame,size(bioTree{nodeFrame}.node,2)+1,1];
            nodeInfo=bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode=1;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.rootInfo=[];
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo=bioTreeNode;
            for i=1:2
                oriTree{oriTree{iNum}.nodeInfo(i)}.bioTreeNode=[nodeInfo(1),nodeInfo(2),i];
            end
        end
        if isempty(oriTree{iNum}.is2Node)
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.is2Node=0;
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.nodeInfo=[];
            if ~isempty(oriTree{iNum}.hasEnd)
                leafFrame=size(oriTree{iNum}.pixelIdxList,1)+bioTreeNode(1)-1;
            else
                leafFrame=frameNum;
            end
            bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.leafInfo=[leafFrame,size(bioTree{leafFrame}.leavies,2)+1];
            leafInfo=bioTree{bioTreeNode(1)}.node{bioTreeNode(2)}.Out{bioTreeNode(3)}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=1;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=bioTreeNode;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=oriTree{iNum}.pixelIdxList{end};
        end
    end
end
end