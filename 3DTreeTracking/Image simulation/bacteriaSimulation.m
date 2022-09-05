function [oriTree,frameNum,imageSize,dirFile]=bacteriaSimulation(frameNum)
dirFile=uigetdir();
dirImage=strcat(dirFile,'\tiff2matlab');   %put the folder of the save results
dirOriTreeSave=strcat(dirFile,'\bestBioTree');   %put the folder of the save results
mkdir(dirFile,'bioTreeResult');
mkdir(dirFile,'tiff2matlab')
mkdir(dirFile,'bestBioTree')
cd(dirImage);
clc;
rootNum=30;
imageSize=[1050,1000];
image=false(imageSize);
for iRoot=1:rootNum
   length=random('Poisson',20);
   orientation=-90+180*rand;
   isProper=0;
   while isProper==0
       xPos=(imageSize(1)-40)*rand+20;
       yPos=(imageSize(2)-40)*rand+20;
       oriTree{iRoot}.info(1,:)=[length,orientation,xPos,yPos];
       oriTree{iRoot}.isRoot=1;
       oriTree{iRoot}.beginFrame=1;
       oriTree{iRoot}.is2Node=[];
       oriTree{iRoot}.hasEnd=[];  
       [pixelIdxList,~]=getFatterIdx(oriTree{iRoot}.info,imageSize);
       imageOne=false(imageSize);
       imageOne(pixelIdxList)=1;
       if max(max(imageOne&image))==0
           [oriTree{iRoot}.pixelIdxList{1,1},~]=getIdx(oriTree{iRoot}.info,imageSize);
           isProper=1;
       end
   end
   image(oriTree{iRoot}.pixelIdxList{1,1})=1;
end
imageStack(:,:,1)=image;
imageNum=1;
savefile=0;
for i=2:frameNum
    disp(i)
    image=false(imageSize);
    for iRoot=1:size(oriTree,2)
        if isempty(oriTree{iRoot}.hasEnd) && isempty(oriTree{iRoot}.is2Node)
            infoNum=size(oriTree{iRoot}.info,1);
            oriTree{iRoot}.info(infoNum+1,:)=bacteriaMovement(oriTree{iRoot}.info(infoNum,:));
            [IdxList,linkBorder]=getFatterIdx(oriTree{iRoot}.info(infoNum+1,:),imageSize);
            [oriTree{iRoot}.pixelIdxList{infoNum+1,1},~]=getIdx(oriTree{iRoot}.info(infoNum+1,:),imageSize);
            if linkBorder==1
                oriTree{iRoot}.hasEnd=1;
                infoNew=oriTree{iRoot}.info(end-1,:);
                infoNew(2)=infoNew(2)+180;
                infoNew(3:4)=imageSize-infoNew(3:4);
                oriTree{iRoot}.info(end,:)=[];
                treeSize=size(oriTree,2);
                oriTree{iRoot}.pixelIdxList(end)=[];
                oriTree{treeSize+1}.info=infoNew;
                oriTree{treeSize+1}.beginFrame=i;
                oriTree{treeSize+1}.isRoot=1;
                oriTree{treeSize+1}.is2Node=[];
                oriTree{treeSize+1}.hasEnd=[];
                [oriTree{treeSize+1}.pixelIdxList{1,1},~]=getIdx(oriTree{treeSize+1}.info,imageSize);
                image(oriTree{treeSize+1}.pixelIdxList{1,1})=1;
%                 oriTree{iRoot}.info(end,:)=oriTree{iRoot}.info(end-1,:);
%                 oriTree{iRoot}.info(end,2)=oriTree{iRoot}.info(end-1,2)+180;
%                 oriTree{iRoot}.pixelIdxList(end)=oriTree{iRoot}.pixelIdxList(end-1);
%                 image(oriTree{iRoot}.pixelIdxList{infoNum+1,1})=1;
                %                 [oriTree,image]=generateAttachingCase(oriTree,image,i,imageSize);
                continue
            end
            imageOne=false(imageSize);            
            imageOne(IdxList)=1;
            imageLeft=createOneLeaveImage(oriTree,iRoot,imageSize);
            isProper=properMovement(imageOne,imageLeft);
            if isProper==0
                oriTree{iRoot}.info(end,:)=oriTree{iRoot}.info(end-1,:);
                oriTree{iRoot}.pixelIdxList{end,1}=oriTree{iRoot}.pixelIdxList{end-1,1};
            end
            [needDivide,detachingCase]=bacteriaFate(oriTree{iRoot}.info(1,1),oriTree{iRoot}.info(end,1));
            if detachingCase==1
                oriTree{iRoot}.hasEnd=1;
                image(oriTree{iRoot}.pixelIdxList{infoNum+1,1})=1;
                continue
            end
            if needDivide==1
                oriTree{iRoot}.is2Node=1;
                oriTree{iRoot}.nodeInfo=[size(oriTree,2)+1,size(oriTree,2)+2];
                oriTree{iRoot}.info(end,:)=[];
                oriTree{iRoot}.pixelIdxList(end)=[];
                afterDivideCase=divisionSimulation(oriTree{iRoot}.info(end,:));
                for divideNum=1:2
                    oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,:)=afterDivideCase(divideNum,:);
                    oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.beginFrame=oriTree{iRoot}.beginFrame+size(oriTree{iRoot}.info,1);
                    oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.isRoot=[];
                    oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.is2Node=[];
                    oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.hasEnd=[];
                    [oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.pixelIdxList{1,1},~]=getIdx(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,:),imageSize);
                    image(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.pixelIdxList{1,1})=1;
                end
                continue
            end
            image(oriTree{iRoot}.pixelIdxList{infoNum+1,1})=1;
        end
    end
    for attachingNum=1:3
        if rand<0.001
           [oriTree,image]=generateAttachingCase(oriTree,image,i,imageSize); 
        end
    end
    imageNum=imageNum+1;
    imageStack(:,:,imageNum)=image;
    if imageNum==1000
%         parfor j=1:imageNum
%             imageStack(:,:,j)=imfill(imageStack(:,:,j),'holes');
%             imageStack(:,:,j)=bwmorph(imageStack(:,:,j),'remove');
%         end
        savefile=savefile+1;
        save(strcat(num2str(savefile),'.mat'),'imageStack');
        bioTree1=generateBioTree(oriTree,i,imageSize);
        save(strcat(dirOriTreeSave,'\bioTree1_',num2str(savefile)),'bioTree1');
        clear imageStack
        clear bioTree1
        imageStack=false(0);
        imageNum=0;
    end
end
if ~isempty(imageStack)
    savefile=savefile+1;
%     parfor j=1:imageNum
%         imageStack(:,:,j)=imfill(imageStack(:,:,j),'holes');
%         imageStack(:,:,j)=bwmorph(imageStack(:,:,j),'remove');
%     end
    save(strcat(num2str(savefile),'.mat'),'imageStack');
    save(strcat(dirOriTreeSave,'\bioTree1_',num2str(savefile)),'bioTree1');
    clear imageStack
end
end
function [oriTree,image]=generateAttachingCase(oriTree,image,i,imageSize)
oriTreeNum=size(oriTree,2);
length=random('Poisson',28);
orientation=-90+180*rand;
% isProper=0;
image2=imdilate(image,ones(5));
image2(1:40,:)=1;
image2(end-39:end,:)=1;
image2(:,1:40)=1;
image2(:,end-39:end)=1;
restImage=1-image2;
cc=regionprops(restImage,'PixelList');
% while isProper==0
numPix=size(cc.PixelList,1);
numPix=fix(rand*numPix);
xPos=cc.PixelList(numPix,2);
yPos=cc.PixelList(numPix,1);
%     xPos=(imageSize(1)-40)*rand+20;
%     yPos=(imageSize(2)-40)*rand+20;
oriTree{oriTreeNum+1}.info(1,:)=[length,orientation,xPos,yPos];
oriTree{oriTreeNum+1}.beginFrame=i;
oriTree{oriTreeNum+1}.isRoot=1;
oriTree{oriTreeNum+1}.is2Node=[];
oriTree{oriTreeNum+1}.hasEnd=[];
[pixelIdxList,~]=getFatterIdx(oriTree{oriTreeNum+1}.info,imageSize);
imageOne=false(imageSize);
imageOne(pixelIdxList)=1;
%     if max(max(imageOne&image))==0
[oriTree{oriTreeNum+1}.pixelIdxList{1,1},~]=getIdx(oriTree{oriTreeNum+1}.info,imageSize);
%         isProper=1;
%     end
% end
image(oriTree{oriTreeNum+1}.pixelIdxList{1,1})=1;
end
function [pixelIdxList,linkBorder]=getIdx(info,imageSize)
model=generateBacteria(info(1),info(2));
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
cc.PixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
if any(cc.PixelList(:,1)<=1)==1 || any(cc.PixelList(:,1)>=imageSize(2))==1 || any(cc.PixelList(:,2)<=1)==1 || any(cc.PixelList(:,2)>=imageSize(1))==1
    linkBorder=1;
else
    linkBorder=0;
end
pixelIdxList=cc.PixelList(:,2)+(cc.PixelList(:,1)-1)*imageSize(1);
end
function [pixelIdxList,linkBorder]=getFatterIdx(info,imageSize)
model=generateFatterBacteria(info(1),info(2));
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
cc.PixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
if any(cc.PixelList(:,1)<=1)==1 || any(cc.PixelList(:,1)>=imageSize(2))==1 || any(cc.PixelList(:,2)<=1)==1 || any(cc.PixelList(:,2)>=imageSize(1))==1
    linkBorder=1;
else
    linkBorder=0;
end
pixelIdxList=cc.PixelList(:,2)+(cc.PixelList(:,1)-1)*imageSize(1);
end
function model=generateFatterBacteria(length,oritation)
% if length<=1
%     p=1
% end
length=ceil(length+2);  
oritation=ceil(oritation);
model=false(11,length);
model(4:8,1)=true;
model(4:8,length)=true;
model(3:9,2)=true;
model(3:9,length-1)=true;
model(2:10,3)=true;
model(2:10,length-2)=true;
model(:,4:length-3)=true;
model=imrotate(model,oritation);
% model=bwmorph(model,'remove');
end
function [needDivide,detachingCase]=bacteriaFate(beginLength,recentLength)
mustDivide=2.2;
unableDivide=1.7;
if beginLength<25
    beginLength=25;
end
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
function newInfo=bacteriaMovement(info)
newInfo(1,1)=info(1)+double(random('Poisson',25))/900;
if rand<0.001
    newInfo(1,2)=info(2)-183+6*rand;
    newInfo(1,3:4)=info(3:4);
    return
end
vProbability=rand;
lamda=2;
velocity=-log(1-vProbability)/lamda;
if velocity<=0.1
    orientationThr=-150*velocity +25;
else
    if velocity<=0.5
        orientationThr=-12*velocity+10;
    else
        orientationThr=3;
    end
end
newInfo(1,2)=info(2)-orientationThr+orientationThr*2*rand;
newInfo(1,3:4)=[info(3),info(4)]+[-1*sin(newInfo(2)/180*pi),1*cos(newInfo(2)/180*pi)]*velocity;
end
function afterDivideCase=divisionSimulation(info)
for i=1:2
    afterDivideCase(i,1)=info(1)/2;
end
afterDivideCase(1,2)=info(2);
afterDivideCase(2,2)=info(2)+180;
afterDivideCase(1,3:4)=[info(3)-(2+info(1)/4)*sin(info(2)/180*pi),info(4)+(2+info(1)/4)*cos(info(2)/180*pi)];
afterDivideCase(2,3:4)=[info(3)+(2+info(1)/4)*sin(info(2)/180*pi),info(4)-(2+info(1)/4)*cos(info(2)/180*pi)];
end
function isProper=properMovement(imageOne,imageAll)
sigma=0.1;
% if ~isequal(size(imageOne),size(imageAll))
%     p=1;
% end
overlapLimit=25;
overlapArea=imageOne&imageAll;
ratio=max(double(numel(overlapArea(overlapArea==1)))/overlapLimit,double(numel(overlapArea(overlapArea==1)))/numel(imageOne(imageOne==1)));
if rand<exp(-ratio/sigma)
    isProper=1;
else
    isProper=0;
end
end
function imageLeft=createOneLeaveImage(oriTree,iNum,imageSize)
info=oriTree{iNum}.info(end,:);
focusRadius=2*info(1);
imageLeft=zeros(imageSize);
for iRoot=1:size(oriTree,2)
    if isempty(oriTree{iRoot}.hasEnd) && isempty(oriTree{iRoot}.is2Node)
        if iRoot<iNum
            info2=oriTree{iRoot}.info(end-1,3:4);
            distance=(sum((info2-info(3:4)).^2))^0.5;
            if distance<=focusRadius
                imageLeft(oriTree{iRoot}.pixelIdxList{end-1})=1;
            end
        end
        if iRoot>iNum
            info2=oriTree{iRoot}.info(end,3:4);
            distance=(sum((info2-info(3:4)).^2))^0.5;
            if distance<=focusRadius
                imageLeft(oriTree{iRoot}.pixelIdxList{end})=1;
            end
        end
    end
end
end
function model=generateBacteria(length,oritation)
% if length<=1
%     p=1
% end
length=ceil(length);
oritation=ceil(oritation);
model=false(7,length);
model(3:5,1)=true;
model(3:5,length)=true;
model(2:6,2)=true;
model(2:6,length-1)=true;
model(:,3:length-2)=true;
model=imrotate(model,oritation);
% model=bwmorph(model,'remove');
end