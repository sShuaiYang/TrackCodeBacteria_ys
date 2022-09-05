function [ultbioTree,frameNum,imageSize,dirFile]=bacteriaSimulationNew(rootNum,imageSize,frameNum)
dirFile=uigetdir();
dirImage=strcat(dirFile,'\tiff2matlab');   %put the folder of the save results
dirOriTreeSave=strcat(dirFile,'\bestBioTree');   %put the folder of the save results
mkdir(dirFile,'bioTreeResult');
mkdir(dirFile,'tiff2matlab')
mkdir(dirFile,'bestBioTree')
cd(dirImage);
clc;
% rootNum=30;
% imageSize=[1000,1000];
largeimageSize=imageSize+200;
image=false(imageSize);
largeimage=false(largeimageSize);
for iRoot=1:rootNum
    length=random('Norm',20,3);
    orientation=-90+180*rand;
    isProper=0;
    while isProper==0
        xPos=(imageSize(1)-40)*rand+20;
        yPos=(imageSize(2)-40)*rand+20;
        oriTree{iRoot}.info(1,:)=[length,orientation,xPos,yPos];
        if rand<0.5
            oriTree{iRoot}.info(1,5)=1;
        else
            oriTree{iRoot}.info(1,5)=0;
        end
        oriTree{iRoot}.isRoot=1;
        oriTree{iRoot}.beginFrame=1;
        oriTree{iRoot}.is2Node=[];
        oriTree{iRoot}.hasEnd=[];
        oriTree{iRoot}.realEnd=[];
        oriTree{iRoot}.climbouttime=[];
        oriTree{iRoot}.entering=[];
        oriTree{iRoot}.rootlength=oriTree{iRoot}.info(1,1);
        oriTree{iRoot}.isbiotree=1;
        oriTree{iRoot}.nodeInfo=[];
        [pixelIdxList,~,~,~]=getFatterIdx(oriTree{iRoot}.info,imageSize,largeimageSize);
        imageOne=false(imageSize);
        imageOne(pixelIdxList)=1;
        if max(max(imageOne&image))==0
            [oriTree{iRoot}.pixelIdxList{1,1},oriTree{iRoot}.largepixelIdxList{1,1},~,~]=getIdx(oriTree{iRoot}.info,imageSize,largeimageSize);
            isProper=1;
        end
    end
    largeimage(oriTree{iRoot}.largepixelIdxList{1,1})=1;
end
image=largeimage(101:largeimageSize(1)-100,101:largeimageSize(2)-100);
% largeimageStack(:,:,1)=largeimage;
imageStack(:,:,1)=image;
imageNum=1;
savefile=0;
for i=2:frameNum
    tic;
    disp(i)
    image=false(imageSize);
    largeimage=false(largeimageSize);
    for iRoot=1:size(oriTree,2)
        if isempty(oriTree{iRoot}.realEnd) && isempty(oriTree{iRoot}.is2Node)
                                                                                                    % tic;
            infoNum=size(oriTree{iRoot}.info,1);
            oriTree{iRoot}.info(infoNum+1,:)=bacteriaMovement(oriTree{iRoot}.info(infoNum,:),~isempty(oriTree{iRoot}.hasEnd),oriTree{iRoot}.info(1,5));
                                                                                                   %  Elapsed time is 0.000307 seconds.                                                                                            
            [IdxList,largeIdxList,linkBorder,Outsideornot]=getFatterIdx(oriTree{iRoot}.info(infoNum+1,:),imageSize,largeimageSize);
                                                                                                   %  Elapsed time is 0.004545 seconds.
            [oriTree{iRoot}.pixelIdxList{infoNum+1,1},oriTree{iRoot}.largepixelIdxList{infoNum+1,1},~,~]=getIdx(oriTree{iRoot}.info(infoNum+1,:),imageSize,largeimageSize);
                                                                                                    %  Elapsed time is 0.008671 seconds.
            if isempty(oriTree{iRoot}.hasEnd) && linkBorder
                oriTree{iRoot}.hasEnd=1;
                treeSize=size(oriTree,2);
                oriTree{iRoot}.climbouttime=i;
                if ~isempty(oriTree{iRoot}.isRoot)
                    if oriTree{iRoot}.climbouttime-oriTree{iRoot}.beginFrame<=1
                        oriTree{iRoot}.isbiotree=[];
                    end
                end
            end
                                                                                                   %Elapsed time is 0.008873 seconds.
            if ~isempty(oriTree{iRoot}.hasEnd) && linkBorder==0
                oriTree{iRoot}.realEnd=1;
                infoNew=oriTree{iRoot}.info(end-1,:);
                oriTree{iRoot}.info(end,:)=[];
                treeSize=size(oriTree,2);
                oriTree{treeSize+1}.info=infoNew;
                oriTree{treeSize+1}.beginFrame=i;
                oriTree{treeSize+1}.isRoot=1;
                oriTree{treeSize+1}.is2Node=[];
                oriTree{treeSize+1}.hasEnd=[];
                oriTree{treeSize+1}.realEnd=[];
                oriTree{treeSize+1}.climbouttime=[];
                oriTree{treeSize+1}.entering=[];
                [oriTree{treeSize+1}.pixelIdxList{1,1},oriTree{treeSize+1}.largepixelIdxList{1,1},~,~]=getIdx(oriTree{treeSize+1}.info,imageSize,largeimageSize);
                oriTree{treeSize+1}.rootlength=oriTree{iRoot}.rootlength;
                oriTree{treeSize+1}.isbiotree=1;
                oriTree{treeSize+1}.nodeInfo=[];
                largeimage(oriTree{treeSize+1}.largepixelIdxList{1,1})=1;
            end
                                                                                                    % Elapsed time is 0.008958 seconds.
            if ~isempty(oriTree{iRoot}.hasEnd) && Outsideornot==1
                oriTree{iRoot}.realEnd=1;
                infoNew=oriTree{iRoot}.info(end-1,:);
                infoNew=newproper(infoNew,imageSize);
                oriTree{iRoot}.info(end,:)=[];
                treeSize=size(oriTree,2);
                oriTree{treeSize+1}.info=infoNew;
                oriTree{treeSize+1}.beginFrame=i;
                oriTree{treeSize+1}.isRoot=1;
                oriTree{treeSize+1}.is2Node=[];
                oriTree{treeSize+1}.hasEnd=1;
                oriTree{treeSize+1}.realEnd=[];
                oriTree{treeSize+1}.climbouttime=[];
                oriTree{treeSize+1}.entering=1;
                [oriTree{treeSize+1}.pixelIdxList{1,1},oriTree{treeSize+1}.largepixelIdxList{1,1},~,~]=getIdx(oriTree{treeSize+1}.info,imageSize,largeimageSize);
                oriTree{treeSize+1}.rootlength=oriTree{iRoot}.rootlength;
                oriTree{treeSize+1}.isbiotree=1;
                oriTree{treeSize+1}.nodeInfo=[];
                largeimage(oriTree{treeSize+1}.largepixelIdxList{1,1})=1;
            end
                                                                                                       % Elapsed time is 0.009033 seconds.
            %             if ~isempty(oriTree{iRoot}.entering) && Outsideornot==1
            %                 oriTree{iRoot}.realEnd=1;
            %                 oriTree{iRoot}.pixelIdxList(end)=[];
            %                 oriTree{iRoot}.largepixelIdxList(end)=[];
            %                 treeSize=size(oriTree,2);
            %                 oriTree{treeSize+1}.info=oriTree{iRoot}.info;
            %                 oriTree{treeSize+1}.isRoot=1;
            %                 oriTree{treeSize+1}.beginFrame=i;
            %                 oriTree{treeSize+1}.is2Node=[];
            %                 oriTree{treeSize+1}.hasEnd=[];
            %                 oriTree{treeSize+1}.realEnd=[];
            %                 oriTree{treeSize+1}.climbouttime=[];
            %                 oriTree{treeSize+1}.entering=[];
            %                 [oriTree{treeSize+1}.pixelIdxList{1,1},oriTree{treeSize+1}.largepixelIdxList{1,1},~,~]=getIdx(oriTree{treeSize+1}.info,imageSize,largeimageSize);
            %                 oriTree{iRoot}.info(end,:)=[];
            %                 largeimage(oriTree{treeSize+1}.largepixelIdxList{1,1})=1;
            %             end
            largeimageOne=false(largeimageSize);
            largeimageOne(largeIdxList)=1;
            largeimageLeft=createOneLeaveImage(oriTree,iRoot,largeimageSize);
            isProper=properMovement(largeimageOne,largeimageLeft);
                                                                                                          % Elapsed time is 0.010746 seconds.
            if isProper==0
                oriTree{iRoot}.info(end,3:5)=oriTree{iRoot}.info(end-1,3:5);
                oriTree{iRoot}.info(end,1)=oriTree{iRoot}.info(end-1,1)+double(random('Norm',25,3))/900;                
                [oriTree{iRoot}.pixelIdxList{end,1},oriTree{iRoot}.largepixelIdxList{end,1},~,~]=getIdx(oriTree{iRoot}.info(end,:),imageSize,largeimageSize);
%                 oriTree{iRoot}.largepixelIdxList{end,1}=oriTree{iRoot}.largepixelIdxList{end-1,1};
            end
            if isempty(oriTree{iRoot}.realEnd) && isempty(oriTree{iRoot}.hasEnd)
                [needDivide,detachingCase]=bacteriaFate(oriTree{iRoot}.rootlength,oriTree{iRoot}.info(end,1),~isempty(oriTree{iRoot}.hasEnd) || isempty(oriTree{iRoot}.isbiotree));
                if detachingCase==1
                    oriTree{iRoot}.hasEnd=1;
                    oriTree{iRoot}.realEnd=1;
                    image(oriTree{iRoot}.largepixelIdxList{infoNum+1,1})=1;
                    continue
                end
                if needDivide==1
                    oriTree{iRoot}.is2Node=1;
                    oriTree{iRoot}.nodeInfo=[size(oriTree,2)+1,size(oriTree,2)+2];
                    oriTree{iRoot}.pixelIdxList(end)=[];
                    oriTree{iRoot}.largepixelIdxList(end)=[];
                    afterDivideCase=divisionSimulation(oriTree{iRoot}.info(end,:),imageSize,largeimageSize);
                    oriTree{iRoot}.info(end,:)=[];
                    for divideNum=1:2
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,:)=afterDivideCase(divideNum,:);
                        if rand<0.5
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,5)=1;
                        else
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,5)=0;
                        end
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.beginFrame=oriTree{iRoot}.beginFrame+size(oriTree{iRoot}.info,1);
                        [oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.pixelIdxList{1,1},oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.largepixelIdxList{1,1},divlinkBorder,divOutsideornot]=getIdxAfterDivision(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,:),imageSize,largeimageSize);
                        oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.nodeInfo=[];  
                        divOutsideornot=0; % new add by jzy
                        divlinkBorder=0;
                        if divOutsideornot==0
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.isRoot=[];
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.is2Node=[];
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.hasEnd=[];
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.realEnd=[];
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.climbouttime=[];
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.entering=[];
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.rootlength=oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.info(1,1);
                            if divlinkBorder==1
                                oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.hasEnd=1;
                                oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.climbouttime=i;
                                oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.isbiotree=[];
                            end
                        end
                        if divOutsideornot==1
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.hasEnd=1;
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.realEnd=1;
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.climbouttime=i;
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.isbiotree=[];
                        else
                            largeimage(oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.largepixelIdxList{1,1})=1;
                            oriTree{oriTree{iRoot}.nodeInfo(divideNum)}.isbiotree=1;
                        end
                    end
                    continue
                end
                                                                                             %  Elapsed time is 0.010958 seconds.
                largeimage(oriTree{iRoot}.largepixelIdxList{infoNum+1,1})=1;
            end
        end
    end
    for attachingNum=1:3
        if rand<0.001
            [oriTree,image,largeimage]=generateAttachingCase(oriTree,image,largeimage,i,imageSize,largeimageSize);
        end
    end
    imageNum=imageNum+1;
    image=largeimage(101:largeimageSize(1)-100,101:largeimageSize(2)-100);
    imageStack(:,:,imageNum)=image;
%     largeimageStack(:,:,i)=largeimage;
                                                                                         %  Elapsed time is 0.102240 seconds.for 10 cells
    if imageNum==1000
        %         parfor j=1:imageNum
        %             imageStack(:,:,j)=imfill(imageStack(:,:,j),'holes');
        %             imageStack(:,:,j)=bwmorph(imageStack(:,:,j),'remove');
        %         end
        savefile=savefile+1;
        save(strcat(num2str(savefile),'.mat'),'imageStack');
        ultbioTree=makebiotree(oriTree);
        bioTree1=generateBioTree(ultbioTree,i,imageSize,oriTree);
        save(strcat(dirOriTreeSave,'\bioTree1_',num2str(savefile)),'bioTree1');
        clear imageStack
        clear bioTree1
%         clear largeimageStack
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
function [pixelIdxList,largepixelIdxList,linkBorder,Outsideornot]=getIdx(info,imageSize,largeimageSize)
model=generateBacteria(info(1),info(2));
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
cc.PixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
if any(cc.PixelList(:,1)<=1)==1 || any(cc.PixelList(:,1)>=imageSize(2))==1 || any(cc.PixelList(:,2)<=1)==1 || any(cc.PixelList(:,2)>=imageSize(1))==1
    linkBorder=1;
else
    linkBorder=0;
end
if all(((cc.PixelList(:,1)<=1)==1) + ((cc.PixelList(:,1)>=imageSize(2))==1) + ((cc.PixelList(:,2)<=1)==1) + ((cc.PixelList(:,2)>=imageSize(1))==1))
    Outsideornot=1;
else
    Outsideornot=0;
end
pixelIdxList=cc.PixelList(:,2)+(cc.PixelList(:,1)-1)*imageSize(1);
largepixelIdxList=100+cc.PixelList(:,2)+(100+(cc.PixelList(:,1)-1))*largeimageSize(1);
end
function [pixelIdxList,largepixelIdxList,linkBorder,Outsideornot]=getIdxAfterDivision(info,imageSize,largeimageSize)
model=generateBacteria(info(1),info(2));
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
cc.PixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
if any(cc.PixelList(:,1)<=1)==1 || any(cc.PixelList(:,1)>=imageSize(2))==1 || any(cc.PixelList(:,2)<=1)==1 || any(cc.PixelList(:,2)>=imageSize(1))==1
    linkBorder=1;
else
    linkBorder=0;
end
if all(((cc.PixelList(:,1)<=1)==1) + ((cc.PixelList(:,1)>=imageSize(2))==1) + ((cc.PixelList(:,2)<=1)==1) + ((cc.PixelList(:,2)>=imageSize(1))==1))
    Outsideornot=1;
else
    Outsideornot=0;
end
pixelIdxList=cc.PixelList(:,2)+(cc.PixelList(:,1)-1)*imageSize(1);
largepixelIdxList=100+cc.PixelList(:,2)+(100+(cc.PixelList(:,1)-1))*largeimageSize(1);
end
function [pixelIdxList,largepixelIdxList,linkBorder,Outsideornot]=getFatterIdx(info,imageSize,largeimageSize)
model=generateFatterBacteria(info(1),info(2));                           %   Elapsed time is 0.002855 seconds.
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
cc.PixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
if any(cc.PixelList(:,1)<=1)==1 || any(cc.PixelList(:,1)>=imageSize(2))==1 || any(cc.PixelList(:,2)<=1)==1 || any(cc.PixelList(:,2)>=imageSize(1))==1
    linkBorder=1;
else
    linkBorder=0;
end
if all(((cc.PixelList(:,1)<=1)==1) + ((cc.PixelList(:,1)>=imageSize(2))==1) + ((cc.PixelList(:,2)<=1)==1) + ((cc.PixelList(:,2)>=imageSize(1))==1))
    Outsideornot=1;
else
    Outsideornot=0;
end
pixelIdxList=cc.PixelList(:,2)+(cc.PixelList(:,1)-1)*imageSize(1);
largepixelIdxList=100+cc.PixelList(:,2)+(100+(cc.PixelList(:,1)-1))*largeimageSize(1);
end
function model=generateBacteria(length,oritation)
% if length<=1
%    length=4;
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
function model=generateFatterBacteria(length,oritation)
% if length<=1
%    length=4;
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

model=imrotate(model,oritation);                %Elapsed time is 0.002908 seconds.

% model=bwmorph(model,'remove');
end
function afterDivideCase=divisionSimulation(info,imageSize,largeimageSize)
for i=1:2
    afterDivideCase(i,1)=info(1)/2;
end
afterDivideCase(1,2)=info(2);
afterDivideCase(2,2)=info(2)+180;
afterDivideCase(1,3:4)=[info(3)-(2+info(1)/4)*sin(info(2)/180*pi),info(4)+(2+info(1)/4)*cos(info(2)/180*pi)];
afterDivideCase(2,3:4)=[info(3)+(2+info(1)/4)*sin(info(2)/180*pi),info(4)-(2+info(1)/4)*cos(info(2)/180*pi)];
afterDivideCase=Divideadjust(afterDivideCase,imageSize,largeimageSize);
end
function adjustedDivideCase=Divideadjust(DivideCase,imageSize,largeimageSize)
index=2;
dividecase(2)=DivideCase(1,2);
dividecase(1)=2*DivideCase(1,1)+2;
dividecase(3:4)=(DivideCase(1,3:4)+DivideCase(1,3:4))/2;
model=generateFatterBacteria(dividecase(1),dividecase(2));  
cc=regionprops(model,'centroid','pixelList');
cc.PixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+dividecase(4));
cc.PixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+dividecase(3));
while any(cc.PixelList(:,1)<=1)==1 || any(cc.PixelList(:,1)>=imageSize(2))==1 || any(cc.PixelList(:,2)<=1)==1 || any(cc.PixelList(:,2)>=imageSize(1))==1
    if any(cc.PixelList(:,1)<=1)==1
        DivideCase(1:2,4)=DivideCase(1:2,4)+index;
        cc.PixelList(:,1)=cc.PixelList(:,1)+index;
    end
    if any(cc.PixelList(:,1)>=imageSize(2))==1
        DivideCase(1:2,4)=DivideCase(1:2,4)-index;
        cc.PixelList(:,1)=cc.PixelList(:,1)-index;
    end
    if any(cc.PixelList(:,2)<=1)==1
        DivideCase(1:2,3)=DivideCase(1:2,3)+index;
        cc.PixelList(:,2)=cc.PixelList(:,2)+index;
    end
    if any(cc.PixelList(:,2)>=imageSize(1))==1
        DivideCase(1:2,3)=DivideCase(1:2,3)-index;
        cc.PixelList(:,2)=cc.PixelList(:,2)-index;
    end
end
adjustedDivideCase=DivideCase;
end
function [oriTree,image,largeimage]=generateAttachingCase(oriTree,image,largeimage,i,imageSize,largeimageSize)
oriTreeNum=size(oriTree,2);
length=random('Norm',25,3);
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
if rand<0.5
    oriTree{oriTreeNum+1}.info(1,5)=1;
else
    oriTree{oriTreeNum+1}.info(1,5)=0;
end
oriTree{oriTreeNum+1}.beginFrame=i;
oriTree{oriTreeNum+1}.isRoot=1;
oriTree{oriTreeNum+1}.is2Node=[];
oriTree{oriTreeNum+1}.hasEnd=[];
oriTree{oriTreeNum+1}.realEnd=[];
oriTree{oriTreeNum+1}.climbouttime=[];
oriTree{oriTreeNum+1}.entering=[];
oriTree{oriTreeNum+1}.rootlength=oriTree{oriTreeNum+1}.info(1,1);
oriTree{oriTreeNum+1}.nodeInfo=[];
[pixelIdxList,largepixelIdxList,~,~]=getFatterIdx(oriTree{oriTreeNum+1}.info,imageSize,largeimageSize);
largeimageOne=false(largeimageSize);
largeimageOne(largepixelIdxList)=1;
%     if max(max(imageOne&image))==0
[oriTree{oriTreeNum+1}.pixelIdxList{1,1},oriTree{oriTreeNum+1}.largepixelIdxList{1,1},~,~]=getIdx(oriTree{oriTreeNum+1}.info,imageSize,largeimageSize);
oriTree{oriTreeNum+1}.isbiotree=1;
%         isProper=1;
%     end
% end
largeimage(oriTree{oriTreeNum+1}.largepixelIdxList{1,1})=1;
image=largeimage(101:largeimageSize(1)-100,101:largeimageSize(2)-100);
end
function newInfo=bacteriaMovement(info,outsideimage,stayOrMove)
if stayOrMove==0
    newInfo(1,1)=info(1);
    if outsideimage==0
%         newInfo(1,1)=newInfo(1,1)+double(random('Norm',25,3))/900;
          newInfo(1,1)=newInfo(1,1)+double(25/random('Unif',600,900));
    end
    if rand<0.001
        newInfo(1,2)=info(2)-183+6*rand;
        newInfo(1,3:4)=info(3:4);
        newInfo(1,5)=0;
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
    newInfo(1,5)=0;
else
    newInfo=info;
%     newInfo(1,1)=newInfo(1,1)+double(random('Norm',25,3))/900;
    newInfo(1,1)=newInfo(1,1)+double(25/random('Unif',600,900));
    newInfo(1,2)=info(2)-183+6*rand;
end
end
function [needDivide,detachingCase]=bacteriaFate(beginLength,recentLength,outsideimage)
mustDivide=2.3;
unableDivide=1.7;
if beginLength<20
    beginLength=25;
end
if beginLength>=30
    beginLength=30;
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
if outsideimage==1
    needDivide=0;
end
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
function imageLeft=createOneLeaveImage(oriTree,iNum,largeimageSize)
info=oriTree{iNum}.info(end,:);
focusRadius=2*info(1);
imageLeft=zeros(largeimageSize);
for iRoot=1:size(oriTree,2)
    if isempty(oriTree{iRoot}.hasEnd) && isempty(oriTree{iRoot}.is2Node)
        if iRoot<iNum
            info2=oriTree{iRoot}.info(end-1,3:4);
            distance=(sum((info2-info(3:4)).^2))^0.5;
            if distance<=focusRadius
                imageLeft(oriTree{iRoot}.largepixelIdxList{end-1})=1;
            end
        end
        if iRoot>iNum
            info2=oriTree{iRoot}.info(end,3:4);
            distance=(sum((info2-info(3:4)).^2))^0.5;
            if distance<=focusRadius
                imageLeft(oriTree{iRoot}.largepixelIdxList{end})=1;
            end
        end
    end
end
end
function properinfo=newproper(info,imageSize)
info(2)=info(2);
info(3:4)=imageSize-info(3:4);
model=generateBacteria(info(1),info(2));
cc=regionprops(model,'centroid','pixelList');
newpixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
newpixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
while all(((newpixelList(:,1)<=1)==1) + ((newpixelList(:,1)>=imageSize(2))==1) + ((newpixelList(:,2)<=1)==1) + ((newpixelList(:,2)>=imageSize(1))==1))
    info=positioncontract(info,imageSize);
    newpixelList(:,1)=ceil(cc.PixelList(:,1)-cc.Centroid(1)+info(4));
    newpixelList(:,2)=ceil(cc.PixelList(:,2)-cc.Centroid(2)+info(3));
end
properinfo=info;
end
function newinfo=positioncontract(info,imageSize)
newinfo=info;
if info(3)>=0.5*imageSize(1)
    newinfo(3)=info(3)-1;
end
if info(3)<0.5*imageSize(1)
    newinfo(3)=info(3)+1;
end
if info(4)>=0.5*imageSize(2)
    newinfo(4)=info(4)-1;
end
if info(4)<0.5*imageSize(2)
    newinfo(4)=info(4)+1;
end
end
function ultbioTree=makebiotree(oribioTree)
bioTreesize=size(oribioTree,2);
numtrans=zeros(1,bioTreesize);
treebranchnum=0;
for iRoot=1:size(oribioTree,2)
    if ~isempty(oribioTree{iRoot}.isbiotree)
        treebranchnum=treebranchnum+1;
        ultbioTree{treebranchnum}.isRoot=oribioTree{iRoot}.isRoot;
        ultbioTree{treebranchnum}.beginFrame=oribioTree{iRoot}.beginFrame;
        ultbioTree{treebranchnum}.is2Node=oribioTree{iRoot}.is2Node;       
        ultbioTree{treebranchnum}.hasEnd=oribioTree{iRoot}.hasEnd;
        ultbioTree{treebranchnum}.climbouttime=oribioTree{iRoot}.climbouttime;
        ultbioTree{treebranchnum}.rootlength=oribioTree{iRoot}.rootlength;
        ultbioTree{treebranchnum}.nodeInfo=oribioTree{iRoot}.nodeInfo;
        ultbioTree{treebranchnum}.pixelIdxList=[];
        numtrans(iRoot)=treebranchnum;
        if ~isempty(oribioTree{iRoot}.climbouttime)
            existime=oribioTree{iRoot}.climbouttime-oribioTree{iRoot}.beginFrame;
            for iframe=1:existime
                ultbioTree{treebranchnum}.pixelIdxList{iframe,1}=oribioTree{iRoot}.pixelIdxList{iframe,1};
                ultbioTree{treebranchnum}.info(iframe,:)=oribioTree{iRoot}.info(iframe,:);
            end            
        else
            ultbioTree{treebranchnum}.pixelIdxList=oribioTree{iRoot}.pixelIdxList;
        end            
    end                
end
for iRoot=1:size(ultbioTree,2)
    if ~isempty(ultbioTree{iRoot}.nodeInfo)        
        ultbioTree{iRoot}.nodeInfo=[numtrans(ultbioTree{iRoot}.nodeInfo(1)) numtrans(ultbioTree{iRoot}.nodeInfo(2))];
    end    
end
end