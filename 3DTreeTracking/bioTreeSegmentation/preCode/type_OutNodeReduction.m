function bioTree=type_OutNodeReduction(bioTree)
for iframe=1:size(bioTree,2)
    % % for iframe=1:156
%     if iframe==1732
%         p=4210;
%     end
%     disp(strcat(num2str(iframe),'_'))
% fprintf('\b\n')
dispFrame(iframe)
if ~isempty(bioTree{iframe}.node)
    for iNode=1:size(bioTree{iframe}.node,2)
        outNum=size(bioTree{iframe}.node{iNode}.Out,2);
        for iOut=1:outNum
            bioTree=type_OutNodeTracking(bioTree,[iframe,iNode,iOut]);
        end
    end
end
end
fprintf('\n')
end
function dispFrame(iConnect) 
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
function bioTree=type_OutNodeTracking(bioTree,nodeInfo)%further track the N input and one Output node, return N new complete bacterial trace
usefulNode=0;
imageSize=bioTree{1}.imageSize;
[xyMin,BWImage]=idx2Xy(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{1},imageSize);
CC=bwconncomp(BWImage);
regionNum=CC.NumObjects;
if regionNum<2
    return
else
    for i=1:regionNum
        image=false(size(BWImage));
        image(CC.PixelIdxList{i})=1;
        pixelIdxListIn{i}=xy2Idx(xyMin,image,imageSize);
        newTrace{i}.traceInfo.pixelIdxList{1}=pixelIdxListIn{i};
        newTrace{i}.is2Node=[];
        newTrace{i}.isBreakNode=[];
    end
end
if size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,2)>=2
    for iTrace=2:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,2)
        [pixelIdxListNew,canDivideorNot]=basicDivideSolutionforOut(pixelIdxListIn,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList{iTrace},imageSize);
        if canDivideorNot==1
            leafNum=size(bioTree{iTrace+nodeInfo(1)-2}.leavies,2);
            [newTrace,pixelIdxListIn,isBreak,error]=getNewTrace(newTrace,pixelIdxListNew,iTrace+nodeInfo(1)-2,leafNum);
            if error==1
                usefulNode=1;
                newTrace=[];
                return
            end
            if isBreak==true
                for iIn=1:size(newTrace,2)
                    if  isempty(newTrace{iIn}.is2Node)
                        if iTrace<size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,2)
                            newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(iTrace+1:end)];
                        end
                    end
                end
                break;
            end
        else
            if canDivideorNot==0
                if iTrace==1
                    usefulNode=1;
                    return
                else
                    for iIn=1:size(newTrace,2)
                        if isempty(newTrace{iIn}.is2Node)
                            newTrace{iIn}.isBreakNode=1;
                            newTrace{iIn}.breakNodeInfo=[nodeInfo(1)+iTrace-1,0];
                        end
                    end
                    break;
                end
            end
        end
    end
end
testis2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node;
if testis2Node==true
    testnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo;
    inNum=size(bioTree{testnodeInfo(1)}.node{testnodeInfo(2)}.In,2);
end
realInNode=0;
realInLeaf=0;
for iIn=1:size(newTrace,2)
    if ~isempty(newTrace{iIn}.isBreakNode)
        if isempty(newTrace{iIn}.is2Node)
            newTrace{iIn}.is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node;
            newTrace{iIn}.afterBreakTrace.pixelIdxList=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(iTrace:end);
            if newTrace{iIn}.is2Node==1
                newTrace{iIn}.nodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo;
            else
                if newTrace{iIn}.is2Node==0
                    newTrace{iIn}.leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo;
                end
            end
        end
    else
        if isempty(newTrace{iIn}.is2Node)
            newTrace{iIn}.is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node;
            if newTrace{iIn}.is2Node==true
                newTrace{iIn}.nodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo;
                realInNode=realInNode+1;
                if realInNode~=1
                    newTrace{iIn}.nodeInfo(1,3)=inNum+realInNode-1;
                end
            end
            if newTrace{iIn}.is2Node==false
                newTrace{iIn}.leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo;
                if realInLeaf~=0
                    newTrace{iIn}.leafInfo(1,2)=size(bioTree{newTrace{iIn}.leafInfo(1)}.leavies,2)+realInLeaf;
                end
                realInLeaf=realInLeaf+1;
            end
        end
    end
end
beforeNodeNum=0;
newNodeIn=0;
num2Leaf=0;
for i=1:size(newTrace,2)
    beforeNodeNum=beforeNodeNum+1;    
    if isempty(newTrace{i}.isBreakNode)
        if newTrace{i}.is2Node==1
            if beforeNodeNum==1
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=newTrace{i}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=newTrace{i}.traceInfo.pixelIdxList;
                nextNode=newTrace{i}.nodeInfo;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=nodeInfo;
                bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
            else
                if beforeNodeNum>=2
                    breakNodeInNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.is2Node=1;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.nodeInfo=newTrace{i}.nodeInfo;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.leafInfo=[];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.traceInfo.pixelIdxList=newTrace{i}.traceInfo.pixelIdxList;
                    nextNode=newTrace{i}.nodeInfo;
                    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
                    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),breakNodeInNum+1];
                    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
                end
            end
        end
        if newTrace{i}.is2Node==0;
            if beforeNodeNum==1
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=0;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=newTrace{i}.leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=newTrace{i}.traceInfo.pixelIdxList;
                nextLeaf=newTrace{i}.leafInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=nodeInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.leaviesPixelDetail=newTrace{i}.traceInfo.pixelIdxList{end};
            else
                if beforeNodeNum>=2
                    breakNodeInNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.is2Node=0;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.leafInfo=newTrace{i}.leafInfo;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.nodeInfo=[];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.traceInfo.pixelIdxList=newTrace{i}.traceInfo.pixelIdxList;
                    nextLeaf=newTrace{i}.leafInfo;
                    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
                    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[nodeInfo(1),nodeInfo(2),breakNodeInNum+1];
                    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
                    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.leaviesPixelDetail=newTrace{i}.traceInfo.pixelIdxList{end};
                end
            end
        end
    else
        if newTrace{i}.isBreakNode==1
            newNodeIn=newNodeIn+1;
            if newNodeIn==1
                newNode=newTrace{i}.breakNodeInfo;
                nodeNum=size(bioTree{newNode(1)}.node,2)+1;
                bioTree{newNode(1)}.node{nodeNum}.Out{1}.is2Node=newTrace{i}.is2Node;
                if bioTree{newNode(1)}.node{nodeNum}.Out{1}.is2Node==1
                    bioTree{newNode(1)}.node{nodeNum}.Out{1}.nodeInfo=newTrace{i}.nodeInfo;
                    bioTree{newNode(1)}.node{nodeNum}.Out{1}.leafInfo=[];
                    bioTree{newNode(1)}.node{nodeNum}.Out{1}.traceInfo.pixelIdxList=newTrace{i}.afterBreakTrace.pixelIdxList;
                    nextNode=newTrace{i}.nodeInfo;
                    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
                    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[newNode(1),nodeNum,1];
                    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
                else if bioTree{newNode(1)}.node{nodeNum}.Out{1}.is2Node==0
                        bioTree{newNode(1)}.node{nodeNum}.Out{1}.leafInfo=newTrace{i}.leafInfo;
                        bioTree{newNode(1)}.node{nodeNum}.Out{1}.nodeInfo=[];
                        bioTree{newNode(1)}.node{nodeNum}.Out{1}.traceInfo.pixelIdxList=newTrace{i}.afterBreakTrace.pixelIdxList;
                        nextLeaf=newTrace{i}.leafInfo;
                        bioTree{nextLeaf(1)}.leaf{nextLeaf(2)}.is2Node=1;
                        bioTree{nextLeaf(1)}.leaf{nextLeaf(2)}.nodeInfo=[newNode(1),nodeNum,1];
                        bioTree{nextLeaf(1)}.leaf{nextLeaf(2)}.rootInfo=[];
                    end
                end
            end
            if beforeNodeNum==1
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[newNode(1),nodeNum,newNodeIn];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[];
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=newTrace{i}.traceInfo.pixelIdxList;
                bioTree{newNode(1)}.node{nodeNum}.In{newNodeIn}.isNode=1;
                bioTree{newNode(1)}.node{nodeNum}.In{newNodeIn}.rootInfo=[];
                bioTree{newNode(1)}.node{nodeNum}.In{newNodeIn}.nodeInfo=nodeInfo;
            else
                if beforeNodeNum>=2
                    breakNodeInNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.is2Node=1;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.nodeInfo=[newNode(1),nodeNum,newNodeIn];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.leafInfo=[];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{breakNodeInNum+1}.traceInfo.pixelIdxList=newTrace{i}.traceInfo.pixelIdxList;
                    bioTree{newNode(1)}.node{nodeNum}.In{newNodeIn}.isNode=1;
                    bioTree{newNode(1)}.node{nodeNum}.In{newNodeIn}.rootInfo=[];
                    bioTree{newNode(1)}.node{nodeNum}.In{newNodeIn}.nodeInfo=[nodeInfo(1),nodeInfo(2),breakNodeInNum+1];
                end
            end
        end
    end
end
end
function [newTrace,pixelIdxListIn,isBreak,error]=getNewTrace(newTrace,pixelIdxListNew,iframe,leafNum)
emptyCount=0;
nleaf=0;
for iIn=1:size(newTrace,2)
    if ~isempty(pixelIdxListNew{iIn});
        if ~isempty(newTrace{iIn}.is2Node)
            error=1;
            isBreak=[];
            pixelIdxListIn=[];
            return
        end
        newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,pixelIdxListNew(iIn)];
    else
        if isempty(newTrace{iIn}.is2Node)
        nleaf=nleaf+1;
        newTrace{iIn}.is2Node=false;
        newTrace{iIn}.leafInfo=[iframe,leafNum+nleaf];
        end
        emptyCount= emptyCount+1;
    end
end
if size(newTrace,2)-emptyCount==1
    isBreak=true;
else
    isBreak=false;
end
error=0;
pixelIdxListIn=pixelIdxListNew;
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=pixelIdxList-(yresult-1)*xSize;
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function [afterDivide,canDivideorNot]=basicDivideSolutionforOut(eachInfo,maskImage,pictureSize)
needDivide=preDivideSolution(eachInfo,maskImage,pictureSize);
if needDivide==0
    canDivideorNot=1;
    afterDivide=eachInfo;
    return
end
inputNum=numel(eachInfo);
realIn=0;
inImage=[];
for i=1:inputNum
    [eachNum,xyMin,regionImage]=findRegionNum(eachInfo{i},pictureSize);
    if eachNum>=2
        canDivideorNot=0;
        afterDivide=[];
    end
    in{i}.xyMin=xyMin;
    in{i}.BWImage=regionImage{1};
    if ~isempty(in{i}.BWImage)
        realIn=realIn+1;
    end
    inImage=[inImage;eachInfo{i}];
end
for i=1:size(in,2)
    image=imfill(in{i}.BWImage,'holes');
    if numel(image(image==1))>700
        afterDivide=[];
        canDivideorNot=0;
        return
    end
end
[outputNum,xyMin,regionImage]=findRegionNum(maskImage,pictureSize);
for i=1:outputNum
    out{i}.xyMin=xyMin;
    out{i}.BWImage=regionImage{i};
end
if realIn<outputNum
    canDivideorNot=false;
    afterDivide=[];
    return
end
if realIn==outputNum
    [afterDivide,canDivideorNot]=n2mNodeSolution(in,out,pictureSize);
    return
end
if realIn>outputNum
    [inArea,minArea]=getSumArea(inImage,in,inputNum,pictureSize);
    [outArea,~]=getSumArea(maskImage,[],outputNum,pictureSize);
    if (inArea-outArea)/(realIn-outputNum)/minArea>=0.6
        [afterDivide,canDivideorNot]=n2mNodeSolution(in,out,pictureSize);
    else
        for i=1:size(out,2)
            out{1}.BWImage=out{1}.BWImage|out{i}.BWImage;
        end
        [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize);
        memberNum=0;
        if canDivideorNot==1
            for i=1:size(afterDivide,2)
                if ~isempty(afterDivide{i})
                    memberNum=memberNum+1;
                end
            end
            if inArea-outArea<150 && memberNum~=realIn || inArea-outArea>250 && memberNum==realIn
                canDivideorNot=0;
            end
            if memberNum==0
                canDivideorNot=0;
            end
        end
    end
end
end
%% this  four function is used to divide N to 1 Node,by rotating the model and erode the out image
function [afterDivide,canDivideorNot]=n21NodeDivide(in,out,pictureSize)
sita=-30:2:30;
outInfo=out{1};
outImage=outInfo.BWImage;
outImage=imfill(outImage,'holes');
xyMin=outInfo.xyMin;
outEndImage=bwmorph(outImage,'thin',inf);
outEndImage=bwmorph(outEndImage,'endpoints');
[xInfo,yInfo]=find(outEndImage==1);
outEndInfo(:,1)=xInfo;
outEndInfo(:,2)=yInfo;
if size(outEndInfo,1)>=numel(in)
    [afterDivide,canDivideorNot]=ASolution(in,outInfo,pictureSize,outEndInfo);
else
    [afterDivide,canDivideorNot]=BSolution(in,sita,outInfo,outImage,pictureSize,xyMin);
end
end
%% this function is used to divide N to M Node directly,choose the nearest bacteria
function [afterDivide,isN2mNode]=n2mNodeSolution(in,out,pictureSize)
% used to deal with N-to-M node and divide them,here N>=M
threShold=30;
inNum=numel(in);
for i=1:inNum
    if ~isempty(in{i}.BWImage)
    Info=regionprops(in{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=in{i}.xyMin;
    inInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
    else
        inInfo(i,:)=zeros(1,4);
    end
end
outNum=numel(out);
for i=1:numel(out)
    Info=regionprops(out{i}.BWImage,'Centroid','MajorAxisLength','MinorAxisLength');
    xyMin=out{i}.xyMin;
    outInfo(i,:)=[Info.Centroid(2)+xyMin(1)-1,Info.Centroid(1)+xyMin(2)-1,Info.MajorAxisLength,Info.MinorAxisLength];
end
distanceMat=zeros(inNum);
distanceMat(1:inNum,1:outNum)=pdist2(inInfo,outInfo);
possibleSol=perms(1:inNum);
caculateDist=zeros(size(possibleSol,1),1);
for i=1:size(possibleSol,1)
    caculateDist(i)=trace(distanceMat(:,possibleSol(i,:)));
end
[minDist,Num]=min(caculateDist);
if minDist>threShold*outNum
    isN2mNode=false;
    afterDivide=[];
else
    isN2mNode=true;
    correctVec=possibleSol(Num,:);
    for i=1:inNum
        if correctVec(i)>outNum
            afterDivide{i}=[];
        else
            afterDivide{i}=xy2Idx(out{correctVec(i)}.xyMin,out{correctVec(i)}.BWImage,pictureSize);
        end
    end
end
end
function [sumArea,minArea]=getSumArea(idxInfo,info,num,pictureSize)
[~,BWImage]=idx2Xy(idxInfo,pictureSize);
BWImage=imfill(BWImage,'holes');
sumArea=numel(BWImage(BWImage==1));
minArea=[];
if ~isempty(info)
    areaInfo=zeros(num,1);
    for i=1:num
        filledImage=imfill(info{i}.BWImage,'holes');
        %     filledImage=info{i}.BWImage;
        areaInfo(i)=numel(filledImage(filledImage==1));
    end
    minArea=min(areaInfo(areaInfo~=0));
end
end
function [regionNum,xyMin,regionImage]=findRegionNum(pixelIdxList,pictureSize)
% find how many regions there are in a input/output
if isempty(pixelIdxList)
    regionNum=0;
    xyMin=zeros(1,2);
    regionImage{1}=[];
else
    [xyMin,BWImage]=idx2Xy(pixelIdxList,pictureSize);
    CC=bwconncomp(BWImage);
    regionNum=CC.NumObjects;
    for i=1:regionNum
        pixelIdxList2=CC.PixelIdxList{i};
        BW=false(CC.ImageSize);
        BW(pixelIdxList2)=1;
        regionImage{i}=BW;
    end
end
end
function needDivide=preDivideSolution(eachInfo,maskImage,pictureSize)
for iIn=2:numel(eachInfo)
    eachInfo{1}=[eachInfo{1};eachInfo{iIn}];
end
imageIn=false(pictureSize);
imageOut=false(pictureSize);
imageIn(eachInfo{1})=1;
imageIn=imfill(imageIn,'holes');
areaIn=numel(imageIn(imageIn==1));
imageOut(maskImage)=1;
imageOut=imfill(imageOut,'holes');
areaOut=numel(imageOut(imageOut==1));
imageOve=imageIn&imageOut;
imageOve=imfill(imageOve,'holes');
areaOve=numel(imageOve(imageOve==1));
if max(areaIn-areaOve,areaOut-areaOve)<=20
    needDivide=0;
else
    needDivide=1;
end
end