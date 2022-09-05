function bioTree=binaryTreeFix(bioTree)
dirFile='D:\troubleNodeSaveResult';
mkdir(dirFile);
cd(dirFile)

% bioTree=initiateDivisionNode(bioTree); % give divisionNode mask

troubleBinaryNode=findTroubleDivisionNode(bioTree);
save('troubleBinaryNode.mat','troubleBinaryNode')
load('troubleBinaryNode')
% load('judgeResult')
for iNode=1:numel(troubleBinaryNode)
    clc
    disp(iNode)
    disp(troubleBinaryNode{iNode}.nodeInfo)
    judgeResult{iNode}=judgeNode(troubleBinaryNode{iNode}.pixelIdxListIn,troubleBinaryNode{iNode}.pixelIdxListOut,bioTree{1}.imageSize);
    judgeResult{iNode}.nodeInfo=troubleBinaryNode{iNode}.nodeInfo;
    save('judgeResult.mat','judgeResult')  % save for every 10 node so that you work would not be broken
end
save('judgeResult.mat','judgeResult') 

% bioTree=fakeDivisionNodeMending(bioTree);
bioTree=binaryFix(bioTree,judgeResult);
end
%% find trouble division node
function troubleNodeList=findTroubleDivisionNode(bioTree)
% type 1 针对分裂节点的node的两个细菌的外形进行判断
% troubleNodeList=type1BinaryList(bioTree);

% type 2 根据二叉树的数据结构判断
troubleNodeList=type2BinaryList(bioTree);

% type 3 根据不同branch交汇的情况进行判断
% troubleNodeList=type3BinaryList(bioTree);  %全为1对2时不用
end
function troubleNodeList=type1BinaryList(bioTree)
i=0;
for iframe=1:size(bioTree,2)
    nodeList=troubleNodeinFrame(bioTree,iframe);
    for iNode=1:size(nodeList)
        nodeInfo=nodeList(iNode,:);
        i=i+1;
        troubleNodeList{i}.pixelIdxListIn=getInputMask(bioTree,nodeInfo);
        troubleNodeList{i}.nodeInfo=nodeInfo;
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            troubleNodeList{i}.pixelIdxListOut{iOut}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
        end
    end
end
end
function troubleNodeList=type2BinaryList(bioTree)
minDivisionTime=6;  %% original 800
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
coreBranch=branchList(:,3)==1;
branchList=branchList(coreBranch,:);
nodeList=[];
for ibranch=1:size(branchList,1)
    iList=branchList(ibranch,:);
    allNode=bioTree{iList(1)}.node{iList(2)}.allNode;
    nodeNum=size(allNode,1);
    for iNode=1:size(allNode,1)
        allNode(iNode,4)=morphCorrection(bioTree,allNode(iNode,:));
        bioTree{allNode(iNode,1)}.node{allNode(iNode,2)}.divisionNode=allNode(iNode,4);
    end
    allNode=allNode(allNode(:,4)==1,:);
    for iNode=1:size(allNode,1)
        nodeInfo=allNode(iNode,1:2);
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.isNode==0
            continue
        end
        preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1}.nodeInfo;
        if nodeInfo(1)-preNode(1)<=minDivisionTime || nodeInfo(1)-preNode(1)>=2400
            allNode(iNode,4)=0;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.divisionNode=0;
        end
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
           is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
           if is2Node==1
               nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
               if nextNode(1)-nodeInfo(1)<=minDivisionTime
                   allNode(iNode,4)=0;
                   bioTree{nodeInfo(1)}.node{nodeInfo(2)}.divisionNode=0;
                   break
               end
           end
        end
    end
end
for iframe=1:size(bioTree,2)
    if iframe>=8000
        continue
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        if bioTree{iframe}.node{iNode}.divisionNode==0
            nodeList=[nodeList;iframe,iNode];
        end
    end
end
for iNode=1:size(nodeList)
    nodeInfo=nodeList(iNode,:);
    troubleNodeList{iNode}.pixelIdxListIn=getInputMask(bioTree,nodeInfo);
    troubleNodeList{iNode}.nodeInfo=nodeInfo;
    for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
        troubleNodeList{iNode}.pixelIdxListOut{iOut}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
    end
end
end
function troubleNodeList=type3BinaryList(bioTree)
iList=0;
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        isTroubleNode=0;
        nodeBranch=bioTree{iframe}.node{iNode}.branchIndex;
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            is2Node=bioTree{iframe}.node{iNode}.Out{iOut}.is2Node;
            if is2Node==1
                nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex~=nodeBranch
                    isTroubleNode=1;
                    break
                end
            end
        end
        if isTroubleNode==1
            iList=iList+1;
            nodeInfo=[iframe,iNode];
            troubleNodeList{iList}.nodeInfo=nodeInfo;
            troubleNodeList{iNode}.pixelIdxListIn=getInputMask(bioTree,nodeInfo);
            for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
                troubleNodeList{iNode}.pixelIdxListOut{iOut}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
            end
        end
    end
end
end
function rightNode=morphCorrection(bioTree,nodeInfo)
    pictureSize=bioTree{1}.imageSize;
    if ~(size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)==1 && size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)==2)
        rightNode=0;
        return
    end
    rightNode=1;
end
function nodeList=troubleNodeinFrame(bioTree,iframe)
nodeList=[];
if min(size(bioTree{iframe}.node))~=0
    for iNode=1:size(bioTree{iframe}.node,2)
        if size(bioTree{iframe}.node{iNode}.In,2)==1 && size(bioTree{iframe}.node{iNode}.Out,2)==2
            pixelIdxListIn=getInputMask(bioTree,[iframe,iNode]);
            [~,BWImage]=idx2Xy(pixelIdxListIn{1},bioTree{1}.imageSize);
            pos=regionprops(BWImage,'FilledArea','Centroid','Eccentricity','MajorAxisLength','Orientation','MinorAxisLength');
            InInfo{1}=pos;
            for iOut=1:2
                pixelIdxListOut{iOut}=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{1};
                [~,BWImage]=idx2Xy(pixelIdxListOut{iOut},bioTree{1}.imageSize);
                pos=regionprops(BWImage,'FilledArea','Centroid','Eccentricity','MajorAxisLength','Orientation','MinorAxisLength');
                OutInfo{iOut}=pos;
            end
            if ~(abs(OutInfo{1}.FilledArea-OutInfo{2}.FilledArea)/double(OutInfo{1}.FilledArea)<=0.3 &&...
                    abs(OutInfo{1}.FilledArea-OutInfo{2}.FilledArea)/double(OutInfo{2}.FilledArea)<=0.3 &&...
                    (InInfo{1}.FilledArea-OutInfo{1}.FilledArea)/double(OutInfo{2}.FilledArea)>=0.6 &&...
                    (InInfo{1}.FilledArea-OutInfo{2}.FilledArea)/double(OutInfo{1}.FilledArea)>=0.6 &&...
                    abs(OutInfo{1}.MajorAxisLength-OutInfo{1}.MinorAxisLength)>=5 && abs(OutInfo{2}.MajorAxisLength-OutInfo{2}.MinorAxisLength)>=5)
                nodeList=[nodeList;[iframe,iNode]];
            end
        end
    end
end
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList{end};
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{end};
        end
    end
end
end
function [xyMin,BWImageGain]=idx2Xy(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
xMin=min(xresult);
xMax=max(xresult);
yMin=min(yresult);
yMax=max(yresult);
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=round(xresult2+(yresult2-1)*(xMax-xMin+1));
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
BWImageGain(2:end-1,2:end-1)=BWImage;
xyMin=[xMin,yMin];
end

%% node judge
function result=judgeNode(pixelIdxListIn,pixelIdxListOut,imageSize)
color=colormap(jet(20));
close all
xyMinMax=getAllXyMinMax(pixelIdxListIn,pixelIdxListOut,imageSize);
[~,~,h1]=getColorMarkImage(pixelIdxListIn,color,imageSize,xyMinMax);
movegui(h1,[1,-1]);
[~,~,h2]=getColorMarkImage(pixelIdxListOut,color,imageSize,xyMinMax);
movegui(h2,[-1,-1]);
drawnow;commandwindow;
stateInput=input('‘y’ or ‘n’ _____');
result.state=stateInput;
if strcmp(stateInput,'n')
    matchNum=input('match to which one');
    result.matchNum=matchNum;
else
    result.matchNum=[];
end
end
function xyMinMax=getAllXyMinMax(in,out,imageSize)
pixelAll=[];
for i=1:numel(in)
    pixelAll=[pixelAll;in{i}];
end
for i=1:numel(out)
    pixelAll=[pixelAll;out{i}];
end
[~,~,xyMinMax]=idx2Xy1(pixelAll,imageSize);
end
function [maskImage,xyMin,h]=getColorMarkImage(pixelIdxList,color,imageSize,xyMinMax)
inNum=numel(pixelIdxList);
colorNum=fix(linspace(1,20,inNum));
image=uint8(false(xyMinMax(3)-xyMinMax(1)+3,xyMinMax(4)-xyMinMax(2)+3));
imageAll=cat(3,image,image,image);
maskImage=false(xyMinMax(3)-xyMinMax(1)+3,xyMinMax(4)-xyMinMax(2)+3);
for i=1:inNum
    [in{i}.xyMin,in{i}.BWImage]=idx2Xy1(pixelIdxList{i},imageSize,xyMinMax);
    in{i}.BWImage=bwmorph(in{i}.BWImage,'remove');
    in{i}.Centroid=regionprops(in{i}.BWImage,'centroid');
    image1=imageAll(:,:,1);
    image2=imageAll(:,:,2);
    image3=imageAll(:,:,3);
    image1(in{i}.BWImage==1)=255*color(colorNum(i),1);
    image2(in{i}.BWImage==1)=255*color(colorNum(i),2);
    image3(in{i}.BWImage==1)=255*color(colorNum(i),3);
    imageAll=cat(3,image1,image2,image3);
    maskImage=maskImage | in{i}.BWImage;
end
figure;h=imshow(imageAll);
set(gcf,'Position',[240 195 694 627]);
for i=1:inNum
    text(in{i}.Centroid.Centroid(1),in{i}.Centroid.Centroid(2),num2str(i),'Color',[1,1,1],'FontSize',10);
end
xyMin=in{1}.xyMin;
end
function [xyMin,BWImageGain,xyMinMax]=idx2Xy1(pixelIdxList,pictureSize,xyMinMax)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=round(pixelIdxList-(yresult-1)*xSize);
if nargin==2
    xMin=min(xresult);
    xMax=max(xresult);
    yMin=min(yresult);
    yMax=max(yresult);
    xyMinMax=[xMin,yMin,xMax,yMax];
else
    xMin=xyMinMax(1);
    yMin=xyMinMax(2);
    xMax=xyMinMax(3);
    yMax=xyMinMax(4);
end
xyMin=[xMin,yMin];
yresult2=yresult-yMin+1;
xresult2=xresult-xMin+1;
BWImage=false(xMax-xMin+1,yMax-yMin+1);
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
BWImageGain(2:end-1,2:end-1)=BWImage;
BWImageGain=imfill(BWImageGain,'holes');
end
%% fix bioTree with judgement
function bioTree=binaryFix(bioTree,judgeResult)
nodeList=[];
for iResult=1:numel(judgeResult)
    nodeList=[nodeList;judgeResult{iResult}.nodeInfo];
end
nodeFrame=unique(nodeList(:,1));
for i=1:numel(nodeFrame)
    matchJudge=find(nodeList(:,1)==nodeFrame(i));
    for iNode=1:numel(matchJudge)
        nodeInfo=judgeResult{matchJudge(iNode)}.nodeInfo;
        if strcmp(judgeResult{matchJudge(iNode)}.state,'n')
            bioTree=bioTreeDeleteLeaf(bioTree,nodeInfo,judgeResult{matchJudge(iNode)}.matchNum);
        end        
    end
    bioTree=bioTreeLinkCorrection(bioTree,nodeFrame(i));
end
end
function bioTree=bioTreeDeleteLeaf(bioTree,nodeOut,matchNum)
% generate new root and link it
newRoot=[nodeOut(1),size(bioTree{nodeOut(1)}.root,2)+1];
is2Node=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{3-matchNum}.is2Node;
if is2Node==1
    nextNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{3-matchNum}.nodeInfo;
    bioTree{newRoot(1)}.root{newRoot(2)}.is2Node=1;
    bioTree{newRoot(1)}.root{newRoot(2)}.nodeInfo=nextNode;
    bioTree{newRoot(1)}.root{newRoot(2)}.leafInfo=[];
    bioTree{newRoot(1)}.root{newRoot(2)}.traceInfo.pixelIdxList=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{3-matchNum}.traceInfo.pixelIdxList;
    bioTree{newRoot(1)}.root{newRoot(2)}.rootPixelDetail=bioTree{newRoot(1)}.root{newRoot(2)}.traceInfo.pixelIdxList{1};
    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
    bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=newRoot;
end
if is2Node==0
    nextLeaf=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{3-matchNum}.leafInfo;
    bioTree{newRoot(1)}.root{newRoot(2)}.is2Node=0;
    bioTree{newRoot(1)}.root{newRoot(2)}.nodeInfo=[];
    bioTree{newRoot(1)}.root{newRoot(2)}.leafInfo=nextLeaf;
    bioTree{newRoot(1)}.root{newRoot(2)}.traceInfo.pixelIdxList=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{3-matchNum}.traceInfo.pixelIdxList;
    bioTree{newRoot(1)}.root{newRoot(2)}.rootPixelDetail=bioTree{newRoot(1)}.root{newRoot(2)}.traceInfo.pixelIdxList{1};
    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
    bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=newRoot;
end

% link node pre and node next
preIsNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.In{1}.isNode;
if preIsNode==1
    preNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.In{1}.nodeInfo;
    nextIsNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.is2Node;
    if nextIsNode==1
        nextNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.nodeInfo;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=nextNode;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.traceInfo.pixelIdxList];
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=preNode;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[];
    end
    if nextIsNode==0
        nextLeaf=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.leafInfo;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=0;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[];
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=nextLeaf;
        bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=[bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList,bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.traceInfo.pixelIdxList];
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=preNode;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[];
    end
end
if preIsNode==0
    preRoot=bioTree{nodeOut(1)}.node{nodeOut(2)}.In{1}.rootInfo;
    nextIsNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.is2Node;
    if nextIsNode==1
        nextNode=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.nodeInfo;
        bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
        bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=nextNode;
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.traceInfo.pixelIdxList];
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
        bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=preRoot;
    end
    if nextIsNode==0
        nextLeaf=bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.leafInfo;
        bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=0;
        bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[];
        bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=nextLeaf;
        bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=[bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList,bioTree{nodeOut(1)}.node{nodeOut(2)}.Out{matchNum}.traceInfo.pixelIdxList];
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=preRoot;
    end
end

% delete node
bioTree{nodeOut(1)}.node{nodeOut(2)}=[];
end
function bioTree=bioTreeLinkCorrection(bioTree,iframe)
emptyNode=[];
for iNode=1:numel(bioTree{iframe}.node)
    if isempty(bioTree{iframe}.node{iNode})
        emptyNode=[emptyNode;iNode];
    end
end
bioTree{iframe}.node(emptyNode)=[];
for iNode=1:numel(bioTree{iframe}.node)
    for iIn=1:numel(bioTree{iframe}.node{iNode}.In)
        isNode=bioTree{iframe}.node{iNode}.In{iIn}.isNode;
        if isNode==1
            preNode=bioTree{iframe}.node{iNode}.In{iIn}.nodeInfo;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.leafInfo=[];
            bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[iframe,iNode,iIn];
        end
        if isNode==0
            preRoot=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo;
            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[iframe,iNode,iIn];
            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
        end
    end
    for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
        is2Node=bioTree{iframe}.node{iNode}.Out{iOut}.is2Node;
        if is2Node==1
            nextNode=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=1;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[iframe,iNode,iOut];
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.leafInfo=[];
        end
        if is2Node==0
            nextLeaf=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=1;
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[iframe,iNode,iOut];
            bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.leafInfo=[];
        end
    end
end
end

%% fake division node searching
% function bioTree=fakeDivisionNodeMending(bioTree)
% searchingRange=5;
% for iframe=1:size(bioTree,2)-1
%     preRootSearching=(iframe-searchingRange):(iframe-1);
%     leafList=[];
%     for leafFrame=preRootSearching
%         for iLeaf=1:size(bioTree{leafFrame}.leavies,2)
%             leafList=[leafList;leafFrame,iLeaf];
%         end
%     end
%     if ~isempty(leafList)
%         for iNode=1:size(bioTree{iframe}.node,2)
%             
%         end
%     end
% end