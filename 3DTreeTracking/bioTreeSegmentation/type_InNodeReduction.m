function bioTree=type_InNodeReduction(bioTree,strength)
nodeCount=0;
for iframe=size(bioTree,2):-1:1
%         p=1;
%     end
    dispFrame(iframe)
    if ~isempty(bioTree{iframe}.node)
        nodeList=getAllList(bioTree{iframe},iframe,bioTree,strength);
        if ~isempty(nodeList)
            for iList=1:size(nodeList,1)
                % nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
                bioTree=fullNodeType_InTracking(bioTree,nodeList(iList,:));
                nodeCount=nodeCount+1;
                % disp(nodeCount);
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
function nodeList=getAllList(bioTreeFrame,iframe,bioTree,strength)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if ~(size(bioTreeFrame.node{iNode}.In,2)==1 && size(bioTreeFrame.node{iNode}.Out,2)==2)
            nodeList=[nodeList;[iframe,iNode]];
        else
            if size(bioTreeFrame.node{iNode}.In,2)==1 && size(bioTreeFrame.node{iNode}.Out,2)==2
                pixelIdxListIn=getInputMask(bioTree,[iframe,iNode]);
                if size(pixelIdxListIn{1},2)<30
                    nodeList=[nodeList;[iframe,iNode]];
                else
                    if strength==2
                        preIsNode=bioTreeFrame.node{iNode}.In{1}.isNode;
                        if preIsNode==1
                            preNodeInfo=bioTreeFrame.node{iNode}.In{1}.nodeInfo;
                            if iframe-preNodeInfo(1)<200
                                if min(size(bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.In,2),size(bioTree{preNodeInfo(1)}.node{preNodeInfo(2)}.Out,2))>=2
                                    nodeList=[nodeList;[iframe,iNode]];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end
function bioTree=fullNodeType_InTracking(bioTree,nodeInfo)%further track the N input and one Output node, return N new complete bacterial trace
pixelIdxListIn=getInputMask(bioTree,nodeInfo);
imageSize=bioTree{1}.imageSize;
pixelIdxListNew=[];
if size(pixelIdxListIn,2)==1
    for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
        [regionNum,~,regionImage]=findRegionNum(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1},imageSize);
        if regionNum~=1
            return
        end
        pixelIdxListOut{1,iOut}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList{1};
    end
    [pixelIdxListNew,canDivideorNot]=basicDivideSolution(pixelIdxListOut,pixelIdxListIn{1}{end},imageSize,1,25);
    if canDivideorNot==0
        return
    else
        pixelIdxListIn{1}{end}=[];
        for i=1:size(pixelIdxListNew,2)
            if isempty(pixelIdxListNew{i})
                [bioTree,pixelIdxListNew]=getNewRoot(bioTree,nodeInfo,pixelIdxListNew);
                if size(pixelIdxListNew,2)==1
                    return
                else
                    break
                end
            end
        end
        %             for i=1:size(pixelIdxListNew,2)
        %                 pixelIdxListIn{1}{end}=[pixelIdxListIn{1}{end};pixelIdxListNew{i}];
        %             end
    end
end
for iIn=1:size(pixelIdxListIn,2)
    bioTree=reductionforOneIn(bioTree,nodeInfo,iIn,pixelIdxListIn,pixelIdxListNew);
end
end
function bioTree=reductionforOneIn(bioTree,nodeInfo,iIn,pixelIdxListIn,pixelIdxListNew)
imageSize=bioTree{1}.imageSize;
% if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==0
%     %         p=1;
%     return
% end
traceInPixelIdxList=pixelIdxListIn{iIn};
pixelIdxListIn=[];
traceInPixelIdxList=cellReversal(traceInPixelIdxList);
if isempty(pixelIdxListNew)
    [regionNum,xyMin,regionImage]=findRegionNum(traceInPixelIdxList{1},imageSize);
    if regionNum==1
        return
    end
    for i=1:regionNum
        image=regionImage{i};
        stats=regionprops(image,'Area','MajorAxisLength','MinorAxisLength');
        if abs(stats.MajorAxisLength*stats.MinorAxisLength*pi/4-stats.Area)>=200
            canDivideorNot=0;
            afterDivide=[];
            return
        end
    end
    for i=1:regionNum
        pixelIdxListIn{i}=xy2Idx(xyMin,regionImage{i},imageSize);
        newTrace{i}.traceInfo.pixelIdxList{1}=pixelIdxListIn{i};
        newTrace{i}.isNode=[];
    end
else
    regionNum=size(pixelIdxListNew,2);
    for i=1:regionNum
        pixelIdxListIn{i}=pixelIdxListNew{i};
        newTrace{i}.traceInfo.pixelIdxList{1}=pixelIdxListNew{i};
        newTrace{i}.isNode=[];
    end
end
if size(traceInPixelIdxList,2)>100
    return
end
for iTrace=2:size(traceInPixelIdxList,2)
    [pixelIdxListNew,canDivideorNot]=basicDivideSolutionforIn(pixelIdxListIn,traceInPixelIdxList{iTrace},imageSize);
    if canDivideorNot==1
        [newTrace,pixelIdxListIn,isBreak,error]=getNewTraceforNode(newTrace,pixelIdxListNew,iTrace,nodeInfo);
        if error==1
            return
        end
        if isBreak==true
            for i=1:size(newTrace,2)
                if  isempty(newTrace{i}.isNode)
                    if iTrace<size(traceInPixelIdxList,2)
                        newTrace{i}.traceInfo.pixelIdxList=[newTrace{i}.traceInfo.pixelIdxList,traceInPixelIdxList(iTrace+1:end)];
                    end
                end
            end
            break
        end
    else
        return
    end
end
if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==1
    preNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
    inNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2);
    outNum=size(bioTree{preNode(1)}.node{preNode(2)}.Out,2);
    oriNode=0;
    for iOut=1:regionNum
        if iOut==1
            if isempty(newTrace{iOut}.isNode)
                oriNode=1;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.is2Node=1;
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode=1;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo=preNode;
            else
                if newTrace{iOut}.isNode==0;
                    rootFrame=newTrace{iOut}.rootInfo(1);
                    rootNum=size(bioTree{rootFrame}.root,2);
                    bioTree{rootFrame}.root{rootNum+1}.is2Node=1;
                    bioTree{rootFrame}.root{1,rootNum+1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                    bioTree{rootFrame}.root{1,rootNum+1}.leafInfo=[];
                    bioTree{rootFrame}.root{1,rootNum+1}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                    bioTree{rootFrame}.root{1,rootNum+1}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode=0;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo=[rootFrame,rootNum+1];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo=[];
                end
            end
        end
        if iOut>1
            if isempty(newTrace{iOut}.isNode)
                if oriNode~=0
                    bioTree{preNode(1)}.node{preNode(2)}.Out{1,outNum+oriNode}.is2Node=1;
                    bioTree{preNode(1)}.node{preNode(2)}.Out{1,outNum+oriNode}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut+inNum-1];
                    bioTree{preNode(1)}.node{preNode(2)}.Out{1,outNum+oriNode}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.isNode=1;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.nodeInfo=[preNode(1),preNode(2),outNum+oriNode];
                    oriNode=oriNode+1;
                else
                    if oriNode==0
                        bioTree{preNode(1)}.node{preNode(2)}.Out{1,preNode(3)}.is2Node=1;
                        bioTree{preNode(1)}.node{preNode(2)}.Out{1,preNode(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut+inNum-1];
                        bioTree{preNode(1)}.node{preNode(2)}.Out{1,preNode(3)}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.isNode=1;
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.nodeInfo=preNode;
                        oriNode=1;
                    end
                end
            else
                if newTrace{iOut}.isNode==0;
                    rootFrame=newTrace{iOut}.rootInfo(1);
                    rootNum=size(bioTree{rootFrame}.root,2);
                    bioTree{rootFrame}.root{1,rootNum+1}.is2Node=1;
                    bioTree{rootFrame}.root{1,rootNum+1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut+inNum-1];
                    bioTree{rootFrame}.root{1,rootNum+1}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                    bioTree{rootFrame}.root{1,rootNum+1}.leafInfo=[];
                    bioTree{rootFrame}.root{1,rootNum+1}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.isNode=0;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.rootInfo=[rootFrame,rootNum+1];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.nodeInfo=[];
                end
            end
        end
    end
else if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==0
        preRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        inNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2);
        oriRoot=0;
        for iOut=1:regionNum
            if iOut==1
                if isempty(newTrace{iOut}.isNode)
                    oriRoot=1;
                    bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
                    bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                    bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                    bioTree{preRoot(1)}.root{preRoot(2)}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                    bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode=0;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo=preRoot;
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo=[];
                else
                    if newTrace{iOut}.isNode==0;
                        rootFrame=newTrace{iOut}.rootInfo(1);
                        rootNum=size(bioTree{rootFrame}.root,2);
                        bioTree{rootFrame}.root{rootNum+1}.is2Node=1;
                        bioTree{rootFrame}.root{1,rootNum+1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
                        bioTree{rootFrame}.root{1,rootNum+1}.leafInfo=[];
                        bioTree{rootFrame}.root{1,rootNum+1}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                        bioTree{rootFrame}.root{1,rootNum+1}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode=0;
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo=[rootFrame,rootNum+1];
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo=[];
                    end
                end
            end
            if iOut>1
                if isempty(newTrace{iOut}.isNode)
                    if oriRoot~=0
                        rootNum=size(bioTree{preRoot(1)}.root,2);
                        bioTree{preRoot(1)}.root{1,rootNum+1}.is2Node=1;
                        bioTree{preRoot(1)}.root{1,rootNum+1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut+inNum-1];
                        bioTree{preRoot(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                        bioTree{preRoot(1)}.root{1,rootNum+1}.leafInfo=[];
                        bioTree{preRoot(1)}.root{1,rootNum+1}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.isNode=0;
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.rootInfo=[preRoot(1),rootNum+1];
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.nodeInfo=[];
                    else
                        if oriRoot==0
                            bioTree{preRoot(1)}.root{preRoot(2)}.is2Node=1;
                            bioTree{preRoot(1)}.root{preRoot(2)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut+inNum-1];
                            bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                            bioTree{preRoot(1)}.root{preRoot(2)}.leafInfo=[];
                            bioTree{preRoot(1)}.root{preRoot(2)}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.isNode=0;
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.rootInfo=preRoot;
                            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.nodeInfo=[];
                            oriRoot=1;
                        end
                    end
                else
                    if newTrace{iOut}.isNode==0;
                        rootFrame=newTrace{iOut}.rootInfo(1);
                        rootNum=size(bioTree{rootFrame}.root,2);
                        bioTree{rootFrame}.root{1,rootNum+1}.is2Node=1;
                        bioTree{rootFrame}.root{1,rootNum+1}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut+inNum-1];
                        bioTree{rootFrame}.root{1,rootNum+1}.traceInfo.pixelIdxList=cellReversal(newTrace{iOut}.traceInfo.pixelIdxList);
                        bioTree{rootFrame}.root{1,rootNum+1}.leafInfo=[];
                        bioTree{rootFrame}.root{1,rootNum+1}.rootPixelDetail=newTrace{iOut}.traceInfo.pixelIdxList{end};
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.isNode=0;
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.rootInfo=[rootFrame,rootNum+1];
                        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,iOut+inNum-1}.nodeInfo=[];
                    end
                end
            end
        end
    end
end
end
function pixelIdxListIn=getInputMask(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        pixelIdxListIn{iIn}=bioTree{nodeInfo_pre(1)}.node{nodeInfo_pre(2)}.Out{nodeInfo_pre(3)}.traceInfo.pixelIdxList;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        if rootInfo(1)==nodeInfo(1)
            pixelIdxListIn{iIn}{1}=bioTree{rootInfo(1)}.root{rootInfo(2)}.rootPixelDetail;
        end
        if  rootInfo(1)<nodeInfo(1)
            pixelIdxListIn{iIn}=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList;
        end
    end
end
end
function [newTrace,pixelIdxListIn,isBreak,error]=getNewTraceforNode(newTrace,pixelIdxListNew,iTrace,nodeInfo)
emptyCount=0;
for iIn=size(newTrace,2)+1:size(pixelIdxListNew,2)
    newTrace{iIn}.traceInfo.pixelIdxList=[];
    newTrace{iIn}.isNode=[];
end
for iIn=1:size(newTrace,2)
    if ~isempty(pixelIdxListNew{iIn});
        if ~isempty(newTrace{iIn}.isNode)
            error=1;
            isBreak=[];
            pixelIdxListIn=[];
            return
        end
        newTrace{iIn}.traceInfo.pixelIdxList=[newTrace{iIn}.traceInfo.pixelIdxList,pixelIdxListNew(iIn)];
    else
        if isempty(newTrace{iIn}.isNode)
            newTrace{iIn}.isNode=false;
            newTrace{iIn}.rootInfo=[nodeInfo(1)-iTrace+1,0];
        end
        emptyCount= emptyCount+1;
    end
    if size(newTrace,2)-emptyCount==1
        isBreak=true;
    else
        isBreak=false;
    end
end
pixelIdxListIn=pixelIdxListNew;
error=0;
end
function [afterDivide,canDivideorNot]=basicDivideSolutionforIn(eachInfo,maskImage,pictureSize)
inputNum=numel(eachInfo);
realIn=0;
inImage=[];
for i=1:inputNum
    realNum(i)=0;
    [eachNum,xyMin,regionImage]=findRegionNum(eachInfo{i},pictureSize);
    if eachNum>=2
        canDivideorNot=0;
        afterDivide=[];
    end
    in{i}.xyMin=xyMin;
    in{i}.BWImage=regionImage{1};
    if ~isempty(in{i}.BWImage)
        realIn=realIn+1;
        realNum(i)=1;
    end
    inImage=[inImage;eachInfo{i}];
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
        [afterDivide(realNum==1),canDivideorNot]=n21NodeDivide(in(realNum==1),out,pictureSize);
        if numel(afterDivide)<numel(in)
            afterDivide{numel(in)}=[];
        end
        memberNum=0;
        if canDivideorNot==1
            for i=1:size(afterDivide,2)
                if ~isempty(afterDivide{i})
                    memberNum=memberNum+1;
                end
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
outInfo=out{1};
[afterDivide,canDivideorNot]=ASolution(in,outInfo,pictureSize);
end
%% this function is used to divide N to M Node directly,choose the nearest bacteria
function [afterDivide,isN2mNode]=n2mNodeSolution(in,out,pictureSize)
% used to deal with N-to-M node and divide them,here N>=M
threShold=15;
inNum=numel(in);
if inNum>8
    afterDivide=[];
    isN2mNode=0;
    return
end
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

function [bioTree,pixelIdxListIn]=getNewRoot(bioTree,nodeInfo,pixelIdxListNew)
emptyNum=[];
for i=1:size(pixelIdxListNew,2)
    if isempty(pixelIdxListNew{i})
        emptyNum=[emptyNum,i];
    end
end
if size(emptyNum,2)==size(pixelIdxListNew,2)
    return
end
t=0;
for i=1:size(pixelIdxListNew,2)
    if isempty(pixelIdxListNew{i})
        rootNum=size(bioTree{nodeInfo(1)}.root,2);
        bioTree{nodeInfo(1)}.root{rootNum+1}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{i};
        bioTree{nodeInfo(1)}.root{rootNum+1}.rootPixelDetail=bioTree{nodeInfo(1)}.root{rootNum+1}.traceInfo.pixelIdxList{1};
        if  bioTree{nodeInfo(1)}.root{rootNum+1}.is2Node==1
            nextNode= bioTree{nodeInfo(1)}.root{rootNum+1}.nodeInfo;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.isNode=0;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[];
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.rootInfo=[nodeInfo(1),rootNum+1];
        else
            if bioTree{nodeInfo(1)}.root{rootNum+1}.is2Node==0;
                nextLeaf=bioTree{nodeInfo(1)}.root{rootNum+1}.leafInfo;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.is2Node=0;
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[];
                bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.rootInfo=[nodeInfo(1),rootNum+1];
            end
        end
    else
        t=t+1;
        pixelIdxListIn{t}=pixelIdxListNew{i};
    end
end
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out(emptyNum)=[];
for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==0
        nextLeaf=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
        bioTree{nextLeaf(1)}.leavies{nextLeaf(2)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut];
    else
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node==1
            nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
            bioTree{nextNode(1)}.node{nextNode(2)}.In{nextNode(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iOut];
        end
    end
end
end
%% here are some basic functions
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
pixelIdxList2=xresult2+(yresult2-1)*(xMax-xMin+1);
BWImage(pixelIdxList2)=true;
BWImageGain=false(xMax-xMin+3,yMax-yMin+3);
xyMin=[xMin,yMin];
BWImageGain(2:end-1,2:end-1)=BWImage;
% BWImageGain=imfill(BWImageGain,'holes');
end
function pixelIdxList=xy2Idx(xyMin,BWImageGain,pictureSize)
%this function is used to convert a small BWImage to its original pixelIdx
BWImage=BWImageGain(2:end-1,2:end-1);
% BWImage=bwmorph(BWImage,'remove');
pixelIdxListOri=find(BWImage==1);
smallxSize=size(BWImage,1);
yresult=ceil(pixelIdxListOri/smallxSize);
xresult=pixelIdxListOri-(yresult-1)*smallxSize;
xresult2=xresult+xyMin(1)-1;
yresult2=yresult+xyMin(2)-1;
xSize=pictureSize(1);
pixelIdxList=xresult2+(yresult2-1)*xSize;
end
function [sumArea,minArea]=getSumArea(idxInfo,info,num,pictureSize)
[~,BWImage]=idx2Xy(idxInfo,pictureSize);
sumArea=numel(BWImage(BWImage==1));
minArea=[];
if ~isempty(info)
    areaInfo=zeros(num,1);
    for i=1:num
        filledImage=info{i}.BWImage;
        %     filledImage=info{i}.BWImage;
        areaInfo(i)=numel(filledImage(filledImage==1));
    end
    minArea=min(areaInfo(areaInfo~=0));
end
end
function A=cellReversal(B)
numB=size(B,2);
for i=1:numB
    A{1,i}=B{1,numB+1-i};
end
end