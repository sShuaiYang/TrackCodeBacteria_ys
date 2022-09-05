function bioTree=type3NodeReduction(bioTree,threShold)
nodeCount=0;
for iframe=1:size(bioTree,2)
dispFrame(iframe)
if ~isempty(bioTree{iframe}.node)
    nodeList=type3NodeinFrame(bioTree{iframe},iframe);
    if ~isempty(nodeList)
        for iList=1:size(nodeList,1)
            nodeInfoIn=getInputInfo(bioTree,nodeList(iList,:));
            [traceInfo,inNum,outNum,canDivideorNot]=fullNodeType3Tracking(bioTree,nodeList(iList,:),threShold);
            if canDivideorNot==1
                bioTree=type3NodeLinker(bioTree,nodeInfoIn,traceInfo,inNum,nodeList(iList,:));
                bioTree=removeType3Node(bioTree,nodeList(iList,:),inNum,outNum);
            end
            nodeCount=nodeCount+1;
            %                 disp(nodeCount);
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
function nodeList=type3NodeinFrame(bioTreeFrame,iframe)
nodeList=[];
if min(size(bioTreeFrame.node))~=0
    for iNode=1:size(bioTreeFrame.node,2)
        if (size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==2)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==2 && size(bioTreeFrame.node{iNode}.Out,2)>=3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==4)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==3 && size(bioTreeFrame.node{iNode}.Out,2)==3)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==5)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==5 && size(bioTreeFrame.node{iNode}.Out,2)==4)...
                ||(size(bioTreeFrame.node{iNode}.In,2)==4 && size(bioTreeFrame.node{iNode}.Out,2)==2)
            nodeList=[nodeList;[iframe,iNode]];
        end
    end
end
end
function bioTree=type3NodeLinker(bioTree,nodeInfoIn,traceInfo,inNum,nodeInfo)
if ~(numel(inNum)==1 && inNum==0)
    for iIn=1:numel(inNum)
        if nodeInfoIn{inNum(iIn)}.isNode==false
            rootInfo=nodeInfoIn{inNum(iIn)}.rootInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node=false;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[];
%                 if leafInfo(1)~=rootInfo(1)
                    bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
%                 end
                bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=false;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=rootInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=nodeInfoT;
                bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList=[bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=false;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=[];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.rootInfo=rootInfo;
            end
        end
        if nodeInfoIn{inNum(iIn)}.isNode==true
            nodeInfo=nodeInfoIn{inNum(iIn)}.nodeInfo;
            if traceInfo{iIn}.is2Node==false;
                leafInfo=traceInfo{iIn}.leafInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=false;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[];
%                 if leafInfo(1)~=nodeInfo(1)
                    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
%                 end
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=leafInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node=true;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=nodeInfo;
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[];
                bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail=traceInfo{iIn}.traceInfo.pixelIdxList{end};
            end
            if traceInfo{iIn}.is2Node==true;
                nodeInfoT=traceInfo{iIn}.nodeInfo;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo= nodeInfoT;
                bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList=[bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList,traceInfo{iIn}.traceInfo.pixelIdxList];
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.isNode=true;
                bioTree{nodeInfoT(1)}.node{nodeInfoT(2)}.In{nodeInfoT(3)}.nodeInfo=nodeInfo;
            end
        end
    end
end
if numel(inNum)==1 && inNum==0
    rootNum=size(bioTree{1,nodeInfo(1)}.root,2);
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node=traceInfo{1,1}.is2Node;
    if bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node==0
        leafInfo=traceInfo{1,1}.leafInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.leafInfo=leafInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.nodeInfo=[];
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.rootPixelDetail=traceInfo{1,1}.traceInfo.pixelIdxList{1};
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=traceInfo{1,1}.traceInfo.pixelIdxList;
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.rootInfo=[nodeInfo(1),rootNum+1];
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.nodeInfo=[];
        bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.is2Node=0;
    end
    if bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.is2Node==1
        newNodeInfo=traceInfo{1,1}.nodeInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.nodeInfo=newNodeInfo;
        bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.leafInfo=[];
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.isNode=0;
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.nodeInfo=[];
        bioTree{1,newNodeInfo(1)}.node{1,newNodeInfo(2)}.In{1,newNodeInfo(3)}.rootInfo=[nodeInfo(1),rootNum+1];
    end
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.rootPixelDetail=traceInfo{1,1}.traceInfo.pixelIdxList{1};
    bioTree{1,nodeInfo(1)}.root{1,rootNum+1}.traceInfo.pixelIdxList=traceInfo{1,1}.traceInfo.pixelIdxList;
end
end
function nodeInfoIn=getInputInfo(bioTree,nodeInfo)
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
        nodeInfo_pre=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
        nodeInfoIn{iIn}.isNode=true;
        nodeInfoIn{iIn}.nodeInfo= nodeInfo_pre;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
        rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
        nodeInfoIn{iIn}.isNode=false;
        nodeInfoIn{iIn}.rootInfo= rootInfo;
    end
end
end
function bioTree=removeType3Node(bioTree,nodeInfo,inNum,outNum)
if ~(numel(inNum)==1 && inNum==0)
    bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In(inNum)=[];
    for i=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.isNode==0
            rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo(3)=i;
        end
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.isNode==1
            NewnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{1,i}.nodeInfo;
            bioTree{NewnodeInfo(1)}.node{NewnodeInfo(2)}.Out{1,NewnodeInfo(3)}.nodeInfo(3)=i;
        end
    end
end
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out(outNum)=[];
for i=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.is2Node==0
        leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.leafInfo;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo(3)=i;
    end
    if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.is2Node==1
        NewnodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{1,i}.nodeInfo;
        bioTree{NewnodeInfo(1)}.node{NewnodeInfo(2)}.In{1,NewnodeInfo(3)}.nodeInfo(3)=i;
    end
end
end