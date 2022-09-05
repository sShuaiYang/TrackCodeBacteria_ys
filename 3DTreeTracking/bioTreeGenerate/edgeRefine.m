function [bioTree,edgeNumber]=edgeRefine(bioTree,frameShift) %cut loop edges in the tree
edgeNumber0=edgeCounter(bioTree,frameShift);
% disp(strcat('orignal edge number=',num2str(edgeNumber)));
for iframe=1+frameShift:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            searching=true;
            iOut=1;
            while searching
                edgeInfo=[iframe,iNode,iOut];
                [isCut,loopList,allNodeList]=isTrueEdge(bioTree,edgeInfo);
                if isCut==true
                    if iframe==5010
                        a=1;
                    end
                    edgeList=findPath(bioTree,edgeInfo,loopList,allNodeList);
                    bioTree=cutEdge(bioTree,edgeInfo,edgeList);                    
                    bioTree=edgeCorrect(bioTree,edgeInfo);  
                    iOut=1;
                else
                    iOut=iOut+1;
                end
                if iOut>size(bioTree{iframe}.node{iNode}.Out,2)
                    searching=false;
                end
            end
        end
    end
end
edgeNumber1=edgeCounter(bioTree,frameShift);
% disp(strcat('afterRefine edge number=',num2str(edgeNumber)));
edgeNumber=[edgeNumber0,edgeNumber1];
end
function [isCut,loopList,allNodeList]=isTrueEdge(bioTree,edgeInfo)
if bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.is2Node==false
    isCut=false;
    loopList=[];
    allNodeList=[];
    return;
end
if bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.is2Node==true
    nextNodeInfo=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo;
    [isLoop,loopList,allNodeList]=serachingLoop(bioTree,edgeInfo,nextNodeInfo);
    if isLoop==true
        isCut=true;
        return;
    else
        isCut=false;
        return;
    end
end
end
function [isLoop,loopList,allNodeList]=serachingLoop(bioTree,edgeInfo,nextNodeInfo)
allNodeList=[];
nodeList=parentNodeList(bioTree,edgeInfo,nextNodeInfo);
while ~isempty(nodeList)
    allNodeList=[allNodeList;nodeList];
    tmep1=nodeList(:,1)-edgeInfo(1);
    tmep2=nodeList(:,2)-edgeInfo(2);
    tmep3=nodeList(:,3)-edgeInfo(3);
    loopList=nodeList(tmep1==0 & tmep2==0 & tmep3~=0,:);
    if ~isempty(loopList)
        isLoop=true;        
        return;
    end
    nodeList=listRuduction(nodeList);
    nodeList=parentNodeList(bioTree,edgeInfo,nodeList);
end
isLoop=false;
loopList=[];
end
function nodeList=parentNodeList(bioTree,edgeInfo,nextNodeInfo)
nodeList=[];
for iList=1:size(nextNodeInfo,1)
    nodeInfo=nextNodeInfo(iList,:);
    for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
            parNodeInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
            if  parNodeInfo(1) >= edgeInfo(1)
                nodeList=[nodeList;parNodeInfo];
            end
        end
    end
end
end
function edgeNumber=edgeCounter(bioTree,frameShift)
edgeNumber=0;
for iframe=1+frameShift:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                edgeNumber=edgeNumber+1;
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            edgeNumber=edgeNumber+size(bioTree{iframe}.node{iNode}.Out,2);
        end
    end
end
end
function edgeList=findPath(bioTree,edgeInfo,loopList,allNodeList)
iList=1;
nodeInfo1=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo;
while true
    edgeList(iList,:)=loopList(1,:);
    nodeInfo2=bioTree{edgeList(iList,1)}.node{edgeList(iList,2)}.Out{edgeList(iList,3)}.nodeInfo;
    if nodeInfo1(1)==nodeInfo2(1) &&  nodeInfo1(2)==nodeInfo2(2)
        return;
    else
        loopList=findNode(nodeInfo2,allNodeList);
        iList=iList+1;
    end
end
end
function loopList=findNode(nodeList,allNodeList)
loopList=[];
temp1=nodeList(1);
temp2=nodeList(2);
loopListTemp=allNodeList((allNodeList(:,1)==temp1 & allNodeList(:,2)==temp2),:);
loopList=[loopList;loopListTemp];
end
function  bioTree=cutEdge(bioTree,edgeInfo,edgeList)
nodeInfo=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo;
bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.is2Node=false;
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=false;
pixelInOldEdge=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.traceInfo.pixelIdxList;
measurmentInOldEdge=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.traceInfo.measurment;
offsetFrame=0;
for iList=1:size(edgeList,1)
    pathInfo=edgeList(iList,:);
    for iTrace=1:size(bioTree{pathInfo(1)}.node{pathInfo(2)}.Out{pathInfo(3)}.traceInfo.pixelIdxList,2)
        if pathInfo(1)==5994 && iTrace==130
            a=1;
        end
        bioTree{pathInfo(1)}.node{pathInfo(2)}.Out{pathInfo(3)}.traceInfo.pixelIdxList{iTrace}=[bioTree{pathInfo(1)}.node{pathInfo(2)}.Out{pathInfo(3)}.traceInfo.pixelIdxList{iTrace};pixelInOldEdge{iTrace+offsetFrame}];
        bioTree{pathInfo(1)}.node{pathInfo(2)}.Out{pathInfo(3)}.traceInfo.measurment{iTrace}=[bioTree{pathInfo(1)}.node{pathInfo(2)}.Out{pathInfo(3)}.traceInfo.measurment{iTrace};measurmentInOldEdge{iTrace+offsetFrame}];
    end
    offsetFrame=offsetFrame+size(bioTree{pathInfo(1)}.node{pathInfo(2)}.Out{pathInfo(3)}.traceInfo.pixelIdxList,2);
end
end
function bioTree=edgeCorrect(bioTree,edgeInfo)
nodeInfo=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo;
bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}=[];
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}=[];
for iOut=1:size(bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out,2)
    if isempty(bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut}) && iOut < size(bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out,2)
        bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut}=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut+1};
        if bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut}.is2Node==false
            leafInfo=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo=[edgeInfo(1),edgeInfo(2),iOut];
        end
        if bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut}.is2Node==true
            nodeInfo1=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut}.nodeInfo;
            bioTree{nodeInfo1(1)}.node{nodeInfo1(2)}.In{nodeInfo1(3)}.nodeInfo=[edgeInfo(1),edgeInfo(2),iOut];
        end
        bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{iOut+1}=[];
    end
end
for iIn=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
    if isempty(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}) && iIn < size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In,2)
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn+1};
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==false
            rootInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
        end
        if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.isNode==true
            nodeInfo1=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn}.nodeInfo;
            bioTree{nodeInfo1(1)}.node{nodeInfo1(2)}.Out{nodeInfo1(3)}.nodeInfo=[nodeInfo(1),nodeInfo(2),iIn];
        end
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{iIn+1}=[];
    end
end
bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out(end)=[];
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In(end)=[];
end
function Listafter=listRuduction(Listbefore)
Listafter=[];
while ~isempty(Listbefore)
    Listafter=[Listafter;Listbefore(1,:)];
    tempListFrame=Listbefore(:,1)-Listbefore(1,1);
    tempListNode=Listbefore(:,2)-Listbefore(1,2);
    Listbefore=Listbefore(tempListFrame~=0|tempListNode~=0,:);
end
end