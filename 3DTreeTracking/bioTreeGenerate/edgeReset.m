function bioTree=edgeReset(bioTree,frameShift)
frameThreshold=20;
for iframe=1+frameShift:size(bioTree,2)
    %         disp(iframe);
    %         if iframe==9394
    %             a=1;
    %         end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            edgeInfo=[iframe,iNode,iOut];
            [isRest,nodeList]=isShortest(bioTree,edgeInfo,frameThreshold);
            if isRest==true
                bioTree=restEdge(bioTree,edgeInfo,loopList);
                bioTree=edgeCorrect(bioTree,edgeInfo);
            end
        end
    end
end
end
function [isRest,nodeList]=isShortest(bioTree,edgeInfo,frameThreshold)
if bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.is2Node==false
    isRest=false;
    nodeList=[];
    return;
end
if bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.is2Node==true
    nextNodeInfo=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo;
    if nextNodeInfo(1)-edgeInfo(1)>frameThreshold
        isRest=false;
        nodeList=[];
        return;
    end
    if nextNodeInfo(1)-edgeInfo(1)<=frameThreshold
       [isRest,nodeList]=serachingParent(bioTree,edgeInfo,nextNodeInfo);
       return;
    end    
end
end
function [isRest,nodeList]=serachingParent(bioTree,edgeInfo,nextNodeInfo)
nodeList=parentNodeList(bioTree,edgeInfo,nextNodeInfo);
tmep1=nodeList(:,1)-edgeInfo(1);
tmep2=nodeList(:,2)-edgeInfo(2);
tmep3=nodeList(:,3)-edgeInfo(3);
nodeList=nodeList(~(tmep1==0 & tmep2==0 & tmep3==0),:);
if ~isempty(nodeList)
    if size(nodeList,1)==1
        isRest=true;
        return;
    else
        isRest=false;
        nodeList=[];
        return;
    end
end
isRest=false;
nodeList=[];
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
function  bioTree=cutEdge(bioTree,edgeInfo,loopList)
nodeInfo=bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo;
bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.is2Node=false;
% bioTree{edgeInfo(1)}.node{edgeInfo(2)}.Out{edgeInfo(3)}.nodeInfo=[];
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.isNode=false;
% bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.nodeInfo=[];
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