function [linkMatrix,centroidInfo,leafList,leafNum]=generateOneBranchPhytree(bioTree,ibranch)
if ibranch(3)~=0
    allRoot=bioTree{ibranch(1)}.node{ibranch(2)}.allRoot;
else
    allRoot=bioTree{ibranch(1)}.root{ibranch(2)}.allRoot;
end
traceMatrix{1}=allRoot;
i=1;
nodeList=[];
while (size(traceMatrix{i},2)==3 && ~all(traceMatrix{i}(:,3)==0)) || size(traceMatrix{i},2)==2
    i=i+1;
    traceMatrix{i}=[];
    lastList=traceMatrix{i-1};
    for iNum=1:size(lastList,1)
        if i==2
            nextIsNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.is2Node;
            if nextIsNode==0
                nextLeaf=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.leafInfo;
                traceMatrix{i}=[traceMatrix{i};nextLeaf(1),nextLeaf(2),0];
            end
            if nextIsNode==1
                nextNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.nodeInfo;
                traceMatrix{i}=[traceMatrix{i};nextNode(1),nextNode(2),1];
            end
        end
        if i>=3
            if lastList(iNum,3)==0
                traceMatrix{i}=[traceMatrix{i};lastList(iNum,:)];
            end
            if lastList(iNum,3)==1
                outNum=size(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out,2);
                outIsNode=zeros(outNum,1);
                for iOut=1:outNum
                    nextIsNode=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.is2Node;
                    outIsNode(iOut)=nextIsNode;
                    if nextIsNode==0
                        leafInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.leafInfo;
                        traceMatrix{i}=[traceMatrix{i};leafInfo(1),leafInfo(2),0];
                    end
                    if nextIsNode==1
                        nodeInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.nodeInfo;
                        traceMatrix{i}=[traceMatrix{i};nodeInfo(1),nodeInfo(2),1];
                    end
                end
                if any(outIsNode==0)
                    nodeList=[nodeList;lastList(iNum,1:2),1];
                else
                    nodeList=[nodeList;lastList(iNum,1:2),0];
                end
            end
        end
    end
end
allRoot(:,3)=-1;
nodeList=[allRoot;nodeList];
nodeNum=size(nodeList,1);
leafNum=size(traceMatrix{end},1);
% if size(nodeList,1)==leafNum-1
nodeList=nodeList;
leafList=traceMatrix{end};
leafOrder=1:leafNum;
nodeOrder=1:nodeNum;
linkMatrix=zeros(leafNum+nodeNum);
for iNode=1:nodeNum;
    if nodeList(iNode,3)~=-1
        nodeInfo=nodeList(iNode,:);
        outNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
        for iOut=1:outNum
            nextIsNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
            if nextIsNode==0
                leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
                linkMatrix(iNode,nodeNum+orderNum)=1;
            end
            if nextIsNode==1
                nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
                linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=1;
            end
        end
    end
    if nodeList(iNode,3)==-1
        rootInfo=nodeList(iNode,:);
        nextIsNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node;
        if nextIsNode==0
            leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
            orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
            linkMatrix(iNode,nodeNum+orderNum)=1;
        end
        if nextIsNode==1
            nextNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
            orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
            linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=1;
        end
    end
end
% end
if ~isempty(nodeList)
    centroidInfo=[nodeList(:,1);leafList(:,1)];
else
    centroidInfo=leafList(:,1);
end
end