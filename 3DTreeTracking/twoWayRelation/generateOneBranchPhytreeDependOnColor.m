%YS °åµ××·×ÙPSL£¬¼Ì³Ð¹ØÏµ
function [linkMatrix,centroidInfo,leafList,leafNum,ageAll,nodeList]=generateOneBranchPhytreeDependOnColor(bioTree,iBra,ibranch,type)
% traceMatrix means the current set of divide or leaf
% (:,3)=0 means not need to trace, =1 means the opposite
% (:,4)=-1 means root, 0 leaf ,1 node.
if strcmp(type,'normal')
    if ibranch(3)~=0
        allRoot=bioTree{ibranch(1)}.node{ibranch(2)}.allRoot;
    else
        allRoot=bioTree{ibranch(1)}.root{ibranch(2)}.allRoot;
    end
    allRoot(:,3)=1;
    allRoot(:,4)=-1;
    traceMatrix{1}=allRoot;
    i=1;
    nodeList=[];
    while  ~all(traceMatrix{i}(:,3)==0)
        i=i+1;
        %     if i==16
        %         p=1;
        %     end
        traceMatrix{i}=[];
        lastList=traceMatrix{i-1};
        for iNum=1:size(lastList,1)
            if i==2
                nextIsNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.is2Node;
                if nextIsNode==0
                    nextLeaf=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.leafInfo;
                    traceMatrix{i}=[traceMatrix{i};nextLeaf(1),nextLeaf(2),0,0];
                end
                if nextIsNode==1
                    nextNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.nodeInfo;
                    traceMatrix{i}=[traceMatrix{i};nextNode(1),nextNode(2),1,1];
                end
            end
            if i>=3
                if lastList(iNum,3)==0
                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,:)];
                end
                if lastList(iNum,3)==1
                    outNum=size(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out,2);
                    noChangeNum=0;
                    for iOut=1:outNum
                        nextIsNode=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.is2Node;
                        if nextIsNode==0
                            leafInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.leafInfo;
                            if ~isempty(traceMatrix{i}) && any(leafInfo(1)==traceMatrix{i}(:,1) & leafInfo(2)==traceMatrix{i}(:,2) & 0==traceMatrix{i}(:,4))
                                noChangeNum=noChangeNum+1;
                                if noChangeNum==1
                                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                                end
                            else
                                traceMatrix{i}=[traceMatrix{i};leafInfo(1),leafInfo(2),0,0];
                            end
                        end
                        if nextIsNode==1
                            nodeInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.nodeInfo;
                            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex~=iBra
                                noChangeNum=noChangeNum+1;
                                if noChangeNum==1
                                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                                end
                                continue
                            end
                            if ~isempty(traceMatrix{i}) && any(nodeInfo(1)==traceMatrix{i}(:,1) & nodeInfo(2)==traceMatrix{i}(:,2) & 1==traceMatrix{i}(:,4))
                                noChangeNum=noChangeNum+1;
                                if noChangeNum==1
                                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                                end
                            else
                                if isempty(nodeList) || ~isempty(nodeList) && ~(any(nodeList(:,1)==nodeInfo(1) & nodeList(:,2)==nodeInfo(2)))
                                    traceMatrix{i}=[traceMatrix{i};nodeInfo(1),nodeInfo(2),1,1];
                                end
                            end
                        end
                    end
                end
            end
        end
        for iNum=1:size(traceMatrix{i},1)
            if traceMatrix{i}(iNum,3)==1 && traceMatrix{i}(iNum,4)==1
                nodeList=[nodeList;traceMatrix{i}(iNum,1),traceMatrix{i}(iNum,2),1,1];
            end
        end
    end
    nodeList=[allRoot;nodeList];
    [~,sortOrder]=sort(nodeList(:,1));
    nodeList=nodeList(sortOrder,:);
    nodeNum=size(nodeList,1);
    % if size(nodeList,1)==leafNum-1
    nodeList=nodeList;
    leafList=traceMatrix{end}(traceMatrix{end}(:,4)==0,:);
    leafNum=size(leafList,1);
    leafOrder=1:leafNum;
    nodeOrder=1:nodeNum;
    linkMatrix=zeros(leafNum+nodeNum);
    for iNode=1:nodeNum
        if nodeList(iNode,4)==1
            nodeInfo=nodeList(iNode,:);
            outNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
            for iOut=1:outNum                                    
                traceRFP=[];
                for iTrace=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.measurment,2)
                    traceRFP=[traceRFP;bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo.measurment{iTrace}.meanRFP];
                end
                nextIsNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
                if nextIsNode==0
                    leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                    orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
                    linkMatrix(iNode,nodeNum+orderNum)=mean(traceRFP);
                end
                if nextIsNode==1
                    nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                    orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
                    linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=mean(traceRFP);
                end
            end
        end
        if nodeList(iNode,4)==-1
            rootInfo=nodeList(iNode,:);
            nextIsNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node;
            traceRFP=[];
            for iTrace=1:size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment,2)
                traceRFP=[traceRFP;bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.measurment{iTrace}.meanRFP];
            end
            if nextIsNode==0
                leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
                orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
                linkMatrix(iNode,nodeNum+orderNum)=mean(traceRFP);
            end
            if nextIsNode==1
                nextNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
                orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
                linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=mean(traceRFP);
            end
        end
    end
    if ~isempty(nodeList)
        centroidInfo=[nodeList(:,1);leafList(:,1)];
    else
        centroidInfo=leafList(:,1);
    end
    ageAll=[];
end
if strcmp(type,'aging')
    if ibranch(3)~=0
        allRoot=bioTree{ibranch(1)}.node{ibranch(2)}.allRoot;
    else
        allRoot=bioTree{ibranch(1)}.root{ibranch(2)}.allRoot;
    end
    allRoot(:,3)=1;
    allRoot(:,4)=-1;
    traceMatrix{1}=allRoot;
    i=1;
    nodeList=[];
    while  ~all(traceMatrix{i}(:,3)==0)
        i=i+1;
        %     if i==16
        %         p=1;
        %     end
        traceMatrix{i}=[];
        lastList=traceMatrix{i-1};
        ageAll=[];
        for iNum=1:size(lastList,1)
            if i==2
                nextIsNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.is2Node;
                if nextIsNode==0
                    nextLeaf=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.leafInfo;
                    traceMatrix{i}=[traceMatrix{i};nextLeaf(1),nextLeaf(2),0,0];
                end
                if nextIsNode==1
                    nextNode=bioTree{lastList(iNum,1)}.root{lastList(iNum,2)}.nodeInfo;
                    traceMatrix{i}=[traceMatrix{i};nextNode(1),nextNode(2),1,1];
                end
                ageAll=[1,1];
            end
            if i>=3
                if lastList(iNum,3)==0
                    ageAll=[ageAll;ageAllPre(iNum)];
                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,:)];
                end
                if lastList(iNum,3)==1
                    outNum=size(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out,2);
                    noChangeNum=0;
                    age=[];
                    for iOut=1:outNum
                        age(iOut)=max(bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.isOld1,bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.isOld2);
                    end
                    [~,sortOrder]=sort(age);
                    orderedAge=ones(outNum,1)*-inf;
                    orderedAge(sortOrder(end))=1;
                    if outNum>=2
                        orderedAge(sortOrder(end-1))=0;
                    end
                    for iOut=sortOrder
                        nextIsNode=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.is2Node;
                        if nextIsNode==0
                            leafInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.leafInfo;
                            if ~isempty(traceMatrix{i}) && any(leafInfo(1)==traceMatrix{i}(:,1) & leafInfo(2)==traceMatrix{i}(:,2) & 0==traceMatrix{i}(:,4))
                                noChangeNum=noChangeNum+1;
                                if noChangeNum==1
                                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                                    ageAll=[ageAll;orderedAge(iOut)];
                                end
                            else
                                traceMatrix{i}=[traceMatrix{i};leafInfo(1),leafInfo(2),0,0];
                                ageAll=[ageAll;orderedAge(iOut)];
                            end
                        end
                        if nextIsNode==1
                            nodeInfo=bioTree{lastList(iNum,1)}.node{lastList(iNum,2)}.Out{iOut}.nodeInfo;
                            if bioTree{nodeInfo(1)}.node{nodeInfo(2)}.branchIndex~=iBra
                                noChangeNum=noChangeNum+1;
                                if noChangeNum==1
                                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                                    ageAll=[ageAll;orderedAge(iOut)];
                                end
                                continue
                            end
                            if ~isempty(traceMatrix{i}) && any(nodeInfo(1)==traceMatrix{i}(:,1) & nodeInfo(2)==traceMatrix{i}(:,2) & 1==traceMatrix{i}(:,4))
                                noChangeNum=noChangeNum+1;
                                if noChangeNum==1
                                    traceMatrix{i}=[traceMatrix{i};lastList(iNum,1),lastList(iNum,2),0,lastList(iNum,4)];
                                    ageAll=[ageAll;orderedAge(iOut)];
                                end
                            else
                                if isempty(nodeList) || ~isempty(nodeList) && ~(any(nodeList(:,1)==nodeInfo(1) & nodeList(:,2)==nodeInfo(2)))
                                    traceMatrix{i}=[traceMatrix{i};nodeInfo(1),nodeInfo(2),1,1];
                                    ageAll=[ageAll;orderedAge(iOut)];
                                end
                            end
                        end
                    end
                end
            end
        end
        ageAllPre=ageAll;
        for iNum=1:size(traceMatrix{i},1)
            if traceMatrix{i}(iNum,3)==1 && traceMatrix{i}(iNum,4)==1
                nodeList=[nodeList;traceMatrix{i}(iNum,1),traceMatrix{i}(iNum,2),1,1];
            end
        end
    end
    nodeList=[allRoot;nodeList];
    [~,sortOrder]=sort(nodeList(:,1));
    nodeList=nodeList(sortOrder,:);
    nodeNum=size(nodeList,1);
    % if size(nodeList,1)==leafNum-1
    nodeList=nodeList;
    leafList=traceMatrix{end}(traceMatrix{end}(:,4)==0,:);
    ageAll=ageAll(traceMatrix{end}(:,4)==0);
    leafNum=size(leafList,1);
    leafOrder=1:leafNum;
    nodeOrder=1:nodeNum;
    linkMatrix=zeros(leafNum+nodeNum);
    for iNode=1:nodeNum
        if nodeList(iNode,4)==1
            nodeInfo=nodeList(iNode,:);
            outNum=size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2);
            for iOut=1:outNum
                nextIsNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
                if nextIsNode==0
                    leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                    orderNum=leafOrder(leafInfo(1)==leafList(:,1) & leafInfo(2)==leafList(:,2));
                    linkMatrix(iNode,nodeNum+orderNum)=max(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.isOld1,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.isOld2); % debug
                end
                if nextIsNode==1
                    nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                    orderNum=nodeOrder(nextNode(1)==nodeList(:,1) & nextNode(2)==nodeList(:,2));
                    linkMatrix(min(orderNum,iNode),max(orderNum,iNode))=max(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.isOld1,bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.isOld2); % debug
                end
            end
        end
        if nodeList(iNode,4)==-1
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
    if ~isempty(nodeList)
        centroidInfo=[nodeList(:,1);leafList(:,1)];
    else
        centroidInfo=leafList(:,1);
    end
end
end