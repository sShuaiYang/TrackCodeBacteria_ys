function bacteriaTime=gainLengthInfo(bioTree,branchList,divideNum)
nodeOrRoot=branchList(:,3);
branchList=branchList(nodeOrRoot==1,:);
bacteriaNum=0;
bacteriaTime=[];
for i=1:size(branchList,1)
    iBranch=branchList(i,:);
    iAllRoot=bioTree{1,iBranch(1)}.node{1,iBranch(2)}.allRoot;
    iAllLeaf=bioTree{1,iBranch(1)}.node{1,iBranch(2)}.allLeaf;
    iAllNode=bioTree{1,iBranch(1)}.node{1,iBranch(2)}.allNode;
    adjMatrix=createAdjacencyMatrix(bioTree,iAllNode);
    for iRoot=1:size(iAllRoot,1)
        iRoot_info=iAllRoot(iRoot,:);
        if bioTree{1,iAllRoot(iRoot,1)}.root{1,iAllRoot(iRoot,2)}.is2Node==1;
            root_Node=bioTree{1,iAllRoot(iRoot,1)}.root{1,iAllRoot(iRoot,2)}.nodeInfo;
            for iLeaf=1:size(iAllLeaf,1)
                iLeaf_info=iAllLeaf(iLeaf,:);
                if bioTree{1,iAllLeaf(iLeaf,1)}.leavies{1,iAllLeaf(iLeaf,2)}.is2Node==1
                    leaf_Node=bioTree{1,iAllLeaf(iLeaf,1)}.leavies{1,iAllLeaf(iLeaf,2)}.nodeInfo;
                    leaf2NodeNum=find(iAllNode(:,1)==leaf_Node(1)&iAllNode(:,2)==leaf_Node(2));
                    if numel(leaf2NodeNum)>=2
                        leaf2NodeNum(2:end)=[];
                    end
                    root2NodeNum=find(iAllNode(:,1)==root_Node(1));
                    if leaf2NodeNum>=root2NodeNum
                        [~,timeSeries,~]=graphshortestpath(adjMatrix,root2NodeNum,leaf2NodeNum);
                        if ~isempty(timeSeries)
                            pType=whos('timeSeries');
                            if strcmp(pType.class,'cell')
                                timeSeries1=timeSeries{1,1};
                            else
                                timeSeries1=timeSeries;
                            end
                            bacteriaNum=bacteriaNum+1;
                            [bacteriaTime,bacteriaNum]=getlengthInfo(bacteriaTime,bacteriaNum,bioTree,timeSeries1,iAllNode,iRoot_info,iLeaf_info,i,divideNum);
                        end
                    end
                else
                    continue
                end
            end
        else
            continue
        end
    end
end
end
function adjMatrix=createAdjacencyMatrix(bioTree,allNode)
dimension=size(allNode,1);
adjMatrix=zeros(dimension);
for i=1:dimension
    adjMatrix(i,i)=1;
    nodeInfo=allNode(i,:);
    nodeToNode=bioTree{1,nodeInfo(1)}.node{1,nodeInfo(2)}.Out;
    for j=1:numel(nodeToNode)
        if bioTree{1,nodeInfo(1)}.node{1,nodeInfo(2)}.Out{1,j}.is2Node==1
        node2Info=bioTree{1,nodeInfo(1)}.node{1,nodeInfo(2)}.Out{1,j}.nodeInfo;
        n=find(allNode(:,1)==node2Info(1));
        adjMatrix(i,n)=1;
        end
    end
end
adjMatrix=sparse(adjMatrix);
end
function [bacteriaTime,num]=getlengthInfo(bacteriaTime,num,bioTree,timeSeries,iAllNode,rootInfo,leafInfo,ibranchList,divideNum)
allNode=iAllNode(timeSeries,:);
allNode=nodeReconfirm(allNode);
dividePoint=numel(find(allNode(:,4)==1));
if dividePoint>=divideNum
    LengthSeries=[];
    beginframe=rootInfo(1);
    endframe=iAllNode(timeSeries(1),1);
    if beginframe~=endframe
        for i=1:endframe-beginframe
            LengthSeries(i,1)=beginframe+i-1;
            LengthSeries(i,2)=bioTree{1,rootInfo(1)}.root{1,rootInfo(2)}.traceInfo.measurment{1,i}(1,1).MajorAxisLength;
            LengthSeries(i,4)=double(bioTree{1,rootInfo(1)}.root{1,rootInfo(2)}.traceInfo.measurment{1,i}(1,1).FilledArea)/LengthSeries(i,2);
        end
        iframe=i;
    else
        iframe=0;
    end
    if iAllNode(timeSeries(1),4)==1
        LengthSeries(iframe+1,3)=1;
    end
    if numel(timeSeries)-1~=0
        if iframe~=0
        infoBeforeOneNode=bioTree{1,rootInfo(1)}.root{1,rootInfo(2)}.traceInfo.pixelIdxList{1,end};
        else
            infoBeforeOneNode=[];
        end
        for t=1:numel(timeSeries)-1
            beginframe=iAllNode(timeSeries(t),1);
            endframe=iAllNode(timeSeries(t+1),1);
            for ibranch=1:numel(bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out)
                if bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out{1,ibranch}.is2Node==1
                    if bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out{1,ibranch}.nodeInfo(1)==endframe
                        break
                    end
                end
            end
            if numel(bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.In)>1
                infoAfterOneNode=bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out{1,ibranch}.traceInfo.pixelIdxList{1,1};
                isgeneral=whethergeneral(infoBeforeOneNode,infoAfterOneNode);
                if isgeneral==0
                    num=num-1;
                    return
                end
            end
            infoBeforeOneNode=bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out{1,ibranch}.traceInfo.pixelIdxList{1,end};
            for i=1:endframe-beginframe
                LengthSeries(i+iframe,1)=iframe+i+rootInfo(1)-1;
                LengthSeries(i+iframe,2)=bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out{1,ibranch}.traceInfo.measurment{1,i}(1,1).MajorAxisLength;
                LengthSeries(i+iframe,4)=double(bioTree{1,iAllNode(timeSeries(t),1)}.node{1,iAllNode(timeSeries(t),2)}.Out{1,ibranch}.traceInfo.measurment{1,i}(1,1).FilledArea)/LengthSeries(i+iframe,2);
            end
            iframe=iframe+i;
            if iAllNode(timeSeries(t+1),4)==1
                LengthSeries(iframe+1,3)=1;
            end
        end
    else
        t=0;
    end
    beginframe=iAllNode(timeSeries(t+1),1);
    endframe=leafInfo(1);
    leaf2Node=bioTree{1,leafInfo(1)}.leavies{1,leafInfo(2)}.nodeInfo;
    for i=1:endframe-beginframe+1
        LengthSeries(i+iframe,1)=iframe+i+rootInfo(1)-1;
%         disp(numel(bioTree{1,leaf2Node(1)}.node{1,leaf2Node(2)}.Out{1,leaf2Node(3)}.traceInfo.measurment)-(endframe-beginframe+1))
        LengthSeries(i+iframe,2)=bioTree{1,leaf2Node(1)}.node{1,leaf2Node(2)}.Out{1,leaf2Node(3)}.traceInfo.measurment{1,i}(1,1).MajorAxisLength;
        LengthSeries(i+iframe,4)=double(bioTree{1,leaf2Node(1)}.node{1,leaf2Node(2)}.Out{1,leaf2Node(3)}.traceInfo.measurment{1,i}(1,1).FilledArea)/LengthSeries(i+iframe,2);
    end
    if size(LengthSeries,2)==2
        LengthSeries(:,3)=0;
    end
    bacteriaTime(1,num).LengthSeries=LengthSeries;
    bacteriaTime(1,num).allNode=allNode;
    bacteriaTime(1,num).branchIndex=ibranchList;
    bacteriaTime(1,num).dividePoint=allNode(allNode(:,4)==1,:);
else
    num=num-1;
end
end
function allNode=nodeReconfirm(allNode)
dividePoint=find(allNode(:,4)==1);
i=1;
while i<numel(dividePoint)
    if allNode(dividePoint(i+1),1)-allNode(dividePoint(i),1)<=10  %% add by jzy 2014.5.19 (pre 300)
        allNode(dividePoint(i+1),4)=0;
        dividePoint(i+1)=[];
        i=i-1;
    end
    i=i+1;
end
end
function isgeneral=whethergeneral(infoBeforeOneNode,infoAfterOneNode)
if ~isempty(infoBeforeOneNode)
    n=0;
    for i=1:numel(infoBeforeOneNode)
        if ~isempty(find(infoAfterOneNode==infoBeforeOneNode(i)))
            n=n+1;
        end
    end
    if double(n/numel(infoBeforeOneNode))>=0.6||double(n/numel(infoAfterOneNode))>=0.6
        isgeneral=1;
    else
        isgeneral=0;
    end
else
    isgeneral=1;
end
end