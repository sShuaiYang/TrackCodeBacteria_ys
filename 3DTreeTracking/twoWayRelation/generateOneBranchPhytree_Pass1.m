function [linkMatrix,centroidInfo]=generateOneBranchPhytree(bioTree,bacteriaList)
linkNode=bacteriaList(:,1:2);
[~,sortOrder]=sort(linkNode(:,1));
linkNode=linkNode(sortOrder,:);
diffData=diff(linkNode);
diffData=[1,1;diffData];
linkNode(diffData(:,1)==0&diffData(:,2)==0,:)=[];
nodeNum=size(bacteriaList,1)-1;
linkMatrix=zeros(nodeNum*2+1);
nodeCentroidInfo=[];
leafCentroidInfo=[];
while nodeNum~=numel(nodeCentroidInfo)
    iNode=linkNode(end,:);
    outNum=size(bioTree{iNode(1)}.node{iNode(2)}.Out,2);
    leafNum=[];
    for iOut=1:outNum
        if bioTree{iNode(1)}.node{iNode(2)}.Out{iOut}.is2Node==0
            leafNum=[leafNum;iOut];
        end
    end
    if numel(leafNum)==2
        for i=1:2
            leafCentroidInfo=[leafCentroidInfo;iNode(1)+size(bioTree{iNode(1)}.node{iNode(2)}.Out{leafNum(i)}.traceInfo.pixelIdxList,2)];
        end
        nodeCentroidInfo=[nodeCentroidInfo;iNode(1)];
        linkMatrix(min(numel(nodeCentroidInfo),nodeNum+numel(leafCentroidInfo)-1),max(numel(nodeCentroidInfo),nodeNum+numel(leafCentroidInfo)-1))=1;
        linkMatrix(min(numel(nodeCentroidInfo),nodeNum+numel(leafCentroidInfo)),max(numel(nodeCentroidInfo),nodeNum+numel(leafCentroidInfo)))=1;
        preIsNode=bioTree{iNode(1)}.node{iNode(2)}.In{1}.isNode;
        if preIsNode==0
            break
        else
            preNode=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo;
            orderNum=linkNode(:,1)==preNode(1)& linkNode(:,2)==preNode(2);
            if any(orderNum==1)
                linkNode(orderNum,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                if size(linkNode,2)==2
                    linkNode(orderNum,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                    linkNode(orderNum,3)=numel(nodeCentroidInfo);
                else
                    for i=3:size(linkNode,2)
                        if linkNode(orderNum,i)==0;
                            linkNode(orderNum,i)=numel(nodeCentroidInfo);
                        else
                            if i==size(linkNode,2)
                                linkNode(orderNum,i+1)=numel(nodeCentroidInfo);
                            end
                        end
                    end
                end
                linkNode(end,:)=[];
            else
                linkNode(end,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                linkNode(end,3)=numel(nodeCentroidInfo);
            end
        end
        [~,sortOrder]=sort(linkNode(:,1));
        linkNode=linkNode(sortOrder,:);
    end
    if numel(leafNum)==1
        leafCentroidInfo=[leafCentroidInfo;iNode(1)+size(bioTree{iNode(1)}.node{iNode(2)}.Out{leafNum(1)}.traceInfo.pixelIdxList,2)];
        nodeCentroidInfo=[nodeCentroidInfo;iNode(1)];
        linkMatrix(min(numel(nodeCentroidInfo),nodeNum+numel(leafCentroidInfo)),max(numel(nodeCentroidInfo),nodeNum+numel(leafCentroidInfo)))=1;
        for i=3:size(linkNode,2)
            if iNode(i)~=0
                linkMatrix(min(numel(nodeCentroidInfo),iNode(i)),max(numel(nodeCentroidInfo),iNode(i)))=1;
            end
        end
        preIsNode=bioTree{iNode(1)}.node{iNode(2)}.In{1}.isNode;
        if preIsNode==0
            break
        else
            preNode=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo;
            orderNum=linkNode(:,1)==preNode(1)& linkNode(:,2)==preNode(2);
            if any(orderNum==1)
                linkNode(orderNum,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                if size(linkNode,2)==2
                    linkNode(orderNum,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                    linkNode(orderNum,3)=numel(nodeCentroidInfo);
                else
                    for i=3:size(linkNode,2)
                        if linkNode(orderNum,i)==0;
                            linkNode(orderNum,i)=numel(nodeCentroidInfo);
                        else
                            if i==size(linkNode,2)
                                linkNode(orderNum,i+1)=numel(nodeCentroidInfo);
                            end
                        end
                    end
                end
                linkNode(end,:)=[];
            else
                linkNode(end,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                linkNode(end,3)=numel(nodeCentroidInfo);
            end
        end
        [~,sortOrder]=sort(linkNode(:,1));
        linkNode=linkNode(sortOrder,:);
    end
    if numel(leafNum)==0;
        nodeCentroidInfo=[nodeCentroidInfo;iNode(1)];
        for i=3:size(linkNode,2)
            if iNode(i)~=0
                linkMatrix(min(numel(nodeCentroidInfo),iNode(i)),max(numel(nodeCentroidInfo),iNode(i)))=1;
            end
        end
        preIsNode=bioTree{iNode(1)}.node{iNode(2)}.In{1}.isNode;
        if preIsNode==0
            break
        else
            preNode=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo;
            orderNum=linkNode(:,1)==preNode(1)& linkNode(:,2)==preNode(2);
            if any(orderNum==1)
                linkNode(orderNum,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                if size(linkNode,2)==2
                    linkNode(orderNum,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                    linkNode(orderNum,3)=numel(nodeCentroidInfo);
                else
                    for i=3:size(linkNode,2)
                        if linkNode(orderNum,i)==0;
                            linkNode(orderNum,i)=numel(nodeCentroidInfo);
                        else
                            if i==size(linkNode,2)
                                linkNode(orderNum,i+1)=numel(nodeCentroidInfo);
                            end
                        end
                    end
                end
                linkNode(end,:)=[];
            else
                linkNode(end,1:2)=bioTree{iNode(1)}.node{iNode(2)}.In{1}.nodeInfo(1:2);
                linkNode(end,3)=numel(nodeCentroidInfo);
            end
        end
        [~,sortOrder]=sort(linkNode(:,1));
        linkNode=linkNode(sortOrder,:);
    end
end
centroidInfo=[nodeCentroidInfo(end:-1:1);leafCentroidInfo];
linkMatrix(1:nodeNum,:)=linkMatrix(nodeNum:-1:1,:);
end

