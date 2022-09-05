function allData=generateBinaryTreeAllData(bioTree)
[correctionRateEach,~]=bioTreeCorrectRate(bioTree);
branchNum=1:size(correctionRateEach,1);
focusNum=branchNum(correctionRateEach(:,1)>=6 & (correctionRateEach(:,1)-correctionRateEach(:,2))<=2 & correctionRateEach(:,3)>=0.9);
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
for i=1:numel(focusNum)
    allData(i)=generateAllData(bioTree,focusNum(i),branchList(focusNum(i),:));
end
end
function allData=generateAllData(bioTree,iBra,ibranch)
[linkMatrix,nodeLeafList]=generateOneBranchPhytree(bioTree,iBra,ibranch);
nPiece=0;
pieceLinkInfo=[];
for i=1:size(nodeLeafList.nodeList,1);
    if i==1
        nPiece=nPiece+1;
        rootInfo=nodeLeafList.nodeList(i,:);
        traceInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo;
        pieceInfoDetail=getDetailInfoForPiece(traceInfo);
        pieceInfoDetail.isDetaching=0;
        nextNode=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
        pieceLinkInfo(nPiece,:)=[rootInfo(1),rootInfo(2),-1,nextNode(1),nextNode(2),1];
        allData.pieceDetail(nPiece)=pieceInfoDetail;
    end
    if i~=1
        nodeInfo=nodeLeafList.nodeList(i,:);
        for iOut=1:size(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out,2)
            isNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.is2Node;
            if isNode==1
                nextNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.nodeInfo;
                if bioTree{nextNode(1)}.node{nextNode(2)}.branchIndex==iBra
                    nPiece=nPiece+1;
                    traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo;
                    pieceInfoDetail=getDetailInfoForPiece(traceInfo);
                    pieceInfoDetail.isDetaching=0;
                    pieceLinkInfo(nPiece,:)=[nodeInfo(1),nodeInfo(2),1,nextNode(1),nextNode(2),1];
                    allData.pieceDetail(nPiece)=pieceInfoDetail;
                end
            end
            if isNode==0
                nextLeaf=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.leafInfo;
                nPiece=nPiece+1;
                traceInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{iOut}.traceInfo;
                pieceInfoDetail=getDetailInfoForPiece(traceInfo);
                pieceInfoDetail.isDetaching=1;
                pieceLinkInfo(nPiece,:)=[nodeInfo(1),nodeInfo(2),1,nextLeaf(1),nextLeaf(2),0];
                allData.pieceDetail(nPiece)=pieceInfoDetail;
            end
        end
    end
end
allData.pieceLinkInfo=pieceLinkInfo;
allData.linkMatrix=linkMatrix;
nodeGeneral(1)=1;
for i=1:size(linkMatrix,1)
    for j=i+1:size(linkMatrix,1);
        if linkMatrix(i,j)==1
            nodeGeneral(j)=nodeGeneral(i)+1;
        end
    end
end
allData.allList=cat(1,nodeLeafList.nodeList(:,[1,2,4]),nodeLeafList.leafList(:,[1,2,4]));
for i=1:size(allData.pieceLinkInfo,1)
    matchNode=allData.allList(:,1)==allData.pieceLinkInfo(i,1) & allData.allList(:,2)==allData.pieceLinkInfo(i,2) & allData.allList(:,3)==allData.pieceLinkInfo(i,3);
    edgeGeneral(i)=nodeGeneral(matchNode);
end
allData.edgeGeneral=edgeGeneral;
allData.branchIndex=iBra;
end
function [linkMatrix,nodeLeafList]=generateOneBranchPhytree(bioTree,iBra,ibranch)
% traceMatrix means the current set of divide or leaf
% (:,3)=0 means not need to trace, =1 means the opposite
% (:,4)=-1 means root, 0 leaf ,1 node.
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
nodeLeafList.nodeList=nodeList;
nodeLeafList.leafList=leafList;
end
function pieceInfoDetail=getDetailInfoForPiece(traceInfo)
traceDetail=traceInfo.measurment;
if numel(traceDetail)>=2
    for i=1:numel(traceDetail)
        pieceInfoDetail.positionInfo(i,:)=traceDetail{i}.Centroid;
        pieceInfoDetail.orientation(i,:)=traceDetail{i}.Orientation;
        pieceInfoDetail.lengthInfo(i,:)=traceDetail{i}.MajorAxisLength;
    end
    pieceInfoDetail.positionInfo=pieceInfoDetail.positionInfo(:,end:-1:1);
    pieceInfoDetail.positionInfo(:,1)=wden(pieceInfoDetail.positionInfo(:,1),'sqtwolog','s','sln',3,'sym10');
    pieceInfoDetail.positionInfo(:,2)=wden(pieceInfoDetail.positionInfo(:,2),'sqtwolog','s','sln',3,'sym10');
    
    % caculate MSD
    frameFrequency=1/3;
    timeDelay=logspace(1,log10(2400),100);
    [msd,tao]=getMSD(pieceInfoDetail.positionInfo,timeDelay,frameFrequency);
    p=polyfit(log10(tao(~isnan(msd))),log10(msd(~isnan(msd))),1);
    pieceInfoDetail.MSDInfo=[tao',msd'];
    pieceInfoDetail.MSDslope=p(1);
    
    % basicInfo
    pieceInfoDetail.retentionTime=numel(traceDetail);
    pieceInfoDetail.lengthInfo=wden(pieceInfoDetail.lengthInfo,'heursure','s','one',5,'sym8');
    pieceInfoDetail.averageLength=mean(pieceInfoDetail.lengthInfo);
    pieceInfoDetail.aveGrowthRate=abs(pieceInfoDetail.lengthInfo(end)-pieceInfoDetail.lengthInfo(1))/numel(pieceInfoDetail.lengthInfo);
    
    pieceInfoDetail.velocityInfo=diff(pieceInfoDetail.positionInfo);
    pieceInfoDetail.aveVelocity=mean((pieceInfoDetail.velocityInfo(:,1).^2+pieceInfoDetail.velocityInfo(:,2).^2).^0.5);
end
if  numel(traceDetail)==1
    pieceInfoDetail.positionInfo=traceDetail{1}.Centroid(end:-1:1);
    pieceInfoDetail.orientation=traceDetail{1}.Orientation;   
    pieceInfoDetail.lengthInfo=traceDetail{1}.MajorAxisLength;    
    pieceInfoDetail.MSDInfo=[];
    pieceInfoDetail.MSDslope=[];
    pieceInfoDetail.retentionTime=numel(traceDetail);
    pieceInfoDetail.averageLength=mean(pieceInfoDetail.lengthInfo);
    pieceInfoDetail.aveGrowthRate=[];
    pieceInfoDetail.velocityInfo=[];
    pieceInfoDetail.aveVelocity=[];
end
end
function [msd,tao]=getMSD(position,timeDelay,frameFrequency)
tao=[1,2,3,4,5,6,7,8,9,fix(timeDelay(timeDelay>=10))];
msd=[];
for iTao=1:size(tao,2)
    pos_pre=position(1:end-tao(iTao),:);
    pos_next=position(1+tao(iTao):end,:);
    dataTemp=(pos_next-pos_pre).^2;
    msd=[msd,mean(dataTemp(:,1)+dataTemp(:,2))];
end
% change tao unit to second
tao = tao./frameFrequency;
end
