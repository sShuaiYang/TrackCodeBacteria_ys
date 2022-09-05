function result=test_2015_3_30(bioTree,bacteriaFrameInfo)
% function [branchMatrix]=test_2015_3_30(clusterTree,generateRelation,branchList,result)
% function linkMatrix=test_2015_3_30(linkMatrix)
% % 每次连接都有一定的概率，发生transfer
% % 每次连接都有一定的概率，发生transfer
% 多次模拟求概率

aveLen=60;
multi=1;
step=200;  %存储时所用步长
frameStep=4;
characteristicTime=60*3;
% characteristicTime=0;
% dirFile=uigetdir();
% cd(dirFile)
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
branchNum=size(branchList,1);
clusterTree=initializeMatrix(bioTree,bacteriaFrameInfo,aveLen*multi,branchNum);
% 
% % % 新的想法，由于每个link都需要计算继承矩阵以及对应关系，先把对应的关系算好后，每次引用该数据结构
generateRelation=matrixGenerateRule(bioTree,clusterTree,frameStep);
% % 需要存储的是branchList,generateRelation和clusterTree
% 
% % 根据算好的对应关系得到生成新的矩阵，得到记录两个细菌何时相遇的信息的大矩阵
result=accumulateLinkTree(clusterTree,step,generateRelation,size(branchList,1),frameStep);
% 
% % 从这个矩阵中，根据推测的特征的持续时间，得到两个细菌可以发生transfer的时间段
for i=1:numel(result)
    result(i).linkMatrix=getTransferTimer(result(i).finalMatrix,characteristicTime,frameStep);
end

% 根据对应的linkMatrix可以生成bipartite网络，行数据为细菌，列数据为不同的表型

% 1.先根据linkMatrix得到细菌在不同时刻的连接总表
for i=1:numel(result)
    result(i).linkInfo=getLinkInfoFromMatrix(result(i).linkMatrix);
end

% 2.根据linkInfo计算新的联系的网络
for i=1:numel(result)
    result(i).branchMatrix=getBranchMatrix(result(i).linkInfo,clusterTree(result(i).frameNum).nodeInfo,size(branchList,1),frameStep,characteristicTime);
end
end
function clusterTree=initializeMatrix(bioTree,bacteriaFrameInfo,length,branchNum)
%　找到所有leaf的坐标并列出[frame,xPosition,yPosition]
leafAll=[];
for iframe=1:numel(bioTree)-1
    for iLeaf=1:numel(bioTree{iframe}.leavies)
        leafAll=[leafAll;iframe,iLeaf,-1,1,bioTree{iframe}.leavies{iLeaf}.branchIndex,bioTree{iframe}.leavies{iLeaf}.leafMeasurment(1).Centroid(1),bioTree{iframe}.leavies{iLeaf}.leafMeasurment(1).Centroid(2)];
    end
end

% 初始化clusterTree的数据结构，一个linkMatrix配上之前frame的leaf，一个空的linkMatrix,一个用于记录所有叠加信息的finalMatrix,还有一列记录节点信息，坐标
% nodeInfo信息（frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition）

% attention由于每帧都对finalTree进行存储，内存太大，因此得到的cluserTree的连接矩阵需要有个步长step
disp('initialize clusterTree structure')
for iframe=1:numel(bioTree)
    clusterTree(iframe).nodeInfo=cat(2,bacteriaFrameInfo{iframe}.bacteriaInfo,bacteriaFrameInfo{iframe}.centroidInfo);
    realNode=size(clusterTree(iframe).nodeInfo,1);
    clusterTree(iframe).nodeInfo=[clusterTree(iframe).nodeInfo;leafAll(leafAll(:,1)<iframe,:)];
    treeLinkMatrix=zeros(size(clusterTree(iframe).nodeInfo,1));
    linkMatrix=pdist2(clusterTree(iframe).nodeInfo(1:realNode,6:7),clusterTree(iframe).nodeInfo(1:realNode,6:7));
    linkMatrix(linkMatrix>length)=10000;
    linkMatrix(linkMatrix~=10000)=1;
    linkMatrix(linkMatrix~=1)=0;
    for i=1:size(linkMatrix,1)
        linkMatrix(i,i)=0;
    end
    linkMatrix=logical(linkMatrix);
    treeLinkMatrix(1:realNode,1:realNode)=linkMatrix;
    [A,B]=find(treeLinkMatrix==1);
    if ~isempty(A)
        linkNode(:,2:3)=cat(2,A(A<B),B(A<B));  % 重要改动。。应该用于所有的网络，bug
        linkNode(:,1)=iframe;
        clusterTree(iframe).currentLink=linkNode;
        linkNode=[];
    end
end
end
function result=accumulateLinkTree(clusterTree,step,generateRelation,branchNum,frameStep)
% 核心程序
% 第一步先判定连接的两个细菌是不是同一个branch的,如果是则不会传递信息
% 前面是供体，后面是受体，传送方向单一且只能从高的向低的流动
% p=0.0001;
n=0;
for iframe=frameStep:frameStep:numel(clusterTree)
    % 第一步，要将preMatrix的信息传递到下一个下一个frame,若供体分裂，则所有信息贡献减半(根据out的数目劈裂)，若受体分裂，则复制信息
    if iframe==frameStep
        newMatrix=cell(size(clusterTree(iframe).nodeInfo,1),size(clusterTree(iframe).nodeInfo,1));
    else
        newMatrix=generateNew(clusterTree(iframe).nodeInfo,preMatrix,generateRelation(iframe),branchNum);
    end
%     newMatrix=newMatrix | full(clusterTree(iframe).branchMatrix);
    % 第二步，建立新的连接
    newMatrix=updataNewMatrix(newMatrix,clusterTree(iframe).currentLink,iframe);
    preMatrix=newMatrix;
    if iframe==numel(generateRelation) || ismember(iframe,200:200:numel(clusterTree))
        n=n+1;
        result(n).frameNum=iframe;
        result(n).finalMatrix=newMatrix;
        result(n).nodeInfo=clusterTree(iframe).nodeInfo;
    end
end
end
function generateRelation=matrixGenerateRule(bioTree,clusterTree,frameStep)
% 计算下一个矩阵与前一个矩阵顺承关系的函数，避免接下来的循环反复计算
% nodeInfo信息（frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition）
% matchOrder=[对应前一个frame的第几个]
% newLink.info1=[第几个,有几个同类]
% newLink.info2=[继承谁] 由于可能从不同的细菌身上继承（F7 数据）,继承数可以有多个
for iframe=2*frameStep:frameStep:numel(bioTree)
    preNode=clusterTree(iframe-frameStep).nodeInfo;
    newNode=clusterTree(iframe).nodeInfo;
    matchOrder=zeros(size(newNode,1),1);
    newLink.info1=[];
    newLink.info2=[];
    newSize=size(preNode,1);
    sequenceNum=1:size(newNode,1);
    preSequence=1:size(preNode,1);
    for iNode=1:size(preNode,1);
        if preNode(iNode,3)==-1
            order=newNode(:,1)==preNode(iNode,1) & newNode(:,2)==preNode(iNode,2) & newNode(:,3)==preNode(iNode,3);
            matchOrder(sequenceNum(order))=iNode;
        end
        if preNode(iNode,3)==0
            rootInfo=preNode(iNode,1:4);
            traceSize=numel(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList);
            if traceSize>=rootInfo(4)+frameStep
                nextInfo=[rootInfo(1:3),rootInfo(4)+frameStep];
                order=newNode(:,1)==nextInfo(1) & newNode(:,2)==nextInfo(2) & newNode(:,3)==nextInfo(3) & newNode(:,4)==nextInfo(4);
                matchOrder(sequenceNum(order))=iNode;
            end
            if traceSize<rootInfo(4)+frameStep
                is2Node=bioTree{rootInfo(1)}.root{rootInfo(2)}.is2Node;
                if is2Node==0
                    leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
                    order=newNode(:,1)==leafInfo(1) & newNode(:,2)==leafInfo(2) & newNode(:,3)==-1;
                    matchOrder(sequenceNum(order))=iNode;
                end
                if is2Node==1
                    nextInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeInfo;
                    inOrder=[];
                    for iIn=1:size(bioTree{nextInfo(1)}.node{nextInfo(2)}.In,2)
                        isNode=bioTree{nextInfo(1)}.node{nextInfo(2)}.In{iIn}.isNode;
                        if isNode==0
                            preRoot=bioTree{nextInfo(1)}.node{nextInfo(2)}.In{iIn}.rootInfo;
                            order=preNode(:,1)==preRoot(1) & preNode(:,2)==preRoot(2) & preNode(:,3)==0;
                            inOrder=[inOrder;preSequence(order)];
                        end
                        if isNode==1
                            inNode=bioTree{nextInfo(1)}.node{nextInfo(2)}.In{iIn}.nodeInfo;
                            order=preNode(:,1)==inNode(1) & preNode(:,2)==inNode(2) & preNode(:,3)==inNode(3);
                            inOrder=[inOrder;preSequence(order)];
                        end
                    end
                    for iOut=1:size(bioTree{nextInfo(1)}.node{nextInfo(2)}.Out,2)
                        order=newNode(:,1)==nextInfo(1) & newNode(:,2)==nextInfo(2) & newNode(:,3)==iOut;
                        if iOut==1
                            matchOrder(sequenceNum(order))=iNode;
                            newLink.info1=[newLink.info1;iNode,size(bioTree{nextInfo(1)}.node{nextInfo(2)}.Out,2)];
                            newLink.info2{size(newLink.info1,1)}=inOrder;
                        else
                            newSize=newSize+1;
                            newLink.info1=[newLink.info1;newSize,size(bioTree{nextInfo(1)}.node{nextInfo(2)}.Out,2)];
                            newLink.info2{size(newLink.info1,1)}=inOrder;
                            matchOrder(sequenceNum(order))=newSize;
                        end
                    end
                end
            end
        end
        if preNode(iNode,3)>=1
            nodeInfo=preNode(iNode,1:4);
            traceSize=numel(bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList);
            if traceSize>=nodeInfo(4)+frameStep
                nextInfo=[nodeInfo(1:3),nodeInfo(4)+frameStep];
                order=newNode(:,1)==nextInfo(1) & newNode(:,2)==nextInfo(2) & newNode(:,3)==nextInfo(3) & newNode(:,4)==nextInfo(4);
                matchOrder(sequenceNum(order))=iNode;
            end
            if traceSize<nodeInfo(4)+frameStep
                is2Node=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node;
                if is2Node==0
                    leafInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo;
                    order=newNode(:,1)==leafInfo(1) & newNode(:,2)==leafInfo(2) & newNode(:,3)==-1;
                    matchOrder(sequenceNum(order))=iNode;
                end
                if is2Node==1
                    nextInfo=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo;
                    inOrder=[];
                    for iIn=1:size(bioTree{nextInfo(1)}.node{nextInfo(2)}.In,2)
                        isNode=bioTree{nextInfo(1)}.node{nextInfo(2)}.In{iIn}.isNode;
                        if isNode==0
                            preRoot=bioTree{nextInfo(1)}.node{nextInfo(2)}.In{iIn}.rootInfo;
                            order=preNode(:,1)==preRoot(1) & preNode(:,2)==preRoot(2) & preNode(:,3)==0;
                            inOrder=[inOrder;preSequence(order)];
                        end
                        if isNode==1
                            inNode=bioTree{nextInfo(1)}.node{nextInfo(2)}.In{iIn}.nodeInfo;
                            order=preNode(:,1)==inNode(1) & preNode(:,2)==inNode(2) & preNode(:,3)==inNode(3);
                            inOrder=[inOrder;preSequence(order)];
                        end
                    end
                    for iOut=1:size(bioTree{nextInfo(1)}.node{nextInfo(2)}.Out,2)
                        order=newNode(:,1)==nextInfo(1) & newNode(:,2)==nextInfo(2) & newNode(:,3)==iOut;
                        if iOut==1
                            matchOrder(sequenceNum(order))=iNode;
                            newLink.info1=[newLink.info1;iNode,size(bioTree{nextInfo(1)}.node{nextInfo(2)}.Out,2)];
                            newLink.info2{size(newLink.info1,1)}=inOrder;
                        else
                            newSize=newSize+1;
                            newLink.info1=[newLink.info1;newSize,size(bioTree{nextInfo(1)}.node{nextInfo(2)}.Out,2)];
                            newLink.info2{size(newLink.info1,1)}=inOrder;
                            matchOrder(sequenceNum(order))=newSize;
                        end
                    end
                end
            end
        end
    end
    generateRelation(iframe).newLink=newLink;
    generateRelation(iframe).preMatch=matchOrder;
    generateRelation(iframe).preMatch(generateRelation(iframe).preMatch==0)=[];
    generateRelation(iframe).matchOrder=matchOrder~=0;
end
end
function newMatrix=generateNew(newNode,preMatrix,generateInfo,branchNum)
% matrix继承的函数--建立对应的下个frame的继承矩阵
% nodeInfo信息（frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition）
% matchOrder=[对应前一个frame的第几个] 0就是新出现的
% newLink.info1=[第几个,有几个同类]
% newLink.info2=[继承谁] 由于可能从不同的细菌身上继承（F7 数据）,继承数可以有多个
matchOrder=generateInfo.matchOrder;
preMatch=generateInfo.preMatch;
newLink=generateInfo.newLink;
for i=1:size(newLink.info1,1)
    preMatrix(newLink.info1(i,1),:)=preMatrix(newLink.info2{i}(1),:);
    preMatrix(:,newLink.info1(i,1))=preMatrix(:,newLink.info2{i}(1));
end
newMatrix=cell(size(newNode,1),size(newNode,1));
newMatrix(matchOrder,matchOrder)=preMatrix(preMatch,preMatch);
end
function newMatrix=updataNewMatrix(newMatrix,currentLink,iframe)
for iLink=1:size(currentLink,1)
    newMatrix{currentLink(iLink,2),currentLink(iLink,3)}=[newMatrix{currentLink(iLink,2),currentLink(iLink,3)};iframe];
    newMatrix{currentLink(iLink,3),currentLink(iLink,2)}=[newMatrix{currentLink(iLink,3),currentLink(iLink,2)};iframe];
end
end
function linkMatrix=getTransferTimer(finalMatrix,characteristicTime,frameStep)
emptyMatrix=cellfun(@isempty,finalMatrix);
noEmptyIndex=find(emptyMatrix==0);
linkMatrix=cell(size(finalMatrix));
for i=1:numel(noEmptyIndex)
    iCell=finalMatrix{noEmptyIndex(i)};
    chaNum=characteristicTime/(frameStep*3);
    if numel(iCell)>=chaNum
        diffCell=diff(iCell);
        diffCell=diffCell==frameStep;
        cc=bwconncomp(diffCell);
        for iObj=1:cc.NumObjects
            pixel=cc.PixelIdxList{iObj};
            if numel(pixel)>=chaNum-1
                if chaNum<=1
                    linkNum=1:numel(iCell);
                else
                    linkNum=pixel(chaNum-1:end)+1;
                end
                linkMatrix{noEmptyIndex(i)}=[linkMatrix{noEmptyIndex(i)};iCell(linkNum)];
            end
        end
    end
end
end
function linkInfo=getLinkInfoFromMatrix(linkMatrix)
emptyMatrix=cellfun(@isempty,linkMatrix);
emptyMatrix=~emptyMatrix;
emptyMatrix=tril(emptyMatrix,-1);
noEmptyIndex=find(emptyMatrix==1);
matrixSize=size(emptyMatrix,1);
linkInfo=[];
for i=1:numel(noEmptyIndex)
    SmallLinkInfo=linkMatrix{noEmptyIndex(i)};
    SmallLinkInfo(:,2)=fix((noEmptyIndex(i)-1)/matrixSize)+1;
    SmallLinkInfo(:,3)=noEmptyIndex(i)-(SmallLinkInfo(:,2)-1)*matrixSize;
    linkInfo=[linkInfo;SmallLinkInfo];
end
[~,sortOrder]=sort(linkInfo(:,1));
linkInfo=linkInfo(sortOrder,:);
end
function branchMatrix=getBranchMatrix(linkInfo,nodeInfo,branchSize,frameStep,characteristicTime)
preMatrix=ones(size(nodeInfo,1),branchSize)*100000;   % 写成100000用于表示未产生的所有链接
for i=1:size(nodeInfo,1)
    preMatrix(i,nodeInfo(i,5))=-10000;  %　写成-10000用来表示原本的branch
end
for iframe=frameStep:frameStep:max(linkInfo(:,1))
    disp(iframe)
    branchMatrix=preMatrix;
    smallLink=linkInfo(linkInfo(:,1)==iframe,2:3);
    for iLink=1:size(smallLink,1)
        aLink=smallLink(iLink,1);
        bLink=smallLink(iLink,2);
        aLinkInfo=preMatrix(aLink,:);
        bLinkInfo=preMatrix(bLink,:);
        aLinkNew=aLinkInfo;
        bLinkNew=bLinkInfo;
        aLinkNew(aLinkNew>iframe-characteristicTime/(3))=100000;
        aLinkNew(aLinkNew<=iframe-characteristicTime/(3))=iframe;
        bLinkNew(bLinkNew>iframe-characteristicTime/(3))=100000;
        bLinkNew(bLinkNew<=iframe-characteristicTime/(3))=iframe;
        aLinkInfo(aLinkInfo==100000)=bLinkNew(aLinkInfo==100000);
        bLinkInfo(bLinkInfo==100000)=aLinkNew(bLinkInfo==100000);
        branchMatrix(aLink,:)=aLinkInfo;
        branchMatrix(bLink,:)=bLinkInfo;
    end 
    preMatrix=branchMatrix;
end
end