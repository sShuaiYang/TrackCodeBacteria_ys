function clusterTree=newNetwork1(bioTree,bacteriaFrameInfo)
% 从HGT中得到的灵感，细菌之间的靠近提供了一个相互作用的概率，用于表征基因改变的率，有以下要点：
% 1.对每张frame，计算细菌之间距离，小于某个距离的，就赋予某个率表示联系（双向）
% 2.对于新的连接，若a和b发生了联系p,那么他们也会将和别人的联系以p*pi传递给对方
% 3.对于分裂的事件，分裂后的细菌依旧保持他们之前的联系，但是对于之前联系的细菌的贡献减半。
% 4.联系P为线性累加。

% 适合该网络的程序框架
% 1.得到bacteriaFrameInfo和bioTree
% 2.对与存在过的所有细菌，排列为matrix，并对于当前frame存在的细菌算出连接矩阵
% 3.找到每个frame的细菌之间的所有连接，每个连接也都分作两次计算a-b,a对b作用一次，b对a作用一次
% 4.作用完之后开始继续传递，若P1>P2则继续传递，传递遵循简单的加法原则，每次的matrix都叠加到对应的总matrix之上。

aveLen=60;
multi=1;
step=100;  %存储时所用步长
[clusterTree,linkInfo]=initializeMatrix(bioTree,bacteriaFrameInfo,aveLen*multi,step);

% 新的想法，由于每个link都需要计算继承矩阵以及对应关系，先把对应的关系算好后，每次引用该数据结构
generateRelation=matrixGenerateRule(bioTree,clusterTree);

% 根据算好的对应关系得到生成新的矩阵，并建立新的关系
for iLink=1:size(linkInfo,1)
    disp(iLink)
    clusterTree=accumulateLinkTree(bioTree,clusterTree,linkInfo(iLink,:),step,generateRelation,linkInfo);
end
end
function [clusterTree,linkInfo]=initializeMatrix(bioTree,bacteriaFrameInfo,length,step)
%　找到所有leaf的坐标并列出[frame,xPosition,yPosition]
leafAll=[];
for iframe=1:numel(bioTree)-1
    for iLeaf=1:numel(bioTree{iframe}.leavies)
        leafAll=[leafAll;iframe,iLeaf,-1,1,bioTree{iframe}.leavies{iLeaf}.branchIndex,bioTree{iframe}.leavies{iLeaf}.leafMeasurment.Centroid(1),bioTree{iframe}.leavies{iLeaf}.leafMeasurment.Centroid(1)];
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
    clusterTree(iframe).linkMatrix=zeros(size(clusterTree(iframe).nodeInfo,1));
    linkMatrix=pdist2(clusterTree(iframe).nodeInfo(1:realNode,6:7),clusterTree(iframe).nodeInfo(1:realNode,6:7));
    linkMatrix(linkMatrix>length)=10000;
    linkMatrix(linkMatrix~=10000)=1;
    linkMatrix(linkMatrix~=1)=0;
    for i=1:size(linkMatrix,1)
        linkMatrix(i,i)=0;
    end
    clusterTree(iframe).linkMatrix(1:realNode,1:realNode)=linkMatrix;
    clusterTree(iframe).linkMatrix=sparse(logical(clusterTree(iframe).linkMatrix));
    if mod(iframe,step)==0 || iframe==numel(bioTree)
        clusterTree(iframe).finalMatrix=sparse(zeros(size(clusterTree(iframe).linkMatrix)));
    end
end

disp('find all link')
% 找到所有的link的情况并记录，A和B连接，设定为 [A -> B] [frame, nodeNumA, nodeNumB] [B -> A]
linkInfo=[];
for iframe=1:numel(bioTree)
    [A,B]=find(clusterTree(iframe).linkMatrix==1);
    if ~isempty(A)
        linkNode(:,2:3)=cat(2,A,B);
        linkNode(:,1)=iframe;
        linkInfo=[linkInfo;linkNode];
        linkNode=[];
        clusterTree(iframe).currentLink=linkNode;
    end
end
disp(strcat('there are',num2str(size(linkInfo,1)),' links'))
end
function clusterTree=accumulateLinkTree(bioTree,clusterTree,linkInfo,step,generateRelation,linkList)
% 核心程序
% 第一步先判定连接的两个细菌是不是同一个branch的,如果是则不会传递信息
% 前面是供体，后面是受体，传送方向单一且只能从高的向低的流动
nodeInfo=clusterTree(linkInfo(1)).nodeInfo;
if nodeInfo(linkInfo(2),5)==nodeInfo(linkInfo(3),5)
    return
end
donorBranch=nodeInfo(linkInfo(2),5);
p=0.001;
preMatrix=sparse(zeros(size(nodeInfo,1)));
preMatrix(linkInfo(2),linkInfo(3))=p;
for iframe=linkInfo(1)+1:numel(bioTree)
    % 第一步，要将preMatrix的信息传递到下一个下一个frame,若供体分裂，则所有信息贡献减半(根据out的数目劈裂)，若受体分裂，则复制信息
    newMatrix=generateNew(clusterTree(iframe).nodeInfo,preMatrix,generateRelation(iframe));
    
    % 第二步，建立新的连接
    newMatrix=updataNewMatrix(newMatrix,full(clusterTree(iframe).linkMatrix),clusterTree(iframe).nodeInfo,donorBranch,p,clusterTree(iframe).currentLink);
    preMatrix=newMatrix;
    if mod(iframe,step)==0 || iframe==numel(bioTree)
        clusterTree(iframe).finalMatrix=sparse(full(clusterTree(iframe).finalMatrix)+preMatrix);
    end
end
end
function generateRelation=matrixGenerateRule(bioTree,clusterTree)
% 计算下一个矩阵与前一个矩阵顺承关系的函数，避免接下来的循环反复计算
% nodeInfo信息（frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition）
% matchOrder=[对应前一个frame的第几个]
% newLink.info1=[第几个,有几个同类]
% newLink.info2=[继承谁] 由于可能从不同的细菌身上继承（F7 数据）,继承数可以有多个
for iframe=2:numel(bioTree)
    preNode=clusterTree((iframe-1)).nodeInfo;
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
            if traceSize>rootInfo(4)
                nextInfo=[rootInfo(1:3),rootInfo(4)+1];
                order=newNode(:,1)==nextInfo(1) & newNode(:,2)==nextInfo(2) & newNode(:,3)==nextInfo(3) & newNode(:,4)==nextInfo(4);
                matchOrder(sequenceNum(order))=iNode;
            end
            if traceSize==rootInfo(4)
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
            if traceSize>nodeInfo(4)
                nextInfo=[nodeInfo(1:3),nodeInfo(4)+1];
                order=newNode(:,1)==nextInfo(1) & newNode(:,2)==nextInfo(2) & newNode(:,3)==nextInfo(3) & newNode(:,4)==nextInfo(4);
                matchOrder(sequenceNum(order))=iNode;
            end
            if traceSize==nodeInfo(4)
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
function newMatrix=generateNew(newNode,preMatrix,generateInfo)
% matrix继承的函数--建立对应的下个frame的继承矩阵
% nodeInfo信息（frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition）
% matchOrder=[对应前一个frame的第几个] 0就是新出现的
% newLink.info1=[第几个,有几个同类]
% newLink.info2=[继承谁] 由于可能从不同的细菌身上继承（F7 数据）,继承数可以有多个
preMatrix=full(preMatrix);
newLink=generateInfo.newLink;
matchOrder=generateInfo.matchOrder;
preMatch=generateInfo.preMatch;
newMatrix=sparse(zeros(size(newNode,1)));
for iNew=1:size(newLink.info1,1)
    sumIn1=sum(preMatrix(newLink.info2{iNew},:),1);
    sumIn2=sum(preMatrix(:,newLink.info2{iNew}),2);
    preMatrix(newLink.info1(iNew,1),:)=sumIn1/newLink.info1(iNew,2);
    preMatrix(1:size(sumIn2,1),newLink.info1(iNew,1))=sumIn2;
    preMatrix(newLink.info1(iNew,1),newLink.info1(iNew,1))=0;
end
preMatrix=sparse(preMatrix);
newMatrix(matchOrder,matchOrder)=preMatrix(preMatch,preMatch);
end
function newMatrix=updataNewMatrix(newMatrix,linkMatrix,nodeInfo,donorBranch,p,currentLink)
addMatrix=zeros(size(linkMatrix));
for i=1:size(currentLink,1)
    if nodeInfo(currentLink(i,3),5)~=donorBranch
        donorLine=newMatrix(:,currentLink(i,2));
        receptorLine=newMatrix(:,currentLink(i,3));
        if max(donorLine)>max(receptorLine)
            addMatrix(:,currentLink(i,3))=receptorLine+donorLine*p;
            if max(donorLine)>max(addMatrix(:,currentLink(i,3)));
                addMatrix(:,currentLink(i,3))=donorLine;
            end
        end
    end
end
newMatrix=addMatrix+newMatrix;
% [A,B]=find(newMatrix~=0);
% linkSequence=1:size(newMatrix,1);
% for i=1:size(A,1)
%     BLine=linkMatrix(:,B(i));
%     linkNum=linkSequence(BLine);
%     for iLink=1:numel(linkNum)
%         if nodeInfo(linkNum(iLink),5)~=donorBranch && newMatrix(A(i),B(i))>newMatrix(A(i),linkNum(iLink))
%             newMatrix(A(i),linkNum(iLink))=newMatrix(A(i),linkNum(iLink))+newMatrix(A(i),B(i))*p;
%             if newMatrix(A(i),linkNum(iLink))>newMatrix(A(i),linkNum(iLink))
%                 newMatrix(A(i),linkNum(iLink))=newMatrix(A(i),linkNum(iLink));
%             end                
%         end
%     end
% end
end