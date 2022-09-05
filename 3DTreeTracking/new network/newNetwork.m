function result=newNetwork(bioTree,bacteriaFrameInfo)
% function result=newNetwork(clusterTree,generateRelation,branchNum)
% ÿ�����Ӷ���һ���ĸ��ʣ�����transfer
% ���ģ�������

aveLen=60;
multi=1;
step=200;  %�洢ʱ���ò���
frameStep=20;
dirFile=uigetdir();
cd(dirFile)
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
branchNum=size(branchList,1);
clusterTree=initializeMatrix(bioTree,bacteriaFrameInfo,aveLen*multi,branchNum);

% �µ��뷨������ÿ��link����Ҫ����̳о����Լ���Ӧ��ϵ���ȰѶ�Ӧ�Ĺ�ϵ��ú�ÿ�����ø����ݽṹ
generateRelation=matrixGenerateRule(bioTree,clusterTree,frameStep);

% ������õĶ�Ӧ��ϵ�õ������µľ��󣬲������µĹ�ϵ
n=0;
number=0;
for i=1:40000
    n=n+1;
    result{n}=accumulateLinkTree(clusterTree,step,generateRelation,branchNum,frameStep);
    if mod(n,20)==0
        number=number+1;
        save(strcat(num2str(number),'.mat'),'result');
        clear result
        n=0;
    end
end
end
function clusterTree=initializeMatrix(bioTree,bacteriaFrameInfo,length,branchNum)
%���ҵ�����leaf�����겢�г�[frame,xPosition,yPosition]
leafAll=[];
for iframe=1:numel(bioTree)-1
    for iLeaf=1:numel(bioTree{iframe}.leavies)
        leafAll=[leafAll;iframe,iLeaf,-1,1,bioTree{iframe}.leavies{iLeaf}.branchIndex,bioTree{iframe}.leavies{iLeaf}.leafMeasurment(1).Centroid(1),bioTree{iframe}.leavies{iLeaf}.leafMeasurment(1).Centroid(2)];
    end
end

% ��ʼ��clusterTree�����ݽṹ��һ��linkMatrix����֮ǰframe��leaf��һ���յ�linkMatrix,һ�����ڼ�¼���е�����Ϣ��finalMatrix,����һ�м�¼�ڵ���Ϣ������
% nodeInfo��Ϣ��frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition��

% attention����ÿ֡����finalTree���д洢���ڴ�̫����˵õ���cluserTree�����Ӿ�����Ҫ�и�����step
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
        linkNode(:,2:3)=cat(2,A,B);
        linkNode(:,1)=iframe;
        clusterTree(iframe).currentLink=linkNode(1:2:end,:);
        linkNode=[];
    end
    branchMatrix=false(size(clusterTree(iframe).nodeInfo,1),branchNum);
    for i=1:size(clusterTree(iframe).nodeInfo,1)
        branchMatrix(i,clusterTree(iframe).nodeInfo(i,5))=1;
    end
    clusterTree(iframe).branchMatrix=sparse(branchMatrix);
end

% disp('find all link')
% �ҵ����е�link���������¼��A��B���ӣ��趨Ϊ [A -> B] [frame, nodeNumA, nodeNumB] [B -> A]
% for iframe=1:numel(bioTree)
%     [A,B]=find(full(clusterTree(iframe).linkMatrix)==1);
%     if ~isempty(A)
%         linkNode(:,2:3)=cat(2,A,B);
%         linkNode(:,1)=iframe;
%         clusterTree(iframe).currentLink=linkNode(1:2:end,:);
%         linkNode=[];
%     end
% end
end
function result=accumulateLinkTree(clusterTree,step,generateRelation,branchNum,frameStep)
% ���ĳ���
% ��һ�����ж����ӵ�����ϸ���ǲ���ͬһ��branch��,������򲻻ᴫ����Ϣ
% ǰ���ǹ��壬���������壬���ͷ���һ��ֻ�ܴӸߵ���͵�����
% p=0.0001;
p=1;
n=0;
for iframe=frameStep:frameStep:numel(clusterTree)
    % ��һ����Ҫ��preMatrix����Ϣ���ݵ���һ����һ��frame,��������ѣ���������Ϣ���׼���(����out����Ŀ����)����������ѣ�������Ϣ
    if iframe==frameStep
        newMatrix=false(size(clusterTree(iframe).nodeInfo,1),branchNum);
    else
        newMatrix=generateNew(clusterTree(iframe).nodeInfo,preMatrix,generateRelation(iframe),branchNum);
    end
    newMatrix=newMatrix | full(clusterTree(iframe).branchMatrix);
    % �ڶ����������µ�����
    newMatrix=updataNewMatrix(newMatrix,clusterTree(iframe).currentLink,p,clusterTree(iframe).nodeInfo);
    preMatrix=newMatrix;
    if mod(iframe,step)==0 || iframe==numel(clusterTree)
       n=n+1;
       result(n).frameNum=iframe;
       result(n).finalMatrix=sparse(logical(newMatrix-full(clusterTree(iframe).branchMatrix)));
    end
end
end
function generateRelation=matrixGenerateRule(bioTree,clusterTree,frameStep)
% ������һ��������ǰһ������˳�й�ϵ�ĺ����������������ѭ����������
% nodeInfo��Ϣ��frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition��
% matchOrder=[��Ӧǰһ��frame�ĵڼ���]
% newLink.info1=[�ڼ���,�м���ͬ��]
% newLink.info2=[�̳�˭] ���ڿ��ܴӲ�ͬ��ϸ�����ϼ̳У�F7 ���ݣ�,�̳��������ж��
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
% matrix�̳еĺ���--������Ӧ���¸�frame�ļ̳о���
% nodeInfo��Ϣ��frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition��
% matchOrder=[��Ӧǰһ��frame�ĵڼ���] 0�����³��ֵ�
% newLink.info1=[�ڼ���,�м���ͬ��]
% newLink.info2=[�̳�˭] ���ڿ��ܴӲ�ͬ��ϸ�����ϼ̳У�F7 ���ݣ�,�̳��������ж��
matchOrder=generateInfo.matchOrder;
preMatch=generateInfo.preMatch;
newLink=generateInfo.newLink;
for i=1:size(newLink.info1,1)
    preMatrix(newLink.info1(i,1),:)=preMatrix(newLink.info2{i}(1),:);
end
newMatrix=false(size(newNode,1),branchNum);
newMatrix(matchOrder,:)=preMatrix(preMatch,:);
end
function newMatrix=updataNewMatrix(newMatrix,currentLink,p,nodeInfo)
addMatrix=false(size(newMatrix));
for iLink=1:size(currentLink,1)
    line1=newMatrix(currentLink(iLink,2),:);
    line2=newMatrix(currentLink(iLink,3),:);
    for i=1:size(line1,2)
        if rand<=p
            if line1(i)<line2(i) 
                addMatrix(currentLink(iLink,2),i)=1;
            end
            if line1(i)>line2(i)
                addMatrix(currentLink(iLink,3),i)=1;
            end
        end
    end
end
newMatrix=newMatrix | addMatrix;
end