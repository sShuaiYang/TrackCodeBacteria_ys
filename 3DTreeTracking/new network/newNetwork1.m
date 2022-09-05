function clusterTree=newNetwork1(bioTree,bacteriaFrameInfo)
% ��HGT�еõ�����У�ϸ��֮��Ŀ����ṩ��һ���໥���õĸ��ʣ����ڱ�������ı���ʣ�������Ҫ�㣺
% 1.��ÿ��frame������ϸ��֮����룬С��ĳ������ģ��͸���ĳ���ʱ�ʾ��ϵ��˫��
% 2.�����µ����ӣ���a��b��������ϵp,��ô����Ҳ�Ὣ�ͱ��˵���ϵ��p*pi���ݸ��Է�
% 3.���ڷ��ѵ��¼������Ѻ��ϸ�����ɱ�������֮ǰ����ϵ�����Ƕ���֮ǰ��ϵ��ϸ���Ĺ��׼��롣
% 4.��ϵPΪ�����ۼӡ�

% �ʺϸ�����ĳ�����
% 1.�õ�bacteriaFrameInfo��bioTree
% 2.������ڹ�������ϸ��������Ϊmatrix�������ڵ�ǰframe���ڵ�ϸ��������Ӿ���
% 3.�ҵ�ÿ��frame��ϸ��֮����������ӣ�ÿ������Ҳ���������μ���a-b,a��b����һ�Σ�b��a����һ��
% 4.������֮��ʼ�������ݣ���P1>P2��������ݣ�������ѭ�򵥵ļӷ�ԭ��ÿ�ε�matrix�����ӵ���Ӧ����matrix֮�ϡ�

aveLen=60;
multi=1;
step=100;  %�洢ʱ���ò���
[clusterTree,linkInfo]=initializeMatrix(bioTree,bacteriaFrameInfo,aveLen*multi,step);

% �µ��뷨������ÿ��link����Ҫ����̳о����Լ���Ӧ��ϵ���ȰѶ�Ӧ�Ĺ�ϵ��ú�ÿ�����ø����ݽṹ
generateRelation=matrixGenerateRule(bioTree,clusterTree);

% ������õĶ�Ӧ��ϵ�õ������µľ��󣬲������µĹ�ϵ
for iLink=1:size(linkInfo,1)
    disp(iLink)
    clusterTree=accumulateLinkTree(bioTree,clusterTree,linkInfo(iLink,:),step,generateRelation,linkInfo);
end
end
function [clusterTree,linkInfo]=initializeMatrix(bioTree,bacteriaFrameInfo,length,step)
%���ҵ�����leaf�����겢�г�[frame,xPosition,yPosition]
leafAll=[];
for iframe=1:numel(bioTree)-1
    for iLeaf=1:numel(bioTree{iframe}.leavies)
        leafAll=[leafAll;iframe,iLeaf,-1,1,bioTree{iframe}.leavies{iLeaf}.branchIndex,bioTree{iframe}.leavies{iLeaf}.leafMeasurment.Centroid(1),bioTree{iframe}.leavies{iLeaf}.leafMeasurment.Centroid(1)];
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
% �ҵ����е�link���������¼��A��B���ӣ��趨Ϊ [A -> B] [frame, nodeNumA, nodeNumB] [B -> A]
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
% ���ĳ���
% ��һ�����ж����ӵ�����ϸ���ǲ���ͬһ��branch��,������򲻻ᴫ����Ϣ
% ǰ���ǹ��壬���������壬���ͷ���һ��ֻ�ܴӸߵ���͵�����
nodeInfo=clusterTree(linkInfo(1)).nodeInfo;
if nodeInfo(linkInfo(2),5)==nodeInfo(linkInfo(3),5)
    return
end
donorBranch=nodeInfo(linkInfo(2),5);
p=0.001;
preMatrix=sparse(zeros(size(nodeInfo,1)));
preMatrix(linkInfo(2),linkInfo(3))=p;
for iframe=linkInfo(1)+1:numel(bioTree)
    % ��һ����Ҫ��preMatrix����Ϣ���ݵ���һ����һ��frame,��������ѣ���������Ϣ���׼���(����out����Ŀ����)����������ѣ�������Ϣ
    newMatrix=generateNew(clusterTree(iframe).nodeInfo,preMatrix,generateRelation(iframe));
    
    % �ڶ����������µ�����
    newMatrix=updataNewMatrix(newMatrix,full(clusterTree(iframe).linkMatrix),clusterTree(iframe).nodeInfo,donorBranch,p,clusterTree(iframe).currentLink);
    preMatrix=newMatrix;
    if mod(iframe,step)==0 || iframe==numel(bioTree)
        clusterTree(iframe).finalMatrix=sparse(full(clusterTree(iframe).finalMatrix)+preMatrix);
    end
end
end
function generateRelation=matrixGenerateRule(bioTree,clusterTree)
% ������һ��������ǰһ������˳�й�ϵ�ĺ����������������ѭ����������
% nodeInfo��Ϣ��frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition��
% matchOrder=[��Ӧǰһ��frame�ĵڼ���]
% newLink.info1=[�ڼ���,�м���ͬ��]
% newLink.info2=[�̳�˭] ���ڿ��ܴӲ�ͬ��ϸ�����ϼ̳У�F7 ���ݣ�,�̳��������ж��
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
% matrix�̳еĺ���--������Ӧ���¸�frame�ļ̳о���
% nodeInfo��Ϣ��frame, nNode/ nLeaf/ nRoot, 1~n/ -1/ 0, which trace, branchIndex, xPosition, yPosition��
% matchOrder=[��Ӧǰһ��frame�ĵڼ���] 0�����³��ֵ�
% newLink.info1=[�ڼ���,�м���ͬ��]
% newLink.info2=[�̳�˭] ���ڿ��ܴӲ�ͬ��ϸ�����ϼ̳У�F7 ���ݣ�,�̳��������ж��
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