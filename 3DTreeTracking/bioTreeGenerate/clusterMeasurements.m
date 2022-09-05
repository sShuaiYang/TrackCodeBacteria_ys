function [bioTree,connection]=clusterMeasurements(bioTree)
[bioTree,startNodeList,endNodeList,nodeIndex]=nodeTypeCounter(bioTree);
disp(strcat('Start Node is ',num2str(size(startNodeList,1))));
disp(strcat('End Node is ',num2str(size(endNodeList,1))));
disp(strcat('Node Number ',num2str(nodeIndex)));
connection=conectMatrix(bioTree,nodeIndex);
% bioTree=setClusterXY(bioTree,endNodeList);
% plotClusterConnection(bioTree);
end
function [bioTree,startNodeList,endNodeList,nodeIndex]=nodeTypeCounter(bioTree)
startNodeList=[];
endNodeList=[];
nodeIndex=1;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.reduce==0
                startType=bioTree{iframe}.node{iNode}.nodeType(1);
                endType=bioTree{iframe}.node{iNode}.nodeType(2);
                if startType==0 && endType==0
                    bioTree{iframe}.node{iNode}.isCluster=false;
                end
                if startType==0 && endType>0
                    bioTree{iframe}.node{iNode}.isCluster=true;
                    bioTree{iframe}.node{iNode}.nodeIndex=nodeIndex;
                    nodeIndex=nodeIndex+1;
                    startNodeList=[startNodeList;[iframe,iNode,endType]];
                    bioTree{iframe}.node{iNode}.xPosition=[];
                end
                if startType>0 && endType==0
                    bioTree{iframe}.node{iNode}.isCluster=true;
                    bioTree{iframe}.node{iNode}.nodeIndex=nodeIndex;
                    nodeIndex=nodeIndex+1;
                    endNodeList=[endNodeList;[iframe,iNode,startType]];
                    bioTree{iframe}.node{iNode}.xPosition=[];
                end
                if startType>0 && endType>0
                    bioTree{iframe}.node{iNode}.isCluster=true;
                    bioTree{iframe}.node{iNode}.nodeIndex=nodeIndex;
                    nodeIndex=nodeIndex+1;
                    bioTree{iframe}.node{iNode}.xPosition=[];
                end
            end
        end
    end
end
end
function plotClusterConnection(bioTree)
figure;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.reduce==0
                if bioTree{iframe}.node{iNode}.isCluster==true
                    startType=bioTree{iframe}.node{iNode}.nodeType(1);
                    endType=bioTree{iframe}.node{iNode}.nodeType(2);
                    
                    x1=bioTree{iframe}.node{iNode}.xPosition;
                    y1=bioTree{iframe}.node{iNode}.yPosition;
                    if startType==0 && endType>0
                        plot(x1,y1,'MarkerFaceColor',[1 0 0],'MarkerSize',10,'Marker','o','LineStyle','none');hold on;
                    end
                    if startType>0 && endType>0
                        plot(x1,y1,'MarkerFaceColor',[0 1 0],'MarkerSize',5,'Marker','o','LineStyle','none');hold on;
                    end
                    if startType>0 && endType==0
                        plot(x1,y1,'MarkerFaceColor',[0 0 1],'MarkerSize',10,'Marker','o','LineStyle','none');hold on;
                    end
                    
                 
                end
            end
        end
    end
end
end
function bioTree=setClusterXY(bioTree,endNodeList)
for iList=1:size(endNodeList,1)
    frame=endNodeList(iList,1);
    node=endNodeList(iList,2);
    bioTree{frame}.node{node}.xPosition=iList;
    bioTree{frame}.node{node}.yPosition=log(9991-frame);
end
searching=true;
i=1;
while searching        
    [bioTreeAfter,endNodeListAfter,searching]=findParentCluster(bioTree,endNodeList);
    endNodeList=endNodeListAfter;
    bioTree=bioTreeAfter;
end
end
function [bioTree,endNodeListAfter,searching]=findParentCluster(bioTree,endNodeListBefore)
endNodeListAfter=[];
for iList=1:size(endNodeListBefore,1)
    xOld=bioTree{endNodeListBefore(iList,1)}.node{endNodeListBefore(iList,2)}.xPosition;
    for jCluster=1:endNodeListBefore(iList,3)
        nodeInfo=bioTree{endNodeListBefore(iList,1)}.node{endNodeListBefore(iList,2)}.In{jCluster}.nodeInfo;
        parentIn=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeType(1);
        xNew=parentX(xOld,endNodeListBefore(iList,3));
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.xPosition=mean([bioTree{nodeInfo(1)}.node{nodeInfo(2)}.xPosition,xNew(jCluster)]);
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.yPosition=(9991-nodeInfo(2));
        endNodeListAfter=[endNodeListAfter;[nodeInfo(1:2),parentIn]];
    end
end
endNodeListAfter=listRuduction(endNodeListAfter);
if ~isempty(endNodeListAfter)
    searching=true;
    for iList=1:size(endNodeListAfter,1)       
        frame=endNodeListAfter(iList,1);
        node=endNodeListAfter(iList,2);
        bioTree{frame}.node{node}.xPosition=mean(bioTree{frame}.node{node}.xPosition);
        bioTree{frame}.node{node}.yPosition=(9991-frame);
    end
else
    searching=false;
end
end
function xNew=parentX(xOld,sNum)
sepration=0.1;
xNew=zeros(1,sNum);
if sNum/2==fix(sNum/2)
    count=1;
    for i=1:sNum/2
        xNew(i)=xOld-sepration*count;
        count=count+1;
    end
    count=1;
    for i=sNum/2+1:sNum
        xNew(i)=xOld+sepration*count;
        count=count+1;
    end
else
    xNew(fix(sNum/2)+1)=xOld;
    count=1;
    for i=1:fix(sNum/2)
        xNew(i)=xOld-sepration*count;
        count=count+1;
    end
    count=1;
    for i=fix(sNum/2)+2:sNum
        xNew(i)=xOld+sepration*count;
        count=count+1;
    end   
end
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
function connection=conectMatrix(bioTree,nodeIndex)
connection=zeros(nodeIndex,nodeIndex);
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.reduce==0
                if bioTree{iframe}.node{iNode}.isCluster==true
                    nodeIndex1=bioTree{iframe}.node{iNode}.nodeIndex;
                    endType=bioTree{iframe}.node{iNode}.nodeType(2);
                    if endType>0
                        for nextNode=1:endType
                            nodeInfo=bioTree{iframe}.node{iNode}.Out{nextNode}.nodeInfo;
                            nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                            connection(nodeIndex1,nodeIndex2)=1;
                        end
                    end
                          
                end
            end
        end
    end
end

end