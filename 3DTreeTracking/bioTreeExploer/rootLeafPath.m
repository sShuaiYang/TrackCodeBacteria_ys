function [rootList,leafList]=rootLeafPath(treeGraph) % exlopr the connectivety in the attaching event
rootList=[];
leafList=[];
for iNode=1:size(treeGraph.Nodes,1)
    if strcmp(treeGraph.Nodes(iNode).Label,'root')
        rootList=[rootList,iNode];
    end
    if strcmp(treeGraph.Nodes(iNode).Label,'leaf')
        leafList=[leafList,iNode];
    end
end
pathM=zeros(size(rootList,2),size(rootList,2));
for iroot=1:size(rootList,2)
    disp(iroot);
    xNode=rootList(iroot);
    yNodeList=rootList;
    distList=shortestPath(treeGraph,xNode,yNodeList);
    pathM(iroot,:)=distList;
end
image(pathM');
end
function distList=shortestPath(treeGraph,xNode,yNodeList)
distList=zeros(1,size(yNodeList,2));
for iNode=1:size(yNodeList,2)
    distList(iNode)= shortestpath(treeGraph,xNode,yNodeList(iNode),'Directed',false);
end
end
