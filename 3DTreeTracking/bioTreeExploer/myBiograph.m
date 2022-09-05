function [bioTree,treeGraph]=myBiograph(bioTree,directed) %creat the biograph and the treeplot; true mean dirtied graph; false mean undirected graph and labeld the biograph and the tree
rootColor=[1,0,0];
leafColor=[0,1,0];
endColor=[0,1,0];
nodeColor=[1,1,0];
lineColor=[0.8,0.8,0.8];
MakerSize=3;
[bioTree,cMatrix]=getCMatrix(bioTree,directed);
% [bioTree,cMatrix]=getReductionCMatrix(bioTree,directed);
treeGraph=biograph(cMatrix);
set(treeGraph, 'LayoutScale',1);
% [bioTree,treeGraph]=lableRootandLeaf(bioTree,treeGraph,rootColor,leafColor,endColor,nodeColor);
% dolayout(treeGraph, 'Pathsonly', false);
% plotClusterConnection(bioTree,treeGraph,rootColor,leafColor,endColor,nodeColor,MakerSize,lineColor);
% plotReductionClusterConnection(bioTree,treeGraph,rootColor,leafColor,endColor,nodeColor,MakerSize,lineColor)
% dolayout(treeGraph, 'Paths', true);
% view(treeGraph);
end

function [bioTree,treeGraph]=lableRootandLeaf(bioTree,treeGraph,rootColor,leafColor,endColor,nodeColor)
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.isCluster==true
                nodeIndex=bioTree{iframe}.root{iroot}.nodeIndex;
                set(treeGraph.Nodes(nodeIndex),'Color',rootColor);
                set(treeGraph.Nodes(nodeIndex),'LineColor',rootColor);
                set(treeGraph.Nodes(nodeIndex),'UserData',iframe);
                set(treeGraph.Nodes(nodeIndex),'Label','root');         
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.isCluster==true
                nodeIndex=bioTree{iframe}.leavies{ileaf}.nodeIndex;               
                if iframe < size(bioTree,2)
                    set(treeGraph.Nodes(nodeIndex),'Color',leafColor);
                    set(treeGraph.Nodes(nodeIndex),'LineColor',leafColor);
                    set(treeGraph.Nodes(nodeIndex),'UserData',iframe);
                    set(treeGraph.Nodes(nodeIndex),'Label','leaf');
                end
                if iframe == size(bioTree,2)
                    set(treeGraph.Nodes(nodeIndex),'Color',endColor);
                    set(treeGraph.Nodes(nodeIndex),'LineColor',endColor);
                    set(treeGraph.Nodes(nodeIndex),'UserData',iframe);
                    set(treeGraph.Nodes(nodeIndex),'Label','endLeaf');
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.isCluster==true
                nodeIndex=bioTree{iframe}.node{iNode}.nodeIndex;
                set(treeGraph.Nodes(nodeIndex),'Color',nodeColor);
                set(treeGraph.Nodes(nodeIndex),'LineColor',nodeColor);
                set(treeGraph.Nodes(nodeIndex),'UserData',iframe);
                set(treeGraph.Nodes(nodeIndex),'Label','Node');
            end
        end
    end
end
end
function plotClusterConnection(bioTree,treeGraph,rootColor,leafColor,endColor,nodeColor,MakerSize,lineColor)
[xMin,xMax]=finfXedgy(treeGraph);
stepScale=1200;
offSet=1;
angle=180;
figure;
xAxis=xMin:(xMax-xMin)/500:xMax;
yAxis=1:stepScale:size(bioTree,2);
for iAxis=1:size(yAxis,2)
    y=yAxis(iAxis);
    [x1,y1]=rePosition(xAxis,y,xMin,xMax,offSet,angle);
    plot(x1,y1);hold on;
end
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.isCluster==true
                nodeIndex1=bioTree{iframe}.root{iroot}.nodeIndex;
                nodeInfo=bioTree{iframe}.root{iroot}.nodeInfo;
                xy1=get(treeGraph.Nodes(nodeIndex1),'Position');
                x1=xy1(1);
                y1=iframe;
                nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                xy2=get(treeGraph.Nodes(nodeIndex2),'Position');
                x2=xy2(1);
                y2=nodeInfo(1);
                [x1,y1]=rePosition(x1,y1,xMin,xMax,offSet,angle);
                [x2,y2]=rePosition(x2,y2,xMin,xMax,offSet,angle);
                plot(x1,y1,'MarkerFaceColor',rootColor,'MarkerSize',MakerSize,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                plot([x1,x2],[y1,y2],'Color',lineColor);hold on;
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.isCluster==true
                nodeIndex1=bioTree{iframe}.leavies{ileaf}.nodeIndex;
                xy1=get(treeGraph.Nodes(nodeIndex1),'Position');
                x1=xy1(1);
                y1=iframe;
                [x1,y1]=rePosition(x1,y1,xMin,xMax,offSet,angle);
                if iframe < size(bioTree,2)
                   plot(x1,y1,'MarkerFaceColor',leafColor,'MarkerSize',MakerSize,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                end
                if iframe == size(bioTree,2)
                   plot(x1,y1,'MarkerFaceColor',endColor,'MarkerSize',MakerSize,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                end
            end
        end
    end
    
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.isCluster==true
                nodeIndex1=bioTree{iframe}.node{iNode}.nodeIndex;
                xy1=get(treeGraph.Nodes(nodeIndex1),'Position');
                x1=xy1(1);
                y1=iframe;
                [x1,y1]=rePosition(x1,y1,xMin,xMax,offSet,angle);
                plot(x1,y1,'MarkerFaceColor',nodeColor,'MarkerSize',MakerSize+2,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                        nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                        nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                        xy2=get(treeGraph.Nodes(nodeIndex2),'Position');
                        x2=xy2(1);
                        y2=nodeInfo(1);
                        [x2,y2]=rePosition(x2,y2,xMin,xMax,offSet,angle);
                        plot([x1,x2],[y1,y2],'Color',lineColor);hold on;
                    end
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                        leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                        nodeIndex2=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeIndex;
                        xy2=get(treeGraph.Nodes(nodeIndex2),'Position');
                        x2=xy2(1);
                        y2=leafInfo(1);
                        [x2,y2]=rePosition(x2,y2,xMin,xMax,offSet,angle);
                        plot([x1,x2],[y1,y2],'Color',lineColor);hold on;
                    end
                end
            end
        end
    end
end
end
function plotReductionClusterConnection(bioTree,treeGraph,rootColor,leafColor,endColor,nodeColor,MakerSize,lineColor)
[xMin,xMax]=finfXedgy(treeGraph);
stepScale=1200;
offSet=1;
angle=180;
figure;
xAxis=xMin:(xMax-xMin)/500:xMax;
yAxis=1:stepScale:size(bioTree,2);
for iAxis=1:size(yAxis,2)
    y=yAxis(iAxis);
    [x1,y1]=rePosition(xAxis,y,xMin,xMax,offSet,angle);
    plot(x1,y1);hold on;
end
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                nodeIndex1=bioTree{iframe}.root{iroot}.nodeIndex;
                nodeInfo=bioTree{iframe}.root{iroot}.nodeInfo;
                xy1=get(treeGraph.Nodes(nodeIndex1),'Position');
                x1=xy1(1);
                y1=iframe;
                nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                xy2=get(treeGraph.Nodes(nodeIndex2),'Position');
                x2=xy2(1);
                y2=nodeInfo(1);
                [x1,y1]=rePosition(x1,y1,xMin,xMax,offSet,angle);
                [x2,y2]=rePosition(x2,y2,xMin,xMax,offSet,angle);
                plot(x1,y1,'MarkerFaceColor',rootColor,'MarkerSize',MakerSize,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                plot([x1,x2],[y1,y2],'Color',lineColor);hold on;
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.is2Node==true
                nodeIndex1=bioTree{iframe}.leavies{ileaf}.nodeIndex;
                xy1=get(treeGraph.Nodes(nodeIndex1),'Position');
                x1=xy1(1);
                y1=iframe;
                [x1,y1]=rePosition(x1,y1,xMin,xMax,offSet,angle);
                if iframe < size(bioTree,2)
                   plot(x1,y1,'MarkerFaceColor',leafColor,'MarkerSize',MakerSize,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                end
                if iframe == size(bioTree,2)
                   plot(x1,y1,'MarkerFaceColor',endColor,'MarkerSize',MakerSize,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                end
            end
        end
    end
    
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.reduce==0
                nodeIndex1=bioTree{iframe}.node{iNode}.nodeIndex;
                xy1=get(treeGraph.Nodes(nodeIndex1),'Position');
                x1=xy1(1);
                y1=iframe;
                [x1,y1]=rePosition(x1,y1,xMin,xMax,offSet,angle);
                plot(x1,y1,'MarkerFaceColor',nodeColor,'MarkerSize',MakerSize+2,'Marker','o','MarkerEdgeColor','none','LineStyle','none');hold on;
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                        nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                        nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                        xy2=get(treeGraph.Nodes(nodeIndex2),'Position');
                        x2=xy2(1);
                        y2=nodeInfo(1);
                        [x2,y2]=rePosition(x2,y2,xMin,xMax,offSet,angle);
                        plot([x1,x2],[y1,y2],'Color',lineColor);hold on;
                    end
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                        for iList=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo,1)
                            nodeIndex2=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo(iList,3);
                            xy2=get(treeGraph.Nodes(nodeIndex2),'Position');
                            x2=xy2(1);
                            y2=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo(iList,1);
                            [x2,y2]=rePosition(x2,y2,xMin,xMax,offSet,angle);
                            plot([x1,x2],[y1,y2],'Color',lineColor);hold on;
                        end
                    end
                end
            end
        end
    end
end
end
function [xMin,xMax]=finfXedgy(treeGraph)
xList=[];
for i=1:size(treeGraph.Nodes,1)
    xy=get(treeGraph.Nodes(i),'Position');
    xList=[xList,xy(1)];
end
xMin=min(xList);
xMax=max(xList);
end
function  [newX,newY]=rePosition(oldX,oldY,xMin,xMax,offset,angle)
criticalFrame=11000;
radius=offset+log(criticalFrame)-log(criticalFrame-oldY);
startAngle=(180-angle)/2;
arcAngle=((angle/(xMax-xMin)).*oldX-angle/(xMax-xMin)+startAngle)*pi/180;
newX=cos(arcAngle).*radius;
newY=sin(arcAngle).*radius;
end
function [bioTree,cMatrix]=getCMatrix(bioTree,directed) %if dirreted is ture return dirted graph, otherwise retun undireted graph
[bioTree,nodeNum]=nodeIndex(bioTree);
disp(nodeNum);
cMatrix=false(nodeNum,nodeNum);
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.isCluster==true
                nodeIndex1=bioTree{iframe}.root{iroot}.nodeIndex;
                nodeInfo=bioTree{iframe}.root{iroot}.nodeInfo;
                nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                if directed==true
                    cMatrix(nodeIndex1,nodeIndex2)=1;
                end
                if directed==false
                    cMatrix(nodeIndex1,nodeIndex2)=1;
                    cMatrix(nodeIndex2,nodeIndex1)=1;
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.isCluster==true
                nodeIndex1=bioTree{iframe}.node{iNode}.nodeIndex;
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                        leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
                        nodeIndex2=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeIndex;
                        if directed==true
                            cMatrix(nodeIndex1,nodeIndex2)=1;
                        end
                        if directed==false
                            cMatrix(nodeIndex1,nodeIndex2)=1;
                            cMatrix(nodeIndex2,nodeIndex1)=1;
                        end
                    end
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                        nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                        nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                        if directed==true
                            cMatrix(nodeIndex1,nodeIndex2)=1;
                        end
                        if directed==false
                            cMatrix(nodeIndex1,nodeIndex2)=1;
                            cMatrix(nodeIndex2,nodeIndex1)=1;
                        end
                    end
                end
            end
        end
    end
end
end

function [bioTree,cMatrix]=getReductionCMatrix(bioTree,directed) %if dirreted is ture return dirted graph, otherwise retun undireted graph
[bioTree,nodeNum]=nodeIndexReduction(bioTree);
disp(nodeNum);
cMatrix=false(nodeNum,nodeNum);
for iframe=1:size(bioTree,2)    
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2) 
            if bioTree{iframe}.node{iNode}.reduce==0
                nodeIndex1=bioTree{iframe}.node{iNode}.nodeIndex;
                for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
                    if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                        for iList=1:size(bioTree{iframe}.node{iNode}.In{iIn}.rootInfo,1)
                            nodeIndex2=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo(iList,3);
                            if directed==true
                                cMatrix(nodeIndex2,nodeIndex1)=1;
                            end
                            if directed==false
                                cMatrix(nodeIndex1,nodeIndex2)=1;
                                cMatrix(nodeIndex2,nodeIndex1)=1;
                            end
                        end
                    end
                end
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                        for iList=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo,1)
                            nodeIndex2=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo(iList,3);
                            if directed==true
                                cMatrix(nodeIndex1,nodeIndex2)=1;
                            end
                            if directed==false
                                cMatrix(nodeIndex1,nodeIndex2)=1;
                                cMatrix(nodeIndex2,nodeIndex1)=1;
                            end
                        end
                    end
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                        nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                        nodeIndex2=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.nodeIndex;
                        if directed==true
                            cMatrix(nodeIndex1,nodeIndex2)=1;
                        end
                        if directed==false
                            cMatrix(nodeIndex1,nodeIndex2)=1;
                            cMatrix(nodeIndex2,nodeIndex1)=1;
                        end
                    end
                end
            end
        end
    end
end
end
function [bioTree,nodeNum]=nodeIndex(bioTree)
nodeIndex=1;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                bioTree{iframe}.root{iroot}.isCluster=true;
                bioTree{iframe}.root{iroot}.nodeIndex=nodeIndex;
                nodeIndex=nodeIndex+1;
            else
                bioTree{iframe}.root{iroot}.isCluster=false;
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.is2Node==true
                bioTree{iframe}.leavies{ileaf}.isCluster=true;
                bioTree{iframe}.leavies{ileaf}.nodeIndex=nodeIndex;
                nodeIndex=nodeIndex+1;
            else
                bioTree{iframe}.leavies{ileaf}.isCluster=false;
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            bioTree{iframe}.node{iNode}.isCluster=true;
            bioTree{iframe}.node{iNode}.nodeIndex=nodeIndex;
            nodeIndex=nodeIndex+1;
        end
    end
end
nodeNum=nodeIndex-1;
end
function [bioTree,nodeNum]=nodeIndexReduction(bioTree)
nodeIndex=1;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            if bioTree{iframe}.node{iNode}.reduce==0
                bioTree{iframe}.node{iNode}.isCluster=true;
                for iIn=1:size(bioTree{iframe}.node{iNode}.In,2)
                    if bioTree{iframe}.node{iNode}.In{iIn}.isNode==false
                        for iList=1:size(bioTree{iframe}.node{iNode}.In{iIn}.rootInfo,1)
                            rootInfo=bioTree{iframe}.node{iNode}.In{iIn}.rootInfo(iList,:);
                            bioTree{rootInfo(1)}.root{rootInfo(2)}.isCluster=true;
                            bioTree{rootInfo(1)}.root{rootInfo(2)}.nodeIndex=nodeIndex;
                            bioTree{iframe}.node{iNode}.In{iIn}.rootInfo(iList,3)=nodeIndex;
                            nodeIndex=nodeIndex+1;
                        end
                    end
                end
                bioTree{iframe}.node{iNode}.nodeIndex=nodeIndex;
                nodeIndex=nodeIndex+1;
                for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                    if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==false
                        for iList=1:size(bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo,1)
                            leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo(iList,:);
                            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.isCluster=true;
                            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeIndex=nodeIndex;
                            bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo(iList,3)=nodeIndex;
                            nodeIndex=nodeIndex+1;
                        end
                    end
                end
            end
        end
    end
end
nodeNum=nodeIndex-1;
end

