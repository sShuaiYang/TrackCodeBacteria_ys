[linkMatrix,centroidInfo,leafList,leafNum]=generateOneBranchPhytree(bioTree,10,branchList(10,:));
plotLinkTree(linkMatrix,centroidInfo,1:leafNum,[32,178,170]/255);
lNum=leafNum;
[linkMatrix,centroidInfo,leafList,leafNum]=generateOneBranchPhytree(bioTree,14,branchList(14,:));
hold on;plotLinkTree(linkMatrix,centroidInfo,lNum:lNum+leafNum-1,[153,50,204]/255);
view(90,-90)
xlim(gca,[0 12460/1200]);
set(gca,'FontWeight','demi','FontSize',12,'YTick',zeros(1,0),'Position',[0.324479166666667,0.382139148494289,0.380520833333333,0.385254413291796]);
box(gca,'on');
hold(gca,'all');
% Create xlabel
xlabel('time/h','FontWeight','demi','FontSize',20);

[treeModel]=generatePDBsphere(oriStructure,bacteriaFrameInfo,bioTree);
h=molviewer(treeModel);
 evalrasmolscript(h,'select mercury;color atoms OPAQUE[0,255,0];spacefill 749')
evalrasmolscript(h,'select copper;color atoms OPAQUE[255,0,0];spacefill 749')
evalrasmolscript(h,'select silver;color atoms OPAQUE[255,255,0];spacefill 749')
evalrasmolscript(h,'select helium;color atoms OPAQUE[32,178,170];spacefill 300')
evalrasmolscript(h,'select lithium;color atoms OPAQUE[153,50,204];spacefill 300')