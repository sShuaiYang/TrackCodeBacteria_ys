function treeModelWithPlot3(bioTree,cutNum)
% 用plot3函数画出不同branch的bioTree的程序
% orderBranchList=[3,4,6];
orderBranchList=[4];
n=0;
xyDiv=3;
zDiv=40;
traceInterval=5;
traceEmptyInterval=0;
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,branchList,allList]=divisionFinder(bioTree,branchList);
% add differentNode;
nodeWidth=[10,9,8,7,6,5,4,3,2,1,1,1,1,1,1];
rotateNum=0;
% for cutNum=[30:50:14930,14930*ones(1,60)]
for cutNum=[14930*ones(1,20)]
    rotateNum=rotateNum+1;
    figure1 = figure('InvertHardcopy','off','Color',[0 0 0]);
    scrsz=get(0,'ScreenSize');
    set (gcf,'Position',[10 10 scrsz(3) scrsz(4)-80]);
    set(gcf, 'PaperPositionMode', 'auto');
%     % axes on
%     axes1 = axes('Parent',figure1,'ZColor',[0 0 1],'YColor',[0 0 1],...
%         'XColor',[0 0 1],'LineWidth',3,'CameraTarget',[100,75,6.5],...
%         'Position',[0.263020833333333 0.122461059190031 0.5175 0.818348909657321],...
%         'FontSize',14,...
%         'Color',[0 0 0]);
%     view(axes1,[-143+rotateNum*1 14]);
%     box(axes1,'on');
    axes1 = axes('Parent',figure1,'ZColor',[0 0 0],'YColor',[0 0 0],...
        'XColor',[0 0 0],'CameraTarget',[100,75,6.5],'CameraPosition',[-934.7,106.82,33.7],'DataAspectRatio',[1,1,0.1],...
        'Color',[0 0 0]);
    xlim(axes1,[0 200]);
    ylim(axes1,[0 150]);
    zlim(axes1,[0 13]);
%     view(axes1,[-143+rotateNum*1 14]);
    view(axes1,[-143 14+rotateNum*7.6]);
    if rotateNum>=11
        view(axes1,[-143+(rotateNum-10)*14.3 90]);
    end
    set(gca,'CameraTarget',[100,75,6.5]);
    box(axes1,'on');
    
    hold(axes1,'all')
%     xlabel('x(μm)','HorizontalAlignment','right','FontSize',14,'Color',[0 0 1]);
%     ylabel('y(μm)','FontSize',14,'Color',[0 0 1],'HorizontalAlignment','left');
%     zlabel('time(h)','HorizontalAlignment','right','FontSize',14,'Color',[0 0 1]);
    bNum=0;
    for iBranch=1:numel(orderBranchList)
        bNum=bNum+1;
        focusBranch=orderBranchList(iBranch);
        aimNode=branchList(focusBranch,:);
        bioTree=initializeNodeGeneration(bioTree,aimNode);  % give order for each generation
        allRoot=bioTree{aimNode(1)}.node{aimNode(2)}.allRoot;
        centroidInfo=bioTree{allRoot(1)}.root{allRoot(2)}.rootMeasurment.Centroid;
        plot3(centroidInfo(1)*0.064,centroidInfo(2)*0.064,allRoot(1)/1200,'Marker','o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none','MarkerSize',8);
        rootTrace=bioTree{allRoot(1)}.root{allRoot(2)}.traceInfo.measurment;
        if numel(rootTrace)+allRoot(1)-1>cutNum
            rootTrace=rootTrace(1:cutNum-allRoot(1)+1);
        end
        rootTraceCentroid=[];
        for iTrace=1+traceEmptyInterval:traceInterval:numel(rootTrace)-traceEmptyInterval
            rootTraceCentroid=[rootTraceCentroid;rootTrace{iTrace}.Centroid];
        end
        if allRoot(1)<=cutNum
            plot3(rootTraceCentroid(:,1)*0.064,rootTraceCentroid(:,2)*0.064,(allRoot(1)+traceEmptyInterval:traceInterval:allRoot(1)+numel(rootTrace)-1-traceEmptyInterval)/1200,'Color',[211-(bNum-1)*20,85,211]/255,'lineWidth',2);
        end
        allLeaf=bioTree{aimNode(1)}.node{aimNode(2)}.allLeaf;
        for iLeaf=1:size(allLeaf,1)
            eachLeaf=allLeaf(iLeaf,:);
            if eachLeaf(1)<cutNum
                centroidInfo=bioTree{eachLeaf(1)}.leavies{eachLeaf(2)}.leafMeasurment.Centroid;
                plot3(centroidInfo(1)*0.064,centroidInfo(2)*0.064,eachLeaf(1)/1200,'Marker','o','MarkerFaceColor',[127,255,0]/255,'MarkerEdgeColor','none','MarkerSize',8);
            end
        end
        allNode=bioTree{aimNode(1)}.node{aimNode(2)}.allNode;
        for iNode=1:size(allNode,1)
            eachNode=allNode(iNode,:);
            if eachNode(1)<cutNum
                if bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.isNode==1
                    preNode=bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.nodeInfo;
                    centroidInfo=bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}.Centroid;
                end
                if bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.isNode==0
                    preRoot=bioTree{eachNode(1)}.node{eachNode(2)}.In{1}.rootInfo;
                    centroidInfo=bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}.Centroid;
                end
                plot3(centroidInfo(1)*0.064,centroidInfo(2)*0.064,(eachNode(1)-1)/1200,'Marker','o','MarkerFaceColor',[255,185,15]/255,'MarkerEdgeColor','none','MarkerSize',8);
                for iOut=1:size(bioTree{eachNode(1)}.node{eachNode(2)}.Out,2)
                    nodeTrace=bioTree{eachNode(1)}.node{eachNode(2)}.Out{iOut}.traceInfo.measurment;
                    nodeTraceCentroid=[];
                    if numel(nodeTrace)+eachNode(1)-1>cutNum
                        nodeTrace=nodeTrace(1:cutNum-eachNode(1)+1);
                    end
                    for iTrace=1+traceEmptyInterval:traceInterval:numel(nodeTrace)-traceEmptyInterval
                        nodeTraceCentroid=[nodeTraceCentroid;nodeTrace{iTrace}(1).Centroid];
                    end
                    %             plot3([nodeTraceCentroid(:,1)]*0.064,[nodeTraceCentroid(:,2)]*0.064,([eachNode(1)+traceEmptyInterval:traceInterval:eachNode(1)+numel(nodeTrace)-1-traceEmptyInterval])/1200,'Color',[186,85,211-10*eachNode(6)]/255,'lineWidth',2);
                    plot3([centroidInfo(1);nodeTraceCentroid(:,1)]*0.064,[centroidInfo(2);nodeTraceCentroid(:,2)]*0.064,([eachNode(1)-1,eachNode(1)+traceEmptyInterval:traceInterval:eachNode(1)+numel(nodeTrace)-1-traceEmptyInterval])/1200,'Color',[211-(bNum-1)*20,85,211-10*eachNode(6)]/255,'lineWidth',2);
                end
            end
        end
    end
%     saveas(gcf,strcat('C:\Users\jzy\Desktop\to Fan Jin\F1\tree\',num2str(rotateNum),'.tif'));
    close all
end
end
function bioTree=initializeNodeGeneration(bioTree,aimNode)
allNode=bioTree{aimNode(1)}.node{aimNode(2)}.allNode;
[~,sortOrder]=sort(allNode(:,1));
allNode=allNode(sortOrder,:);
allNode(:,6)=0;
allNode(1,6)=2;
while any(allNode(:,6)==0)
    for iNum=2:size(allNode,1)
        if allNode(iNum,6)==0
            preNode=bioTree{allNode(iNum,1)}.node{allNode(iNum,2)}.In{1}.nodeInfo;
            preNodeOrder=(allNode(:,1)==preNode(1) & allNode(:,2)==preNode(2)); 
            allNode(iNum,6)=allNode(preNodeOrder,6)+1;
        end
    end
end
bioTree{aimNode(1)}.node{aimNode(2)}.allNode=allNode;
end