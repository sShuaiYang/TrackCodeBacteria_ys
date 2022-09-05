function longTimeResult=longTimeStatisticForTree(bioTree,dirFile,longTimeResult)
close all
% % % % dirFile=uigetdir();
dirResultSave=strcat(dirFile,'\longTimeResult');
mkdir(dirResultSave)
cd(dirResultSave)
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
[bioTree,~,allList]=divisionFinder(bioTree,branchList);
longTimeResult.attachingNum=getAttachingFigure(bioTree);
longTimeResult.detachingNum=getDetachingFigure(bioTree,longTimeResult);
[attachingToDivision,attachingToDetaching,AtoDtoDi,AtoDtoDe,attach2DetachTime]=findAttachingTransRatio(longTimeResult.attachingNum.realAttachBranch,bioTree,branchList);
[divisionToDivision,divisionToDetaching,divisionToDivision1,divisionToDetaching1,divi2DetachTime]=findDivisionTransRatio(longTimeResult.attachingNum.firstFrameBranch,bioTree,branchList);
longTimeResult.ratio.attachingToDivision=attachingToDivision;
longTimeResult.ratio.attachingToDetaching=attachingToDetaching;
longTimeResult.ratio.AtoDtoDi=AtoDtoDi;
longTimeResult.ratio.AtoDtoDe=AtoDtoDe;
longTimeResult.ratio.divisionToDivision=divisionToDivision;
longTimeResult.ratio.divisionToDetaching=divisionToDetaching;
longTimeResult.ratio.divisionToDivisionFirst=divisionToDivision1;
longTimeResult.ratio.divisionToDetachingFirst=divisionToDetaching1;
bioTreeSize=size(bioTree,2);
longTimeResult.divideFrame=getDivisionFigure(size(bioTree,2),allList,longTimeResult,bioTreeSize);
longTimeResult.bacteriaNum=getBacteriaNum(bioTree,longTimeResult);
h1=createfigure(longTimeResult.bacteriaNum.firstFrame(:,1),longTimeResult.bacteriaNum.firstFrame(:,2),...
    longTimeResult.detachingNum.firstFrame.crawlOut(:,1),longTimeResult.detachingNum.firstFrame.crawlOut(:,2),...
    longTimeResult.detachingNum.firstFrame.realDetach(:,1),longTimeResult.detachingNum.firstFrame.realDetach(:,2),...
    longTimeResult.divideFrame.firstFrame(:,1),longTimeResult.divideFrame.firstFrame(:,2),...
    longTimeResult.bacteriaNum.attachingOne(:,1),longTimeResult.bacteriaNum.attachingOne(:,2),...
    longTimeResult.detachingNum.attachingOne.crawlOut(:,1),longTimeResult.detachingNum.attachingOne.crawlOut(:,2),...
    longTimeResult.detachingNum.attachingOne.realDetach(:,1),longTimeResult.detachingNum.attachingOne.realDetach(:,2),...
    longTimeResult.divideFrame.attachingOne(:,1),longTimeResult.divideFrame.attachingOne(:,2),...
    longTimeResult.attachingNum.realAttach(:,1),longTimeResult.attachingNum.realAttach(:,2),...
    longTimeResult.bacteriaNum.crawlOne(:,1),longTimeResult.bacteriaNum.crawlOne(:,2),...
    longTimeResult.detachingNum.crawlOne.crawlOut(:,1),longTimeResult.detachingNum.crawlOne.crawlOut(:,2),...
    longTimeResult.detachingNum.crawlOne.realDetach(:,1),longTimeResult.detachingNum.crawlOne.realDetach(:,2),...
    longTimeResult.divideFrame.crawlOne(:,1),longTimeResult.divideFrame.crawlOne(:,2),...
    longTimeResult.attachingNum.crawlIn(:,1),longTimeResult.attachingNum.crawlIn(:,2),...
    [longTimeResult.ratio.attachingToDetaching,longTimeResult.ratio.attachingToDivision;longTimeResult.ratio.divisionToDivision,longTimeResult.ratio.divisionToDetaching;...
    longTimeResult.ratio.divisionToDivisionFirst,longTimeResult.ratio.divisionToDetachingFirst;longTimeResult.ratio.AtoDtoDi,longTimeResult.ratio.AtoDtoDe]);
saveas(h1,'communityProfile.fig');
pResult1=easyHist(attach2DetachTime);
pResult2=easyHist(divi2DetachTime);
h2=createRetentionTimeBeforeDetaching(pResult1(:,1),pResult1(:,2),pResult2(:,1),pResult2(:,2));
saveas(h2,'detachingRetentionDistribution.fig')
h3=getRetentionTime(branchList,bioTree,bioTreeSize);
saveas(h3,'retentionTime.fig')
bacteriaTime=gainLengthInfo(bioTree,branchList,2);
getbacteriaTimeGraph(bacteriaTime,dirFile);
save('longTimeResult.mat','longTimeResult');
close all
end
function attachingNum=getAttachingFigure(bioTree)
attachingNum.crawlIn=[];
attachingNum.realAttach=[];
attachingNum.crawlInBranch=[];
attachingNum.realAttachBranch=[];
attachingNum.firstFrameBranch=[];
crawlRoot=0;
attachRoot=0;
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
for iRoot=1:size(bioTree{1}.root,2)
    attachingNum.firstFrameBranch=[attachingNum.firstFrameBranch;bioTree{1}.root{iRoot}.branchIndex];
end
for iframe=2:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iRoot=1:size(bioTree{iframe}.root,2)
            pixelIdxList=bioTree{iframe}.root{iRoot}.rootPixelDetail;
            if any(ismember(limitRegionIdx,pixelIdxList))==1
                crawlRoot=crawlRoot+1;
                attachingNum.crawlInBranch=[attachingNum.crawlInBranch;bioTree{iframe}.root{iRoot}.branchIndex];
            else
                attachRoot=attachRoot+1;
                attachingNum.realAttachBranch=[attachingNum.realAttachBranch;bioTree{iframe}.root{iRoot}.branchIndex];
            end
        end
        attachingNum.crawlIn=[attachingNum.crawlIn;[iframe,crawlRoot]];
        attachingNum.realAttach=[attachingNum.realAttach;[iframe,attachRoot]];
    end
end
% figure;plot(attachingNum.crawlIn(:,1),attachingNum.crawlIn(:,2),'r');
% hold on
% plot(attachingNum.realAttach(:,1),attachingNum.realAttach(:,2),'c');
end
function detachingNum=getDetachingFigure(bioTree,longTimeResult)
detachingNum.firstFrame=getDetachingCase(bioTree,longTimeResult.attachingNum.firstFrameBranch);
detachingNum.crawlOne=getDetachingCase(bioTree,longTimeResult.attachingNum.crawlInBranch);
detachingNum.attachingOne=getDetachingCase(bioTree,longTimeResult.attachingNum.realAttachBranch);
% figure;plot(detachingNum.crawlOut(:,1),detachingNum.crawlOut(:,2),'r');
% hold on
% plot(detachingNum.realDetach(:,1),detachingNum.realDetach(:,2),'c');
end
function detachingBranch=getDetachingCase(bioTree,properBranch)
detachingBranch.crawlOut=[];
detachingBranch.realDetach=[];
crawlLeaf=0;
detachLeaf=0;
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
for iframe=1:size(bioTree,2)-1
    if ~isempty(bioTree{iframe}.leavies)
        for iLeaf=1:size(bioTree{iframe}.leavies,2)
            pixelIdxList=bioTree{iframe}.leavies{iLeaf}.leaviesPixelDetail;
            branchIndex=bioTree{iframe}.leavies{iLeaf}.branchIndex;
            if any(ismember(limitRegionIdx,pixelIdxList))==1 && ismember(branchIndex,properBranch)
                crawlLeaf=crawlLeaf+1;
            end
            if ~(any(ismember(limitRegionIdx,pixelIdxList))==1) && ismember(branchIndex,properBranch)
                detachLeaf=detachLeaf+1;
            end
        end
        detachingBranch.crawlOut=[detachingBranch.crawlOut;[iframe,crawlLeaf]];
        detachingBranch.realDetach=[detachingBranch.realDetach;[iframe,detachLeaf]];
    end
end
end
function divideFrame=getDivisionFigure(frameNum,allList,longTimeResult,bioTreeSize)
nodeList=allList.allNode(:,[1,4,6]);
nodeList(nodeList(:,2)==0,:)=[];
nodeListInfo=nodeList(:,1);
branchInfo=nodeList(:,3);
divideFrame.firstFrame=zeros(bioTreeSize,2);
divideFrame.firstFrame(:,1)=1:bioTreeSize;
divideFrame.crawlOne=zeros(bioTreeSize,2);
divideFrame.crawlOne(:,1)=1:bioTreeSize;
divideFrame.attachingOne=zeros(bioTreeSize,2);
divideFrame.attachingOne(:,1)=1:bioTreeSize;
for iframe=1:frameNum
    divideFrame.firstFrame(iframe,2)=numel(nodeListInfo(nodeListInfo<=iframe & ismember(branchInfo,longTimeResult.attachingNum.firstFrameBranch)));
    divideFrame.crawlOne(iframe,2)=numel(nodeListInfo(nodeListInfo<=iframe & ismember(branchInfo,longTimeResult.attachingNum.crawlInBranch)));
    divideFrame.attachingOne(iframe,2)=numel(nodeListInfo(nodeListInfo<=iframe & ismember(branchInfo,longTimeResult.attachingNum.realAttachBranch)));
end
% figure;plot(divideFrame)
end
function bacteriaNum=getBacteriaNum(bioTree,longTimeResult)
bacteriaNum.firstFrame=zeros(size(bioTree,2),2);
bacteriaNum.firstFrame(:,1)=1:size(bioTree,2);
bacteriaNum.crawlOne=zeros(size(bioTree,2),2);
bacteriaNum.crawlOne(:,1)=1:size(bioTree,2);
bacteriaNum.attachingOne=zeros(size(bioTree,2),2);
bacteriaNum.attachingOne(:,1)=1:size(bioTree,2);
for iframe=1:size(bioTree,2)
    for smallFrame=1:iframe
        if ~isempty(bioTree{smallFrame}.root)
            for iRoot=1:size(bioTree{smallFrame}.root,2)
                traceInfo=bioTree{smallFrame}.root{iRoot}.traceInfo.pixelIdxList;
                if size(traceInfo,2)+smallFrame-1>=iframe
                    branchIndex=bioTree{smallFrame}.root{iRoot}.branchIndex;
                    if ismember(branchIndex,longTimeResult.attachingNum.firstFrameBranch)
                        bacteriaNum.firstFrame(iframe,2)=bacteriaNum.firstFrame(iframe,2)+1;
                    end
                    if ismember(branchIndex,longTimeResult.attachingNum.crawlInBranch)
                        bacteriaNum.crawlOne(iframe,2)=bacteriaNum.crawlOne(iframe,2)+1;
                    end
                    if ismember(branchIndex,longTimeResult.attachingNum.realAttachBranch)
                        bacteriaNum.attachingOne(iframe,2)=bacteriaNum.attachingOne(iframe,2)+1;
                    end
                end
            end
        end
        if ~isempty(bioTree{smallFrame}.node)
            for iNode=1:size(bioTree{smallFrame}.node,2)
                for iOut=1:size(bioTree{smallFrame}.node{iNode}.Out,2);
                    traceInfo=bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList;
                    if size(traceInfo,2)+smallFrame-1>=iframe
                        branchIndex=bioTree{smallFrame}.node{iNode}.branchIndex;
                        if ismember(branchIndex,longTimeResult.attachingNum.firstFrameBranch)
                            bacteriaNum.firstFrame(iframe,2)=bacteriaNum.firstFrame(iframe,2)+1;
                        end
                        if ismember(branchIndex,longTimeResult.attachingNum.crawlInBranch)
                            bacteriaNum.crawlOne(iframe,2)=bacteriaNum.crawlOne(iframe,2)+1;
                        end
                        if ismember(branchIndex,longTimeResult.attachingNum.realAttachBranch)
                            bacteriaNum.attachingOne(iframe,2)=bacteriaNum.attachingOne(iframe,2)+1;
                        end
                    end
                end
            end
        end
    end
end
end
function limitRegionIdx=getLimitRegion(cropInfo)
cropInfo1=bwmorph(cropInfo,'remove');
cropInfo1=imdilate(cropInfo1,true(13));
limitRegion=cropInfo1 & cropInfo;
cc=bwconncomp(limitRegion);
limitRegionIdx=cc.PixelIdxList{1,1};
end
function [attachingToDivision,attachingToDetaching,AtoDtoDi,AtoDtoDe,attach2DetachTime]=findAttachingTransRatio(attachingBranch,bioTree,branchList)
toDivision=0;
toDetaching=0;
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
usefulAttaching=0;
attachingToDivisionBranch=[];
attach2DetachTime=[];
for iBranch=1:size(branchList,1)
    if ismember(iBranch,attachingBranch)
        branchInfo=branchList(iBranch,:);
        if branchInfo(3)==0
            rootInfo=branchInfo([1,2]);
            leafInfo=bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo;
            if leafInfo(1)~=size(bioTree,2)
                pixelIdxList=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail;
                usefulAttaching=usefulAttaching+1;
                if ~(any(ismember(limitRegionIdx,pixelIdxList))==1)
                    attach2DetachTime=[attach2DetachTime;size(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList,2)];
                    toDetaching=toDetaching+1;
                end
            end
        end
        if branchInfo(3)==1
            nodeInfo=branchInfo([1,2]);
            allNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allNode;
            usefulAttaching=usefulAttaching+1;
            if allNode(1,4)==1
                toDivision=toDivision+1;
                attachingToDivisionBranch=[attachingToDivisionBranch;iBranch];
            end
        end
    end
end
attachingToDivision=toDivision/usefulAttaching;
attachingToDetaching=toDetaching/usefulAttaching;
toDivision=0;
toDetaching=0;
divisionNode=0;
for iBranch=1:size(branchList,1)
    if ismember(iBranch,attachingToDivisionBranch)
        branchInfo=branchList(iBranch,:);
        if branchInfo(3)==1
            nodeInfo=branchInfo([1,2]);
            allNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allNode;
            allNode(1,:)=[];
            for iNode=1:size(allNode,1)
                smallNodeInfo=allNode(iNode,:);
                if smallNodeInfo(:,4)==1
                    for iOut=1:2
                        divisionNode=divisionNode+1;
                        if bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.is2Node==0;
                            leafInfo=bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.leafInfo;
                            if leafInfo(1)~=size(bioTree,2)
                                pixelIdxList=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail;
                                if ~(any(ismember(limitRegionIdx,pixelIdxList))==1)
                                    toDetaching=toDetaching+1;
                                else
                                    if numel(bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList)>=600
                                        toDivision=toDivision+1;
                                    end
                                end
                            else
                                divisionNode=divisionNode-1;
                            end
                        end
                        if bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.is2Node==1;
                            nextNode=bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.nodeInfo;
                            if allNode(allNode(:,1)==nextNode(1) & allNode(:,2)==nextNode(2),4)==1;
                                toDivision=toDivision+1;
                            end
                        end
                    end
                end
            end
        end
    end
end
AtoDtoDi=toDivision/divisionNode;
AtoDtoDe=toDetaching/divisionNode;
end
function [divisionToDivision,divisionToDetaching,divisionToDivision1,divisionToDetaching1,divi2detachTime]=findDivisionTransRatio(firstBranch,bioTree,branchList)
toDivision=0;
toDetaching=0;
divisionNode=0;
toDivision1=0;
toDetaching1=0;
divisionNode1=0;
limitRegionIdx=getLimitRegion(bioTree{1}.imageProcessingInfo.cropInfo);
divi2detachTime=[];
for iBranch=1:size(branchList,1)
    branchInfo=branchList(iBranch,:);
    if branchInfo(3)==1
        nodeInfo=branchInfo([1,2]);
        allNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allNode;
        for iNode=1:size(allNode,1)
            smallNodeInfo=allNode(iNode,:);
            if smallNodeInfo(:,4)==1
                for iOut=1:2
                    divisionNode=divisionNode+1;
                    if bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.is2Node==0;
                        leafInfo=bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.leafInfo;
                        if leafInfo(1)~=size(bioTree,2)
                            pixelIdxList=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail;
                            if ~(any(ismember(limitRegionIdx,pixelIdxList))==1)
                                toDetaching=toDetaching+1;
                                divi2detachTime=[divi2detachTime;size(bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList,2)];
                            else
                                if numel(bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList)>=600
                                    toDivision=toDivision+1;
                                end
                            end
                        else
                            divisionNode=divisionNode-1;
                        end
                    end
                    if bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.is2Node==1;
                        nextNode=bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.nodeInfo;
                        if allNode(allNode(:,1)==nextNode(1) & allNode(:,2)==nextNode(2),4)==1;
                            toDivision=toDivision+1;
                        end
                    end
                end
            end
        end
    end
end
for iBranch=1:size(branchList,1)
    if ismember(iBranch,firstBranch)
        branchInfo=branchList(iBranch,:);
        if branchInfo(3)==1
            nodeInfo=branchInfo([1,2]);
            allNode=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allNode;
            for iNode=1:size(allNode,1)
                smallNodeInfo=allNode(iNode,:);
                if smallNodeInfo(:,4)==1
                    for iOut=1:2
                        divisionNode1=divisionNode1+1;
                        if bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.is2Node==0;
                            leafInfo=bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.leafInfo;
                            if leafInfo(1)~=size(bioTree,2)
                                pixelIdxList=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.leaviesPixelDetail;
                                if ~(any(ismember(limitRegionIdx,pixelIdxList))==1)
                                    toDetaching1=toDetaching1+1;
                                else
                                    if numel(bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.traceInfo.pixelIdxList)>=600
                                        toDivision1=toDivision1+1;
                                    end
                                end
                            else
                                divisionNode1=divisionNode1-1;
                            end
                        end
                        if bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.is2Node==1;
                            nextNode=bioTree{smallNodeInfo(1)}.node{smallNodeInfo(2)}.Out{iOut}.nodeInfo;
                            if allNode(allNode(:,1)==nextNode(1) & allNode(:,2)==nextNode(2),4)==1;
                                toDivision1=toDivision1+1;
                            end
                        end
                    end
                end
            end
        end
    end
end
divisionToDivision=toDivision/divisionNode;
divisionToDetaching=toDetaching/divisionNode;
divisionToDivision1=toDivision1/divisionNode1;
divisionToDetaching1=toDetaching1/divisionNode1;
end
function h=createfigure(a1,a2,b1,b2,c1,c2,d1,d2,a3,a4,b3,b4,c3,c4,d3,d4,e3,e4,a5,a6,b5,b6,c5,c6,d5,d6,e5,e6,ymatrix)
xvector=[1,2,3,4];
%CREATEFIGURE(LONGTIMERESULT1,YMATRIX1,LONGTIMERESULT2,YMATRIX2)
%  LONGTIMERESULT1:  vector of x data
%  YMATRIX1:  matrix of y data
%  LONGTIMERESULT2:  vector of x data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 19-Oct-2012 21:53:36

% Create figure
figure1 = figure;
scrsz=get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);

% Create axes1
axes1 = axes('Parent',figure1,...
    'Position',[0.550520833333333 0.557569296375267 0.283854166666667 0.390951877841952]);
box(axes1,'on');
hold on

% Create title
title('first frame bacteria fate','FontSize',20);

% Create xlabel
xlabel('time/frame','FontSize',20);

% Create ylabel
ylabel('bacteriaNum','FontSize',20);

% Create multiple lines using matrix input to plot
plot1 = plot(a1,a2,'Parent',axes1);
set(plot1(1),'LineWidth',3,'Color',[1 0 0],'DisplayName','bacteriaNum');

plot2 = plot(b1,b2,'Parent',axes1);
set(plot2(1),'LineWidth',3,'DisplayName','crawlOut');

% Create multiple lines using matrix input to plot
plot3 = plot(c1,c2,'Parent',axes1);
set(plot3(1),'LineWidth',3,'Color',[0.87058824300766 0.490196079015732 0],...
    'DisplayName','realDetach');

plot4 = plot(d1,d2,'Parent',axes1);
set(plot4(1),'LineWidth',3,...
    'Color',[0 1 1],...
    'DisplayName','divisionNum');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.567317708333335 0.833448619658259 0.0755208333333333 0.109425287356322]);

% Create axes2
axes2 = axes('Parent',figure1,...
    'Position',[0.548958333333333 0.0469722814498909 0.287499999999999 0.408251599147124]);
box(axes2,'on');
hold on

% Create title
title('attaching bacteria fate','FontSize',20);

% Create xlabel
xlabel('time/frame','FontSize',20);

% Create ylabel
ylabel('bacteriaNum','FontSize',20);

% Create multiple lines using matrix input to plot
plot1 = plot(a3,a4,'Parent',axes2);
set(plot1(1),'LineWidth',3,'Color',[1 0 0],'DisplayName','bacteriaNum');

plot2 = plot(b3,b4,'Parent',axes2);
set(plot2(1),'LineWidth',3,'DisplayName','crawlOut');

% Create multiple lines using matrix input to plot
plot3 = plot(c3,c4,'Parent',axes2);
set(plot3(1),'LineWidth',3,'Color',[0.87058824300766 0.490196079015732 0],...
    'DisplayName','realDetach');

plot4 = plot(d3,d4,'Parent',axes2);
set(plot4(1),'LineWidth',3,...
    'Color',[0 1 1],...
    'DisplayName','divisionNum');

plot5 = plot(e3,e4,'Parent',axes2);
set(plot5(1),'LineWidth',3,...
    'Color',[0 1 0],'DisplayName','realAttach');

% Create legend
legend2 = legend(axes2,'show');
set(legend2,...
    'Position',[0.562109375000002 0.334494522461581 0.0755208333333333 0.109425287356322]);

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.1784375 0.0426439232409382 0.293958333333333 0.424307036247335]);
box(axes3,'on');
hold on

% Create title
title('crawling bacteria fate','FontSize',20);

% Create xlabel
xlabel('time/frame','FontSize',20);

% Create ylabel
ylabel('bacteriaNum','FontSize',20);

% Create multiple lines using matrix input to plot
plot1 = plot(a5,a6,'Parent',axes3);
set(plot1(1),'Color',[1 0 0],'DisplayName','bacteriaNum','LineWidth',3);

plot2 = plot(b5,b6,'Parent',axes3);
set(plot2(1),'LineWidth',3,'DisplayName','crawlOut');

% Create multiple lines using matrix input to plot
plot3 = plot(c5,c6,'Parent',axes3);
set(plot3(1),'Color',[0.87058824300766 0.490196079015732 0],'LineWidth',3,...
    'DisplayName','realDetach');

plot4 = plot(d5,d6,'Parent',axes3);
set(plot4(1),'LineWidth',3,...
    'Color',[0 1 1],...
    'DisplayName','divisionNum');

plot5 = plot(e5,e6,'Parent',axes3);
set(plot5(1),'LineWidth',3,...
    'Color',[0 1 0],'DisplayName','crawlIn');

% Create legend
legend3 = legend(axes3,'show');
set(legend3,...
    'Position',[0.192838541666669 0.338034824116492 0.0755208333333332 0.113005664856733]);

% Create axes
axes4 = axes('Parent',figure1,'XTick',zeros(1,0),...
    'Position',[0.189583333333333 0.551365795221191 0.283333333333333 0.390951877841952]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes4,[0 1]);
box(axes4,'on');
hold(axes4,'all');

% Create multiple lines using matrix input to bar
bar(xvector,ymatrix,'FaceColor',[0 0 1],'Parent',axes4);

% Create textbox
% Create textbox
annotation(figure1,'textbox',...
    [0.203604166666667 0.513176144244105 0.0208750000000001 0.0319001386962552],...
    'String',{'A-De'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.224645833333333 0.513183041078701 0.0208750000000001 0.0319001386962552],...
    'String',{'A-Di'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.273291666666666 0.513987479358868 0.0208750000000001 0.0319001386962552],...
    'String',{'Di-Di'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.342770833333333 0.514791917639034 0.0208750000000001 0.0319001386962552],...
    'String',{'1Di-Di'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.293083333333333 0.513987479358868 0.0208750000000001 0.0319001386962552],...
    'String',{'Di-De'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.364125 0.514797994145286 0.0208750000000001 0.0319001386962552],...
    'String',{'1Di-De'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.413083333333333 0.514791917639034 0.0208750000000001 0.0319001386962552],...
    'String',{'ADi-Di'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.434958333333333 0.514797994145286 0.0208750000000001 0.0319001386962552],...
    'String',{'ADi-De'},...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.206208333333334 0.913845201802262 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(1,1))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.228291666666667 0.581687480005985 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(1,2))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.277458333333333 0.798270174420004 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(2,1))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.297458333333333 0.664359988220661 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(2,2))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.3469375 0.825367655252424 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(3,1))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.3688125 0.653407085701494 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(3,2))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.4188125 0.801848537804193 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(4,1))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.439125 0.656217469843028 0.0167083333333335 0.0167890870933898],...
    'String',{num2str(ymatrix(4,2))},...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none');
h=gcf;
end
function h=getRetentionTime(branchList,bioTree,bioTreeSize)
imageSize=bioTree{1}.imageSize;
limitRegionIdx=getLimitRegion(imageSize);
nodeOrRoot=branchList(:,3);
coreBranch=numel(nodeOrRoot(nodeOrRoot==1));
retentionTime=[];
for i=1:size(branchList,1)-coreBranch
    rootInfo=branchList(coreBranch+i,[1,2]);
    retentionTime=[retentionTime;numel(bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList)];
end
for i=1:coreBranch
    nodeInfo=branchList(i,[1,2]);
    allRoot=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allRoot;
    if size(allRoot,1)==1
        allLeaf=bioTree{nodeInfo(1)}.node{nodeInfo(2)}.allLeaf;
        if isempty(allLeaf)
            continue
        end
        maxFrame=max(allLeaf(:,1));
        retentionTime=[retentionTime;maxFrame-allRoot(1)+1];
        for iLeaf=1:size(allLeaf,1)
            leafInfo=allLeaf(iLeaf,:);
            if leafInfo(1)~=bioTreeSize;
                is2Node=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.is2Node;
                if is2Node==0
                    rootInfo=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo;
                    retentionTime=[retentionTime;leafInfo(1)-rootInfo(1)+1];
                end
                if is2Node==1
                    nodeInfo=bioTree{leafInfo(1)}.leavies{leafInfo(2)}.nodeInfo;
                    retentionTime=[retentionTime;leafInfo(1)-nodeInfo(1)+1];
                end
            end
        end
    end
end
binSpace=linspace(0,size(bioTree,2),31);
[countP1,n1]=get1Dhist(retentionTime,binSpace);
figure;h=plot(n1,countP1);
end
function [count n]=get1Dhist(velocityData,binSpace)
nBin=binSpace;
count=[];
n=[];
for i=2:2:size(nBin,2)
    count=[count,numel(velocityData(velocityData>=nBin(i-1)&velocityData<nBin(i+1)))];
    n=[n,nBin(i)];
end
end
function pResult=easyHist(originalData)
p=linspace(1,max(originalData),121);
pResult=[];
for iNum=2:2:numel(p)-1
    pHist=0;
    for i=1:numel(originalData)
        if originalData(i)>=p(iNum-1) && originalData(i)<p(iNum+1)
            pHist=pHist+1;
        end
    end
    pResult=[pResult;p(iNum),pHist];
end
pResult(:,2)=pResult(:,2)/sum(pResult(:,2));
plot(pResult(:,1),pResult(:,2))
end
function h=createRetentionTimeBeforeDetaching(X1, Y1, X2, Y2)
%CREATEFIGURE(X1,Y1,X2,Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 30-Oct-2012 20:21:13

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on',......
    'Position',[0.133704735376045 0.221821086261979 0.762256267409471 0.54814696485623]);
box(axes1,'on');
hold(axes1,'all');

% Create loglog
loglog(X1,Y1,'Parent',axes1,'Marker','o','LineWidth',3,'Color',[1 0 0],...
    'DisplayName','atta2detach');

% Create loglog
loglog(X2,Y2,'Parent',axes1,'Marker','square','LineWidth',3,...
    'DisplayName','divi2detach');

% Create xlabel
xlabel('retention Time','FontSize',20);

% Create ylabel
ylabel('ratio','FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.647131789693594 0.650044060166837 0.175487465181058 0.0670926517571885]);
h=gca;
end