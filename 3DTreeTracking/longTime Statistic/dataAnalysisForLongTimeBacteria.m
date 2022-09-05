function allResult=dataAnalysisForLongTimeBacteria(bioTree)
% 慢扫分析时所用的程序  暂时多出四张图
[bioTree,~,~,~]=myBiograph_new2(bioTree);
allResult=[];

% 得到分裂之前和分裂之后的长度的数据
[allResult,bioTree]=getLengthInfoBeforeAndAfterDivision(bioTree,allResult);

% 得到分裂时间的具体数据
allResult=getDivisionTimeInfo(bioTree,allResult);

% detaching rate随时间的分布
allResult=getDetachingRateInfo(bioTree,allResult);

% 画四张图
allResultSingleMapDrawing(allResult);
end
function [result,bioTree]=getLengthInfoBeforeAndAfterDivision(bioTree,result)
result.afterDivision=[];
result.beforeDivision=[];
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        nodeInfo=bioTree{iframe}.node{iNode};
        isDivisionNode=judgeDivisionNode(nodeInfo);
        if isDivisionNode==1
            bioTree{iframe}.node{iNode}.isDivision=1;
            result.afterDivision=[result.afterDivision;nodeInfo.Out{1}.traceInfo.measurment{1}.MajorAxisLength*bioTree{1}.scaleInfo.scaleInfo];
            result.afterDivision=[result.afterDivision;nodeInfo.Out{2}.traceInfo.measurment{1}.MajorAxisLength*bioTree{1}.scaleInfo.scaleInfo];
            isNode=nodeInfo.In{1}.isNode;
            if isNode==0
                preRoot=nodeInfo.In{1}.rootInfo;
                result.beforeDivision=[result.beforeDivision;bioTree{preRoot(1)}.root{preRoot(2)}.traceInfo.measurment{end}(1).MajorAxisLength*bioTree{1}.scaleInfo.scaleInfo];
            end
            if isNode==1
                preNode=nodeInfo.In{1}.nodeInfo;
                result.beforeDivision=[result.beforeDivision;bioTree{preNode(1)}.node{preNode(2)}.Out{preNode(3)}.traceInfo.measurment{end}(1).MajorAxisLength*bioTree{1}.scaleInfo.scaleInfo];
            end
        else
            bioTree{iframe}.node{iNode}.isDivision=0;
        end
    end
end
end
function isDivisionNode=judgeDivisionNode(nodeInfo)
if ~(size(nodeInfo.In,2)==1 && size(nodeInfo.Out,2)==2)
    isDivisionNode=0;
    return
end
if numel(nodeInfo.Out{1}.traceInfo.measurment{1})>=2 || numel(nodeInfo.Out{2}.traceInfo.measurment{1})>=2
    isDivisionNode=0;
    return
end
if nodeInfo.Out{1}.traceInfo.measurment{1}.MajorAxisLength/nodeInfo.Out{2}.traceInfo.measurment{1}.MajorAxisLength<=(10/7) && nodeInfo.Out{1}.traceInfo.measurment{1}.MajorAxisLength/nodeInfo.Out{2}.traceInfo.measurment{1}.MajorAxisLength>=0.7
    isDivisionNode=1;
else
    isDivisionNode=0;
end
end
function allResult=getDivisionTimeInfo(bioTree,allResult)
allResult.divisionTime=[];
for iRoot=1:size(bioTree{1}.root,2)
    is2Node=bioTree{1}.root{iRoot}.is2Node;
    if is2Node==1
        nodeInfo=bioTree{1}.root{iRoot}.nodeInfo;
        isDivisionNode=judgeDivisionNode(bioTree{nodeInfo(1)}.node{nodeInfo(2)});
        if isDivisionNode==1
            divisionTime=nodeInfo(1)-1;
            if divisionTime*bioTree{1}.scaleInfo.timeInterval>=900
                allResult.divisionTime=[allResult.divisionTime;1*bioTree{1}.scaleInfo.timeInterval/60,divisionTime*bioTree{1}.scaleInfo.timeInterval/60];
            end
        end
    end
end
for iframe=1:size(bioTree,2)
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            is2Node=bioTree{iframe}.node{iNode}.Out{iOut}.is2Node;
            if is2Node==1
                isDivisionNode=judgeDivisionNode(bioTree{iframe}.node{iNode});
                if isDivisionNode==1
                    nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                    isDivisionNode=judgeDivisionNode(bioTree{nodeInfo(1)}.node{nodeInfo(2)});
                    if isDivisionNode==1
                        divisionTime=nodeInfo(1)-iframe;
                        if divisionTime*bioTree{1}.scaleInfo.timeInterval>=900
                            allResult.divisionTime=[allResult.divisionTime;iframe*bioTree{1}.scaleInfo.timeInterval/60,divisionTime*bioTree{1}.scaleInfo.timeInterval/60];
                        end
                    end
                end
            end
        end
    end
end
end
function allResult=getDetachingRateInfo(bioTree,allResult)
bacteriaFrameInfo=getEachBacteriaInFrame(bioTree);
cropInfo=bioTree{1}.imageProcessingInfo.cropInfo;
limitRegionIdx=getLimitRegion(cropInfo);
allResult.detachingResult=[];
% for iframe=1:size(bioTree,2)-1    
%     bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
%     bacNum=size(bacteriaInfo,1);
%     leafNum=size(bioTree{iframe}.leavies,2);
%     for iLeaf=1:leafNum
%         leafInfo=bioTree{iframe}.leavies{iLeaf}.leaviesPixelDetail;
%         if any(ismember(leafInfo,limitRegionIdx))
%             leafNum=leafNum-1;
%         end
%     end
%     if leafNum~=0
%         allResult.detachingResult=[allResult.detachingResult;bacNum,leafNum/bacNum];
%     end
% end
allResult.detachingResult=[];
n=0;
bacPieceInfo=[];
leafPieceInfo=[];
for iframe=1:size(bioTree,2)-1
    n=n+1;
    if n==101
        n=1;
        allResult.detachingResult=[allResult.detachingResult;mean(bacPieceInfo),sum(leafPieceInfo)/mean(bacPieceInfo)];
        bacPieceInfo=[];
        leafPieceInfo=[];
    end
    bacteriaInfo=bacteriaFrameInfo{iframe}.bacteriaInfo;
    bacNum=size(bacteriaInfo,1);
    leafNum=size(bioTree{iframe}.leavies,2);
    for iLeaf=1:leafNum
        leafInfo=bioTree{iframe}.leavies{iLeaf}.leaviesPixelDetail;
        if any(ismember(leafInfo,limitRegionIdx))
            leafNum=leafNum-1;
        end
    end
    bacPieceInfo=[bacPieceInfo;bacNum];
    leafPieceInfo=[leafPieceInfo;leafNum];
end
end
function limitRegionIdx=getLimitRegion(cropInfo)
cropInfo1=bwmorph(cropInfo,'remove');
cropInfo1=imdilate(cropInfo1,true(13));
limitRegion=cropInfo1 & cropInfo;
cc=bwconncomp(limitRegion);
limitRegionIdx=cc.PixelIdxList{1,1};
end
function allResultSingleMapDrawing(allResult)
% 四张图，分裂前长度的分布，分裂后长度的分布，分裂时间的分布，detaching率随细菌数量的分布
allResultNew=load('J:\2014-02-11 F1 longTime(oilwater) Jzy\allResult.mat');
beforeDivision=allResult.beforeDivision;
afterDivision=allResult.afterDivision;
divisionTime=allResult.divisionTime;
detachingResult=allResult.detachingResult;

beforeDivision1=allResultNew.allResult.beforeDivision;
afterDivision1=allResultNew.allResult.afterDivision;
divisionTime1=allResultNew.allResult.divisionTime;
detachingResult1=allResultNew.allResult.detachingResult;

% Create figure
figure1 = figure;
scrsz=get(0,'ScreenSize');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[1 1 scrsz(3) scrsz(4)]);

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.25 0.537540805223069 0.234173505275498 0.406964091403699]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
[n,xout] = hist(afterDivision1,linspace(1,4,50));
n=n./sum(n);
plot(xout,n,'Parent',axes1,'Marker','.','LineStyle','--','Color',[0 0 1]);
xlim([1,4])
[n,xout] = hist(afterDivision,linspace(1,4,50));
n=n./sum(n);
plot(xout,n,'Parent',axes1,'Marker','.','LineStyle','--','Color',[1 0 0]);

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.519270833333326 0.53754080522307 0.234173505275498 0.406964091403699]);
box(axes2,'on');
hold(axes2,'all');

% Create plot
[n,xout] = hist(beforeDivision1,linspace(2,7,50));
n=n./sum(n);
plot(xout,n,'Parent',axes2,'Marker','.','LineStyle','--','Color',[0 0 1]);
xlim([2,7])
[n,xout] = hist(beforeDivision,linspace(2,7,50));
n=n./sum(n);
plot(xout,n,'Parent',axes2,'Marker','.','LineStyle','--','Color',[1 0 0]);

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.25 0.0816104461371087 0.234173505275498 0.406964091403699]);
box(axes3,'on');
hold(axes3,'all');

% Create xlabel
xlabel('Time(min)','FontSize',15);

% Create ylabel
ylabel('Division Time(min)','FontSize',15);

% Create plot
plot(divisionTime1(:,1),divisionTime1(:,2),'Parent',axes3,'Marker','.','LineStyle','none','Color',[0,0,1]);
plot(divisionTime(:,1),divisionTime(:,2),'Parent',axes3,'Marker','.','LineStyle','none','Color',[1,0,0]);

% Create axes
axes4 = axes('Parent',figure1,...
    'Position',[0.519270833333326 0.0816104461371071 0.234173505275498 0.406964091403699]);
box(axes4,'on');
hold(axes4,'all');

% Create xlabel
xlabel('bacteriaNum','FontSize',15);

% Create ylabel
ylabel('detaching rate','FontSize',15);

% Create plot
ylim([0,1])
plot(detachingResult1(:,1),detachingResult1(:,2),'Parent',axes4,'Marker','.','LineStyle','--',...
    'Color',[0,0,1]);
plot(detachingResult(:,1),detachingResult(:,2),'Parent',axes4,'Marker','.','LineStyle','--',...
    'Color',[1,0,0]);

% Create textbox
annotation(figure1,'textbox',...
    [0.530374999999993 0.8648 0.089625 0.0512000000000004],...
    'String',{'length before division'},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.376 0.8656 0.0969166666666667 0.0512000000000004],...
    'String',{'length after division'},...
    'FontSize',20,...
    'FitBoxToText','off',...
    'LineStyle','none');

end