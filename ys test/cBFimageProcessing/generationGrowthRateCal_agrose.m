%% each generation growth rate calculate
function [growthInfo]=generationGrowthRateCal_agrose(bioTree,dirFile)

if isfolder([dirFile,'\growthResult'])
    %     delete ([dirFile,'\growthResult\','*.mat']);
    %     delete ([dirFile,'\growthResult\','*.fig']);
    %     delete ([dirFile,'\growthResult\','*.tif']);
    delete ([dirFile,'\growthResult\','*.*']);
    rmdir([dirFile,'\growthResult']);
end

dirSave=[dirFile,'\growthResult'];
mkdir(dirSave);

growthInfo={};
growthInfo.gOne={};
growthInfo.gTwo={};
growthInfo.gThree={};
if isempty(bioTree)
    disp('bioTree is empty')
    return
end
if isempty(bioTree{1}.root)
    disp('no cells in the 1st frame')
    return
end
g1=0;
figure,
for iBac= 1:numel(bioTree{1}.root)
    % 统计root-node 生长率
    if bioTree{1}.root{iBac}.is2Node==1
        
        if bioTree{1}.root{iBac}.nodeInfo(1,1)>=12&&...
                numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{end})/...
                numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{1})>1.1
            g1=g1+1;
            growthInfo.gOne{g1}.is2Node=1;
            growthInfo.gOne{g1}.traceInfo=bioTree{1}.root{iBac}.traceInfo;
            growthInfo.gOne{g1}.rootInfo=[1,iBac];
            growthInfo.gOne{g1}.nodeInfo=bioTree{1}.root{iBac}.nodeInfo;
            growthInfo.gOne{g1}.leafInfo=[];
            for iframe=1:numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList)
                growthInfo.gOne{g1}.cellArea(iframe)=numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{iframe});
                growthInfo.gOne{g1}.MajorAxisLength(iframe)=bioTree{1}.root{iBac}.traceInfo.measurment{iframe}.MajorAxisLength;
            end
            growthInfo.gOne{g1}.timer=bioTree{1, 1}.bioTreeTimer(1:iframe);%min
            growthInfo.gOne{g1}.divTime=growthInfo.gOne{g1}.timer(end)-growthInfo.gOne{g1}.timer(1);%min
            % gR=ln(L1/L0)/deltT linear fit
            [growthInfo.gOne{g1}.gR,growthInfo.gOne{g1}.fitperiodata,growthInfo.gOne{g1}.fitgof]=prepareCurveData( growthInfo.gOne{g1}.timer(1:2:end), growthInfo.gOne{g1}.cellArea(1:2:end) );
            
            plot(growthInfo.gOne{g1}.fitperiodata, growthInfo.gOne{g1}.timer(1:2:end), log(growthInfo.gOne{g1}.cellArea(1:2:end)) );
            hold on
            saveas(gcf,[dirSave,'\','root',num2str(iBac,'%03d'),'.fig']);
            saveas(gcf,[dirSave,'\','root',num2str(iBac,'%03d'),'.tif']);
            close all;
            
            % 第二代细菌生长信息
            g2_node=bioTree{1}.root{iBac}.nodeInfo();
            [daughterInfo,nodeErrIdx]=daughterCellsGrowthRateFrOneNode(bioTree,[1,iBac],g2_node,dirSave);
            growthInfo.gTwo=[ growthInfo.gTwo,daughterInfo];
            
            % 第三代细菌生长信息
            if nodeErrIdx==0&&~isempty(bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out)&&...
                    numel(bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out)>=g2_node(1,3)*2
                for j=g2_node(1,3)*2-1:g2_node(1,3)*2 %1:numel(bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out)
                    if numel(bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out{j}.traceInfo.pixelIdxList)>12&&...
                            numel(bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out{j}.traceInfo.pixelIdxList{end})/...
                            numel(bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out{j}.traceInfo.pixelIdxList{1})>1.1
                        
                        if bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out{j}.is2Node==1
                            g3_node=bioTree{g2_node(1,1)}.node{g2_node(1,2)}.Out{j}.nodeInfo;
                            [daughterInfo,~]=daughterCellsGrowthRateFrOneNode(bioTree,g2_node,g3_node,dirSave);
                            growthInfo.gThree=[ growthInfo.gThree,daughterInfo];
                            
                        end
                    end
                end
            end
            
        end
    end
    % 统计第一代root-leaf细菌的生长率
    if bioTree{1}.root{iBac}.is2Node==0
        if bioTree{1}.root{iBac}.leafInfo(1,1)>=12&&...
                numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{end})/...
                numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{1})>1.1
            g1=g1+1;
            growthInfo.gOne{g1}.is2Node=0;
            growthInfo.gOne{g1}.traceInfo=bioTree{1}.root{iBac}.traceInfo;
            growthInfo.gOne{g1}.rootInfo=[1,iBac];
            growthInfo.gOne{g1}.leafInfo=bioTree{1}.root{iBac}.leafInfo;
            growthInfo.gOne{g1}.nodeInfo=[];
            for iframe=1:numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList)
                growthInfo.gOne{g1}.cellArea(iframe)=numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{iframe});
                growthInfo.gOne{g1}.MajorAxisLength(iframe)=bioTree{1}.root{iBac}.traceInfo.measurment{iframe}.MajorAxisLength;
            end
            growthInfo.gOne{g1}.timer=bioTree{1, 1}.bioTreeTimer(1:iframe);%min
            % gR=ln(L1/L0)/deltT fit
            [growthInfo.gOne{g1}.gR,growthInfo.gOne{g1}.fitperiodata,growthInfo.gOne{g1}.fitgof]=prepareCurveData( growthInfo.gOne{g1}.timer(1:2:end), growthInfo.gOne{g1}.cellArea(1:2:end) );
            
            figure,
            plot(growthInfo.gOne{g1}.fitperiodata, growthInfo.gOne{g1}.timer(1:2:end), log(growthInfo.gOne{g1}.cellArea(1:2:end)) );
            hold on
            saveas(gcf,[dirSave,'\','root',num2str(iBac,'%03d'),'.fig']);
            saveas(gcf,[dirSave,'\','root',num2str(iBac,'%03d'),'.tif']);
            close all;
            
        end
    end
    
end

plotAllgRfit(growthInfo,dirSave);
plotGrowthRateBasedGeneraration(growthInfo,dirSave);

[growthInfo]=motherDaughterCellgrowthRateGet(growthInfo);
save([dirSave,'\growthInfo.mat'],'growthInfo')
end
%% 从node出发寻找两个子代的生长率
function [daughterInfo,nodeErrIdx]=daughterCellsGrowthRateFrOneNode(bioTree,nodeInfo_In, nodeInfo,dirSave)
%nodeInfo_In 代表当前mode是从那个node或root来的
daughterInfo={};
k=0;
nodeErrIdx=0;
try
    bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.In{nodeInfo(1,3)};
catch
    disp('error node')
    nodeErrIdx=1;
    return
end

if ~isempty(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)&&...
        numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)>=nodeInfo(1,3)*2 %确保1个In 2个Out
    
    for i=nodeInfo(1,3)*2-1:nodeInfo(1,3)*2%numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)
        if numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList)>12&&...
                numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList{end})/...
                numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList{1})>1.1
            k=k+1;
            daughterInfo{k}.nodeInfo_In=nodeInfo_In;
            daughterInfo{k}.nodeInfo=nodeInfo;
            daughterInfo{k}.traceInfo=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo  ;
            for iframe=1:numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList)
                daughterInfo{k}.cellArea(iframe)=numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList{iframe});
                daughterInfo{k}.MajorAxisLength(iframe)=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.measurment{iframe}.MajorAxisLength;
            end
            
            daughterInfo{k}.timer=bioTree{1, 1}.bioTreeTimer(nodeInfo(1,1):nodeInfo(1,1)+iframe-1);%min
            
            if bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.is2Node==1
                daughterInfo{k}.is2Node=1;
                daughterInfo{k}.nodeInfo_Out=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.nodeInfo;
                daughterInfo{k}.leafInfo_Out=[];
                daughterInfo{k}.divTime=daughterInfo{k}.timer(end)-daughterInfo{k}.timer(1);%min
            else
                daughterInfo{k}.is2Node=0;
                daughterInfo{k}.nodeInfo_Out=[];
                daughterInfo{k}.leafInfo_Out=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.leafInfo;
            end
            
            [daughterInfo{k}.gR,daughterInfo{k}.fitperiodata,daughterInfo{k}.fitgof]=prepareCurveData( daughterInfo{k}.timer(1:2:end), daughterInfo{k}.cellArea(1:2:end) );
            figure,
            plot(daughterInfo{k}.fitperiodata, daughterInfo{k}.timer(1:2:end), log(daughterInfo{k}.cellArea(1:2:end)) );
            hold on
            saveas(gcf,[dirSave,'\','node',num2str(nodeInfo(1)),num2str(nodeInfo(2)),num2str(nodeInfo(3)),'_',num2str(i,'%02d'),'.fig']);
            saveas(gcf,[dirSave,'\','node',num2str(nodeInfo(1)),num2str(nodeInfo(2)),num2str(nodeInfo(3)),'_',num2str(i,'%02d'),'.tif']);
            close all;
            
        end
        
    end
end

end
%% growth rate fit
function [gR,fitperiodata,gof]=prepareCurveData( timeSeries, dataPickI )
% row number
timeSeries=timeSeries(:);
dataPickI=dataPickI(:);
% 取对数后 smooth
dataPickI=log(dataPickI);
dataPickI=smooth(dataPickI,'lowess');
dataPickI=smooth(dataPickI,'rlowess');

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Robust = 'Bisquare';
opts.Upper = [Inf Inf];

% Fit model to data.
[fitperiodata,gof] = fit(timeSeries,dataPickI,ft,opts);
gR=fitperiodata.p1;
end

%% plot each geneation growth rate
function plotGrowthRateBasedGeneraration(growthInfo,dirSave)
g1=[];
g2=[];
g3=[];
if numel(growthInfo.gOne)>0
    for i=1: numel(growthInfo.gOne)
        g1(1,i)=growthInfo.gOne{i}.gR;
        g1(2,i)=1;
        g1=g1(:, g1(1,:)>0);
    end
end

if numel(growthInfo.gTwo)>0
    for i=1: numel(growthInfo.gTwo)
        g2(1,i)=growthInfo.gTwo{i}.gR;
        g2(2,i)=2;
        g2=g2(:, g2(1,:)>0);
    end
end

if numel(growthInfo.gThree)>0
    for i=1: numel(growthInfo.gThree)
        g3(1,i)=growthInfo.gThree{i}.gR;
        g3(2,i)=3;
        g3=g3(:, g3(1,:)>0);
    end
end

figure,
if ~isempty(g1)
    scatter(g1(2,:),g1(1,:),50,'filled','MarkerFaceAlpha',0.5)
end
if ~isempty(g2)
    hold on
    scatter(g2(2,:),g2(1,:),50,'filled','MarkerFaceAlpha',0.5)
end
if ~isempty(g3)
    hold on
    scatter(g3(2,:),g3(1,:),50,'filled','MarkerFaceAlpha',0.5)
end


% 创建 ylabel
ylabel('gR (min-1)');

% 创建 xlabel
xlabel('Generation');

% 取消以下行的注释以保留坐标区的 X 范围
xlim(gca,[0.9 3]);
% 设置其余坐标区属性
set(gca,'FontSize',12,'XTick',[1 2 3]);

saveas(gcf,[dirSave,'\','gR_generation','.fig']);
saveas(gcf,[dirSave,'\','gR_generation','.tif']);
close all;
end
%% plot all growth fit curves
function plotAllgRfit(growthInfo,dirSave)
figure,
if numel(growthInfo.gOne)>0
    for i=1: numel(growthInfo.gOne)
        plot(growthInfo.gOne{i}.fitperiodata,growthInfo.gOne{i}.timer(1:2:end),log(growthInfo.gOne{i}.cellArea(1:2:end)));
        legend('off')
        hold on
    end
end

if numel(growthInfo.gTwo)>0
    for i=1: numel(growthInfo.gTwo)
        plot(growthInfo.gTwo{i}.fitperiodata,growthInfo.gTwo{i}.timer(1:2:end),log(growthInfo.gTwo{i}.cellArea(1:2:end)));
        legend('off')
        hold on
    end
end

if numel(growthInfo.gThree)>0
    for i=1: numel(growthInfo.gThree)
        plot(growthInfo.gThree{i}.fitperiodata,growthInfo.gThree{i}.timer(1:2:end),log(growthInfo.gThree{i}.cellArea(1:2:end)));
        legend('off')
        hold on
    end
end
hold off
saveas(gcf,[dirSave,'\','allgRfit','.fig']);
saveas(gcf,[dirSave,'\','allgRfit','.tif']);
close all;
end
%% plot the nth generation gR vs (n+1)th generation growthrate
function [growthInfo]=motherDaughterCellgrowthRateGet(growthInfo)
%gR(n+1) vs gR(n)
if numel(growthInfo.gTwo)>0
    gR_21=zeros(numel(growthInfo.gTwo),2);
    for i=1: numel(growthInfo.gTwo)
        nodeInfo=growthInfo.gTwo{i}.nodeInfo;
        gR_21(i,1)=growthInfo.gTwo{i}.gR;
        for j=1: numel(growthInfo.gOne)
            if isequal(nodeInfo,growthInfo.gOne{j}.nodeInfo)
                gR_21(i,2)=growthInfo.gOne{j}.gR;
            end
        end
    end
    growthInfo.gR_21=gR_21;
else
    growthInfo.gR_21=[];
end
if numel(growthInfo.gThree)>0
    gR_32=zeros(numel(growthInfo.gThree),2);
    for i=1: numel(growthInfo.gThree)
        nodeInfo=growthInfo.gThree{i}.nodeInfo_In;
        gR_32(i,1)=growthInfo.gThree{i}.gR;
        for j=1: numel(growthInfo.gTwo)
            if isequal(nodeInfo,growthInfo.gTwo{j}.nodeInfo)
                gR_32(i,2)=growthInfo.gTwo{j}.gR;
            end
        end
    end
    growthInfo.gR_32=gR_32;
else
    growthInfo.gR_32=[];
end
end

