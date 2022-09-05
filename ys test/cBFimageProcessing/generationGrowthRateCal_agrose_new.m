%% each generation growth rate calculate
function [growthData,growthInfo] = generationGrowthRateCal_agrose_new(bioTree,dirFile)
% comments :Nice code Shuai Yang 2021.08.11
% if isfolder([dirFile,'\growthResult'])
%     %     delete ([dirFile,'\growthResult\','*.mat']);
%     %     delete ([dirFile,'\growthResult\','*.fig']);
%     %     delete ([dirFile,'\growthResult\','*.tif']);
%     delete ([dirFile,'\growthResult\','*.*']);
%     rmdir([dirFile,'\growthResult']);
% end
fieldTag = str2double(dirFile(end-3:end)); 
dirSave = [dirFile,'\growthResult2'];
mkdir(dirSave);
% initialization
growthInfo = {};
growthData = {};

if isempty(bioTree)
    disp('bioTree is empty')
    return
end
if isempty(bioTree{1}.root)
    disp('no cells in the 1st frame')
    return
end

n1 = 0;%第一代细菌数目
traceNum = 3;% 计算gR需要追踪的frame的阈值
for iBac = 1:numel(bioTree{1}.root)
    % 统计root-node 生长率
    g=1;%细菌代数
    if bioTree{1}.root{iBac}.is2Node == 1
        
        if bioTree{1}.root{iBac}.nodeInfo(1,1) >= traceNum %default 12
            cellArea = zeros(1,numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList));
            MajorAxisLength = zeros(1,numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList));
            Centroid = zeros(numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList),2);
            for iFrame = 1:numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList)
                cellArea(iFrame) = numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{iFrame});
                MajorAxisLength(iFrame) = bioTree{1}.root{iBac}.traceInfo.measurment{iFrame}.MajorAxisLength;
                Centroid(iFrame,:) = bioTree{1}.root{iBac}.traceInfo.measurment{iFrame}.Centroid;
            end
            timer = bioTree{1, 1}.bioTreeTimer(1:iFrame);%min
            if isequal(timer(1),timer(2)) %用面积计算生长率，也可以用菌长MajorAxisLength
                [gR,fitperiodata,fitgof] = prepareCurveData( timer(1:2:end), cellArea(1:2:end) );
                %for case: bioTree 生成,1,1,2,2....
            else
                [gR,fitperiodata,fitgof] = prepareCurveData(timer, cellArea);
            end
            
            %             if sum(diff(cellArea))/cellArea(1)>0.1&&gR>0.001&&gR<0.1&&fitgof.rsquare>=0.8    %至少增长原来的10%
            if gR > 0.001 && gR < 0.1 && fitgof.rsquare >= 0.8
                n1 = n1+1;
                growthInfo{g}(n1).fieldTag = fieldTag;
                growthInfo{g}(n1).genrTag = g;
                growthInfo{g}(n1).is2Node = 1;
                growthInfo{g}(n1).rootInfo = [1,iBac];
                growthInfo{g}(n1).nodeInfo = bioTree{1}.root{iBac}.nodeInfo;
                growthInfo{g}(n1).leafInfo = [];
                growthInfo{g}(n1).traceInfo = bioTree{1}.root{iBac}.traceInfo;
                growthInfo{g}(n1).Centroid = Centroid;
                growthInfo{g}(n1).absCentroid = [];
                growthInfo{g}(n1).cellArea = cellArea;
                growthInfo{g}(n1).MajorAxisLength = MajorAxisLength;
                
                growthInfo{g}(n1).timer = timer;%min
                growthInfo{g}(n1).divTime = growthInfo{g}(n1).timer(end)-growthInfo{g}(n1).timer(1);%min
                %             gR=ln(L1/L0)/deltT linear fit
                growthInfo{g}(n1).gR = gR;
                growthInfo{g}(n1).fitperiodata = fitperiodata;
                growthInfo{g}(n1).fitgof = fitgof;
                
                % 第一代之后细菌的生长信息
                g2_node = bioTree{1}.root{iBac}.nodeInfo();
                [growthInfo] = daughterCellsGrowthRateFrOneNode(bioTree,growthInfo,[1,iBac],g2_node,g);
                
            end
        end
    end
    %     统计第一代root-leaf细菌的生长率 分开统计 如果数据不需要 可以直接删除
    if bioTree{1}.root{iBac}.is2Node == 0
        if bioTree{1}.root{iBac}.leafInfo(1,1) >= traceNum
            
            cellArea = zeros(1,numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList));
            MajorAxisLength = zeros(1,numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList));
            Centroid = zeros(numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList),2);
            for iFrame = 1:numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList)
                cellArea(iFrame) = numel(bioTree{1}.root{iBac}.traceInfo.pixelIdxList{iFrame});
                MajorAxisLength(iFrame) = bioTree{1}.root{iBac}.traceInfo.measurment{iFrame}.MajorAxisLength;
                Centroid(iFrame,:) = bioTree{1}.root{iBac}.traceInfo.measurment{iFrame}.Centroid;
            end
            timer = bioTree{1, 1}.bioTreeTimer(1:iFrame);%min
            
            if isequal(timer(1),timer(2)) %用面积计算生长率，也可以用菌长MajorAxisLength
                [gR,fitperiodata,fitgof] = prepareCurveData( timer(1:2:end), cellArea(1:2:end) );
            else
                [gR,fitperiodata,fitgof] = prepareCurveData(timer, cellArea);
            end
            
            %             if sum(diff(cellArea))/cellArea(1)>0.1
            if gR > 0.001 && gR < 0.1 && fitgof.rsquare >= 0.8
                n1=n1+1;
                growthInfo{g}(n1).fieldTag = fieldTag;
                growthInfo{g}(n1).genrTag = g;
                growthInfo{g}(n1).is2Node = 0;               
                growthInfo{g}(n1).rootInfo = [1,iBac];
                growthInfo{g}(n1).leafInfo = bioTree{1}.root{iBac}.leafInfo;
                growthInfo{g}(n1).nodeInfo = [];
                growthInfo{g}(n1).traceInfo = bioTree{1}.root{iBac}.traceInfo;
                growthInfo{g}(n1).Centroid = Centroid;
                growthInfo{g}(n1).absCentroid = [];
                growthInfo{g}(n1).cellArea = cellArea;
                growthInfo{g}(n1).MajorAxisLength = MajorAxisLength;
                
                growthInfo{g}(n1).timer = timer;%min
                growthInfo{g}(n1).divTime = [];%leaf不统计divisionTime          
                %             gR=ln(L1/L0)/deltT linear fit
                growthInfo{g}(n1).gR = gR;
                growthInfo{g}(n1).fitperiodata = fitperiodata;
                growthInfo{g}(n1).fitgof = fitgof;
            end
        end
    end
    
end

plotAllgRfit(growthInfo,dirSave);
% [gRGet] = plotGrowthRateBasedGeneraration(growthInfo,dirSave);
% [gR_kin] = motherDaughterCellgrowthRateGet(growthInfo);

growthData.growthInfo = growthInfo;
% growthData.gRGet = gRGet;
% growthData.gR_kin = gR_kin;
save([dirSave,'\growthData.mat'],'growthData')
end
%% 从node出发寻找两个子代的生长率
function [growthInfo]=daughterCellsGrowthRateFrOneNode(bioTree,growthInfo,nodeInfo_In,nodeInfo,g)
g=g+1;
%获取第g代现有的细菌数目k
try
    growthInfo{1,g};
    k= numel(growthInfo{1,g});
catch
    k=0;
end
%nodeInfo_In 代表当前node是从那个node或root来的
try
    bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.In{nodeInfo(1,3)};
catch
    disp('error node')
    return
end

if ~isempty(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)&&...
        numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)>=nodeInfo(1,3)*2 %确保有Out以及1个In 2个Out
    
    for i=nodeInfo(1,3)*2-1:nodeInfo(1,3)*2%numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)
        if numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList)>12
            cellArea=zeros(1,numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList));
            MajorAxisLength=zeros(1,numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList));
            Centroid=zeros(numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList),2);
            for iframe=1:numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList)
                cellArea(iframe)=numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.pixelIdxList{iframe});
                MajorAxisLength(iframe)=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.measurment{iframe}.MajorAxisLength;
                Centroid(iframe,:)=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo.measurment{iframe}.Centroid;
            end
            timer=bioTree{1, 1}.bioTreeTimer(nodeInfo(1,1):nodeInfo(1,1)+iframe-1);%min
            [gR,fitperiodata,fitgof]=prepareCurveData( timer(1:2:end), cellArea(1:2:end) );
            if gR>0.001&&gR<0.1&&fitgof.rsquare>=0.8
                k=k+1;
                growthInfo{g}(k).fieldTag=growthInfo{1}(1).fieldTag;
                growthInfo{g}(k).genrTag=g;
                growthInfo{g}(k).is2Node=[];
                growthInfo{g}(k).nodeInfo_In=nodeInfo_In;
                growthInfo{g}(k).nodeInfo=nodeInfo;
                growthInfo{g}(k).nodeInfo_Out=[];
                growthInfo{g}(k).leafInfo_Out=[];
                growthInfo{g}(k).traceInfo=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.traceInfo  ;
                growthInfo{g}(k).Centroid=Centroid;
                growthInfo{g}(k).absCentroid=[];
                growthInfo{g}(k).cellArea=cellArea;
                growthInfo{g}(k).MajorAxisLength=MajorAxisLength;
                growthInfo{g}(k).timer=timer;%min
                growthInfo{g}(k).divTime=[];
                growthInfo{g}(k).gR=gR;
                growthInfo{g}(k).fitperiodata=fitperiodata;
                growthInfo{g}(k).fitgof=fitgof;
                if bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.is2Node==1
                    growthInfo{g}(k).is2Node=1;
                    growthInfo{g}(k).nodeInfo_Out=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.nodeInfo;                   
                    growthInfo{g}(k).divTime=growthInfo{g}(k).timer(end)-growthInfo{g}(k).timer(1);%min
                    next_node=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.nodeInfo;
                    [growthInfo]=daughterCellsGrowthRateFrOneNode(bioTree,growthInfo,nodeInfo,next_node,g);
                else
                    growthInfo{g}(k).is2Node=0;
                    growthInfo{g}(k).leafInfo_Out=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{i}.leafInfo;
                end                
            end
        end
        
    end
end

end
%% growth rate fit
function [gR,fitperiodata,gof] = prepareCurveData( timeSeries, dataPickI )
% row number
timeSeries = timeSeries(:);
dataPickI = dataPickI(:);
% 取对数后 smooth
dataPickI = log(dataPickI);
dataPickI = smooth(dataPickI,'lowess');
dataPickI = smooth(dataPickI,'rlowess');

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Robust = 'Bisquare';
opts.Upper = [Inf Inf];
warning('off');
% Fit model to data.
[fitperiodata,gof] = fit(timeSeries,dataPickI,ft,opts);
gR = fitperiodata.p1;
end

%% plot each geneation growth rate
function [gRGet]=plotGrowthRateBasedGeneraration(growthInfo,dirSave)

if numel(growthInfo)==0
    gRGet={};
    return
end

gRGet=cell(1,numel(growthInfo));

for i=1:numel(growthInfo)
    for j=1:numel(growthInfo{i})
        gRGet{i}(1,j)=i;
        gRGet{i}(2,j)=growthInfo{i}(j).gR;
        gRGet{i}(3,j)=growthInfo{i}(j).fitgof.rsquare;
        
    end
    scatter(gRGet{i}(1,:),gRGet{i}(2,:),50,'filled','MarkerFaceAlpha',0.5)
    hold on
end
% 创建 ylabel
ylabel('gR (min-1)');

% 创建 xlabel
xlabel('Generation');

% 取消以下行的注释以保留坐标区的 X 范围
xlim(gca,[0.9 i]);
% 设置其余坐标区属性
ylim(gca,[0.001 0.1]);
set(gca,'FontSize',12,'XTick',(1: i),'YScale','log');

saveas(gcf,[dirSave,'\','gR_generation','.fig']);
saveas(gcf,[dirSave,'\','gR_generation','.tif']);
hold off
close all;
end
%% plot all growth fit curves
function plotAllgRfit(growthInfo,dirSave)

if numel(growthInfo)==0
    return
end
% % 单个作图
for i=1:numel(growthInfo)
    for j=1:numel(growthInfo{i})
        figure,
        set(gcf,'visible','off')
        plot(growthInfo{i}(j).fitperiodata,growthInfo{i}(j).timer(1:2:end),log(growthInfo{i}(j).cellArea(1:2:end)))
        title(['gR=',num2str(growthInfo{i}(j).gR,'%.04f'),32,'rsquare=',num2str(growthInfo{i}(j).fitgof.rsquare,'%.03f' )]);
%         saveas(gcf,[dirSave,'\','g',num2str(i),'_',num2str(j,'%03d'),'.fig']);
        saveas(gcf,[dirSave,'\','g',num2str(i),'_',num2str(j,'%03d'),'.tif']);
        close all;
    end
    
end
% %全部做到一张图
figure,
set(gcf,'visible','off')
for i=1:numel(growthInfo)
    for j=1:numel(growthInfo{i})
        plot(growthInfo{i}(j).fitperiodata,growthInfo{i}(j).timer(1:2:end),log(growthInfo{i}(j).cellArea(1:2:end)))
        legend('off')
        hold on
    end
end
hold off
% saveas(gcf,[dirSave,'\','allgRfit','.fig']);
saveas(gcf,[dirSave,'\','allgRfit','.tif']);
close all;

end
%% plot the nth generation gR vs (n+1)th generation growthrate
function [gR_kin]=motherDaughterCellgrowthRateGet(growthInfo)
%gR(n+1) vs gR(n) son and mother cell growth rate
%gR_kin{1}means [g2,g1];gR_kin{2}means [g3,g2];

if numel(growthInfo)<2
    gR_kin={};
    return
end
gR_kin=cell(1,numel(growthInfo)-1);
for i=2:numel(growthInfo)
    gR_kin{i-1}=zeros(numel(growthInfo{i}),2);
    for j=1:numel(growthInfo{i})
        nodeInfo=growthInfo{i}(j).nodeInfo;
        gR_kin{i-1}(j,1)=growthInfo{i}(j).gR;
        if i==2
            for k=1:numel(growthInfo{i-1})
                if isequal(nodeInfo,growthInfo{i-1}(k).nodeInfo)
                    gR_kin{i-1}(j,2)=growthInfo{i-1}(k).gR;
                end
            end
        else
            for k=1:numel(growthInfo{i-1})
                if isequal(nodeInfo,growthInfo{i-1}(k).nodeInfo_Out)
                    gR_kin{i-1}(j,2)=growthInfo{i-1}(k).gR;
                end
            end
        end
        
    end
    
end

end

