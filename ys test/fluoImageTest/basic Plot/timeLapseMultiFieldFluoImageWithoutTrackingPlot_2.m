function [allDataCollect] = timeLapseMultiFieldFluoImageWithoutTrackingPlot_2(dirFile,fixedChannel)
% 用于对经过[bioInfo] = cellFluoIntensityCorrection(bioInfo)光强进行校正后的作图
% bioInfo 中有bioInfo结构体中有类似intsfGFP_c的field
% 未经过串光校正的荧光强度作图
% allCorrectedFluoName = {'intsfGFP_c','intVenus_c',...
%     'intCyOFP_c','intmScarletI_c','intTDsmURFP_c'};% c represents for corrected
% shuai Yang 2020.06.09
% fluoChannels_c = {'sfGFP','Venus','CyOFP','mScarletI','TDsmURFP'};%无PVD，PVD暂不校正
% update 将每个通道的光强进行筛选要求>0 去掉,2021/10/13
% 只要求fixedChannel的光强>intThreshold,2021/10/13
% 取消了位置的筛选和光强的筛选 2021/10/28
 
% allFluoName_c = {'intsfGFP','intVenus','intCyOFP',...
%     'intmScarletI','intTDsmURFP'};%校正的荧光通道

% scale 获得
try
    load([dirFile,filesep,'mip.mat']);
    scale = mip.Calib.scale ;
catch
    scale = 0.065;%100x pixel to um %default
end

% intThreshold = 5;%已经减过背景 default 5 最小为0 intensity Threshold
% if strcmp(fixedChannel,'BF1')
%     %     xPosLimit = [244,1860];
%     %BF的mask，荧光图像的视野范围小
%     xPosLimit = [0,2048];%2020.12.30 视野照明范围变大了 不用再截了
% else
%     xPosLimit = [0,2048];
% end

dirSave = strcat(dirFile,'\result_basic2');
mkdir(dirSave);

dirTimeBioInfo = strcat(dirFile,'\timeBioInfo');
timeBioInfoList = dir([dirTimeBioInfo,filesep,'bioInfo*']);
allDataCollect = cell(1, length(timeBioInfoList));
k = 1;
for iTime = 1:length(timeBioInfoList)
    if ~strcmp(timeBioInfoList(iTime).name(1:7),'bioInfo')
        continue
    end
    load([dirTimeBioInfo,'\',timeBioInfoList(iTime).name]);
    dataCollect = struct('intsfGFP',[ ],'intmScarletI',[ ],'intVenus',[ ],...
        'intPVD',[ ],'intCyOFP',[ ],'intTDsmURFP',[ ],'MajorAxisLength',[ ],...
        'Centroid',[ ]);

    for i = 1:numel(bioInfo)
        dataCollect.intsfGFP(end+1:end+size(bioInfo(i).intsfGFP,1),:) = bioInfo(i).intsfGFP_c ;
        dataCollect.intmScarletI(end+1:end+size(bioInfo(i).intmScarletI,1),:) = bioInfo(i).intmScarletI_c;
        dataCollect.intVenus(end+1:end+size(bioInfo(i).intVenus,1),:) = bioInfo(i).intVenus_c;
        dataCollect.intPVD(end+1:end+size(bioInfo(i).intPVD,1),:) = bioInfo(i).intPVD - bioInfo(i).BG_PVD; %PVD无校正的数据 直接减背景
        dataCollect.intCyOFP(end+1:end+size(bioInfo(i).intCyOFP,1),:) = bioInfo(i).intCyOFP_c;
        dataCollect.intTDsmURFP(end+1:end+size(bioInfo(i).intTDsmURFP,1),:) = bioInfo(i).intTDsmURFP_c;
        dataCollect.MajorAxisLength(end+1:end+size(bioInfo(i).MajorAxisLength,1),:) = bioInfo(i).MajorAxisLength;
        dataCollect.Centroid(end+1:end+size(bioInfo(i).Centroid,1),:) = bioInfo(i).Centroid;
    end

    % 通过Centroid判断是否有数据点 如果为空说明此时间点没有细菌的数据
    if isempty(dataCollect.Centroid)
        continue
    end

    % 删除为空的结构体dataCollect的field
    fields = fieldnames(dataCollect);%结构体的fields of struct
    for iFds = 1 :numel(fields)%iFds is short for Fields of struct
        if isempty(dataCollect.(fields{iFds}))
            dataCollect = rmfield(dataCollect,fields{iFds});
        end
    end
    %
    %     fields = fieldnames(dataCollect);
    %
    %     posTF = (dataCollect.Centroid(:,1)>= xPosLimit(1)&dataCollect.Centroid(:,1)<= xPosLimit(2));
    %     if strcmp(fixedChannel,'BF1')
    %         % intTF = dataCollect.(fields{1})> intThreshold;%光强的筛选
    %         intTF = true(size(posTF));
    %     else
    %         intTF = dataCollect.(strcat('int',fixedChannel))> intThreshold;
    %     end
    %     TF = posTF & intTF;
    %
    %     for iFds = 1 : numel(fields)
    %         if strlength(fields{iFds})<3% case for string less than 3 letters
    %             continue
    %         end
    %         if strcmp(fields{iFds}(1:3),'int')
    %             if ~ismember(fields{iFds},allFluoName_c)
    %                 %没有校正的荧光通道光强减去100，其他在cellFluoIntensityCorrection已经减过
    %                 dataCollect.(fields{iFds}) =  dataCollect.(fields{iFds})-100;
    %             end
    %             % eachChIntTF = dataCollect.(fields{iFds})> 0;% 对所有荧光通道进行筛选 保证为正值
    %             % TF = TF & eachChIntTF;
    %         end
    %     end
    %
    %     for iFds = 1 : numel(fields)
    %         dataCollect.(fields{iFds}) = dataCollect.(fields{iFds})(TF,:);
    %     end
    %
    dataCollect = rmfield(dataCollect,'Centroid');
    allDataCollect{iTime} = dataCollect;

    if k == 1
        save(strcat(dirSave,'\dataCollect',timeBioInfoList(iTime).name(8:end-4),'.mat'),'dataCollect','-v7.3');
        %plot function
        basicInfoPlot(dataCollect,scale,dirSave)
    end
    k = k+1;
end
save(strcat(dirSave,'\allDataCollect.mat'),'allDataCollect','-v7.3');
intensityChangePlot(allDataCollect,dirSave,dirFile);
% intensityChangePlot_scatterEachPoint(allDataCollect,dirSave,dirFile);
end
%% plot function ：荧光强度和菌长基本信息
function basicInfoPlot(dataCollect,scale,dirSave)
fields = fieldnames(dataCollect);
cellNum = numel(dataCollect.(fields{1}));%细菌数目
channelNum = 0;
for iFds = 1:numel(fields)
    if strlength(fields{iFds})<3
        continue
    end
    if strcmp(fields{iFds}(1:3),'int')
        channelNum = channelNum + 1;
    end
end
fname = 'basicHist';
f1 = figure('Name',fname);
subplot1 = subplot(2,1,1,'Parent',f1);
for iFds = 1:channelNum %除去MajorAxisLength fields,-1
    meanValue = mean(dataCollect.(fields{iFds}));
    stdValue = std(dataCollect.(fields{iFds}));
    dataName = strcat(fields{iFds}(4:end),32, num2str(meanValue,'%.0f'),'±',num2str(stdValue,'%.0f'));
    histogram(dataCollect.(fields{iFds}),'DisplayName',dataName,'EdgeColor','none')
    hold(subplot1,'on')
end
set(subplot1,'FontSize',12);
box(subplot1,'off');
legend(subplot1,'show');
hold(subplot1,'off')
xlabel({'intensity (a.u.)'});
ylabel({'cell Num (a.u.)'});

subplot2 = subplot(2,1,2,'Parent',f1);
meanValue = mean(dataCollect.MajorAxisLength)*scale;
stdValue = std(dataCollect.MajorAxisLength)*scale;
dataName = strcat('cellLength',32, num2str(meanValue,'%.1f'),'±',num2str(stdValue,'%.1f'),'um');
histogram(dataCollect.MajorAxisLength*scale,'DisplayName',dataName,'FaceColor','k','EdgeColor','none')

set(subplot2,'FontSize',12);
box(subplot2,'off');
legend(subplot2,'show');
xlabel({'cell Len (um)'});
ylabel({'cell Num (a.u.)'});
saveas(f1,[dirSave,'\',fname,'.fig']);
saveas(f1,[dirSave,'\',fname,'.tif'])

close all

% 每一个channel光强与其他channel光强对应信息 需要至少两个channel荧光信息
if channelNum == 1
    return
end
fname = 'intensityScatter';
if channelNum>3
    f1 = figure('WindowState','maximized','Name',fname);
else
    f1 = figure('Name',fname);
end

if channelNum >= 2 && channelNum <=3
    plotNum = channelNum*(channelNum-1)/2;%每两两作图
    n = 1;%用于计数
    for iCh = 1:channelNum
        for jCh = iCh+1 :channelNum
            subplt = subplot(1,plotNum,n);
            R = corr2(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}));
            dataName = strcat('R = ',32,num2str(R,'%.2f') );
            if cellNum >10000
                histogram2(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),'BinMethod','scott','DisplayStyle','tile');
            else
                scatter(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),20,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.1);
            end
            %             scatterhist(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),'Kernel','on')
            set(subplt,'FontSize',12);
            box(subplt,'off');
            xlabel(strcat(fields{iCh}(4:end),32,'(a.u.)'));
            ylabel(strcat(fields{jCh}(4:end),32,'(a.u.)'));
            title(dataName)
            hold(subplt,'off');
            n = n + 1;
        end
    end
end
if channelNum >= 4
    plotNum = channelNum*(channelNum-1)/2;
    n = 1;%用于计数
    for iCh = 1:channelNum
        for jCh = iCh+1:channelNum

            subplt = subplot(2,ceil(plotNum/2),n,'Parent',f1);
            R = corr2(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}));
            dataName = strcat('R = ',32,num2str(R,'%.2f') );
            if cellNum >10000
                histogram2(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),'BinMethod','scott','DisplayStyle','tile');
            else
                scatter(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),20,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.1);
            end
            %             scatterhist(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),'Kernel','on')
            set(subplt,'FontSize',12);
            box(subplt,'off');
            xlabel(strcat(fields{iCh}(4:end),32,'(a.u.)'));
            ylabel(strcat(fields{jCh}(4:end),32,'(a.u.)'));
            title(dataName)
            hold(subplt,'off');
            n = n + 1;
        end
    end
end
saveas(f1,[dirSave,'\',fname,'.fig']);
saveas(f1,[dirSave,'\',fname,'.tif'])
close all
end
%%
function intensityChangePlot(allDataCollect,dirSave,dirFile)

%case for no cells in this timepoint 有可能某个时间点数据为空
for iTime = 1:numel(allDataCollect)
    if ~isempty(allDataCollect{iTime})
        fields = fieldnames(allDataCollect{iTime});
        break
    end
end

% case for empty sample
if ~exist('fields','var')
    disp('Empty sample. No Data. ')
    warning('Empty sample. No Data.');
    return
end

channelNum = 0;
for iFds = 1:numel(fields)
    if strlength(fields{iFds})<3
        continue
    end
    if strcmp(fields{iFds}(1:3),'int')
        channelNum = channelNum + 1;
    end
end
if channelNum ==0
    return
end
% 将横坐标转化为时间 非frame ys 2020.08.14
fieldList = dir([dirFile,filesep,'field*']);
load([dirFile,'\',fieldList(1).name,'\Tracking\frameInfo.mat']);
frameInfo = frameInfo(1:numel(allDataCollect),:);
time_lag = abs(etime (frameInfo(1,1:6),frameInfo(:,1:6))/60);

fname = 'intenisytChange';
f1 = figure('Name',fname);
for iFds = 1:channelNum
    fluoInt = zeros(numel(allDataCollect),2);
    for iTime = 1:numel(allDataCollect)
        try
            fluoInt(iTime,1) = mean(allDataCollect{iTime}.(fields{iFds}));
            fluoInt(iTime,2) = std(allDataCollect{iTime}.(fields{iFds}));
        catch
            msg = ['No cell was dectected in frame ',num2str(iTime),';',...
                ' Specify the ',fields{iFds}(4:end),' values as NaN.'];
            warning(msg);
            disp(msg);
            fluoInt(iTime,1) = NaN;
            fluoInt(iTime,2) = NaN;
        end
    end
    subplot1 = subplot(channelNum,1,iFds,'Parent',f1);
    %     errorbar(1:iTime,fluoInt(:,1),fluoInt(:,2),'o')
    %     line(1:iTime,fluoInt(:,1),'Color',[0,0.45,0.74],'LineWidth',2,'Marker','o' );
    %     line(1:iTime,fluoInt(:,1),'Color',[0,0.45,0.74],'LineWidth',1 );
    line(time_lag,fluoInt(:,1),'MarkerSize',4,'Marker','o','Color',[0,0.45,0.74],'LineWidth',1 )
    set(subplot1,'FontSize',12);
    ylabel(strcat(fields{iFds}(4:end),32,'(a.u.)'));
end
xlabel({'time (min)'});
% xlabel({'frame'});

saveas(f1,[dirSave,'\',fname,'.fig']);
saveas(f1,[dirSave,'\',fname,'.tif'])
close all
end
%%
function intensityChangePlot_scatterEachPoint(allDataCollect,dirSave,dirFile)
% 对每个时间点的细菌荧光做scatter图

for iTime = 1:numel(allDataCollect)
    if ~isempty(allDataCollect{iTime}) %case for no cells in this timepoint
        fields = fieldnames(allDataCollect{iTime});
        break
    end
end
% case for empty sample
if ~exist('fields','var')
    disp('Empty sample. No Data. ')
    return
end

channelNum = 0;
for iFds = 1:numel(fields)
    if strlength(fields{iFds})<3
        continue
    end
    if strcmp(fields{iFds}(1:3),'int')
        channelNum = channelNum + 1;
    end
end
if channelNum ==0
    return
end
% 将横坐标转化为时间 非frame ys 2020.08.14
fieldList = dir([dirFile,filesep,'field*']);
load([dirFile,'\',fieldList(1).name,'\Tracking\frameInfo.mat']);
frameInfo = frameInfo(1:numel(allDataCollect),:);
time_lag = abs(etime (frameInfo(1,1:6),frameInfo(:,1:6))/60);

% map = colormap(jet(256));% 自定义颜色数目
map = linspecer(256,'seq');

fname = 'intenisytChange_scatterEachCell';
f1 = figure('Name',fname);
for iFds = 1:channelNum
    fluoInt = zeros(numel(allDataCollect),2);
    eachCellInt = [];
    for iTime = 1:numel(allDataCollect)
        % case for empty
        if isempty(allDataCollect{iTime})
            msg = ['No cell was dectected in frame ',num2str(iTime),';',...
                ' Specify the ',fields{iFds}(4:end),' values as NaN.'];
            warning(msg);
            disp(msg);
            fluoInt(iTime,1) = NaN;
            fluoInt(iTime,2) = NaN;
            eachCellInt = NaN(1,3);
            continue
        end
        tempInt = zeros(size(allDataCollect{iTime}.(fields{iFds}),1),3);
        tempInt(:,1) = allDataCollect{iTime}.(fields{iFds});
        tempInt(:,2) = time_lag(iTime);

        %     omitted repeated mnumbers Shuai Yang
        %  这一步可以先用ksdensity计算概率 再取整不要相同的数据 或者相近的数据
        %  以进一步减少作图的点
        %  Shuai Yang 2022/5/18

        intValue = round(tempInt(:,1));
        [~, ia, ~] = unique(intValue,'stable');
        tempInt = tempInt(ia,:);

        colorIdx = getDensityColor(tempInt(:,1));
        tempInt(:,3) = colorIdx;
        eachCellInt = [eachCellInt;tempInt];

        fluoInt(iTime,1) = mean(allDataCollect{iTime}.(fields{iFds}));
        fluoInt(iTime,2) = std(allDataCollect{iTime}.(fields{iFds}));

    end
    subplot1 = subplot(channelNum,1,iFds,'Parent',f1);
    %     errorbar(1:iTime,fluoInt(:,1),fluoInt(:,2),'o')
    %     line(1:iTime,fluoInt(:,1),'Color',[0,0.45,0.74],'LineWidth',2,'Marker','o' );
    %     line(1:iTime,fluoInt(:,1),'Color',[0,0.45,0.74],'LineWidth',1 );
    line(time_lag,fluoInt(:,1),'MarkerSize',4,'Marker','o','Color',[0,0.45,0.74],'LineWidth',1 );
    hold on

    TF = ~isnan(eachCellInt(:,1));
    eachCellInt = eachCellInt(TF,:);

    %     line(eachCellInt(:,2),eachCellInt(:,1),'MarkerSize',4,'Marker','o','Color',[0,0.45,0.74],'LineWidth',0.5,'LineStyle','none' );
    %     scatter(eachCellInt(:,2),eachCellInt(:,1),15,[0.85,0.33,0.10],'filled','MarkerFaceAlpha',0.05);
    scatter(eachCellInt(:,2),eachCellInt(:,1),15,map(eachCellInt(:,3),:),'filled','MarkerFaceAlpha',0.5);

    %     for k = 1 :size(eachCellInt,1)
    %         line('Xdata',eachCellInt(k,2),'Ydata',eachCellInt(k,1), ...
    %             'LineStyle','none','Color',map(eachCellInt(k,3),:), ...
    %             'Marker','.','MarkerSize',12);
    %     end

    set(subplot1,'FontSize',12);
    ylabel(strcat(fields{iFds}(4:end),32,'(a.u.)'));
end
xlabel({'time (min)'});
% xlabel({'frame'});

saveas(f1,[dirSave,'\',fname,'.fig']);
saveas(f1,[dirSave,'\',fname,'.tif'])
close all
end

function colorIdx = getDensityColor(data)
% based on ksdensity
%列数据
data = data(:);
colorNum = 256;
colorIdx = zeros(size(data,1),1);
for i = 1:size(data,1)
    [f,~] = ksdensity(data(:,1),data(i,1));

    colorIdx(i) = f;
end

colorIdx = rescale(colorIdx,1,colorNum);
colorIdx = round(colorIdx);
end


