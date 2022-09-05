function [dataCollect] = timeLapseMultiFieldFluoImagePlot(dirFile,fixedChannel)
% timeLapseMontageFluoImagePlot
% 未经过串光校正的荧光强度作图
% shuai Yang 2020.05.26
% update 将每个通道的光强进行筛选要求>0 去掉, 2021/10/13
% 只要求fixedChannel的光强>intThreshold,2021/10/13
% 取消了位置的筛选和光强的筛选 2021/10/28

% scale 获得
try
    load([dirFile,filesep,'mip.mat']);
    scale = mip.Calib.scale ; 
catch
    scale = 0.065;%100x pixel to um %default
end

% intThreshold = 105;%没有减背景 最小为100  intensity Threshold
% if strcmp(fixedChannel,'BF1')
%     xPosLimit = [0,2048];
%     %BF的mask，荧光图像的视野范围小
% else
%     xPosLimit = [0,2048];
% end

load([dirFile,'\timeBioInfo\bioInfo_t0000.mat']);
dirSave =strcat(dirFile,'\result_basic');
mkdir(dirSave);
dataCollect = struct('intsfGFP',[ ],'intmScarletI',[ ],'intVenus',[ ],...
    'intPVD',[ ],'intCyOFP',[ ],'intTDsmURFP',[ ],'MajorAxisLength',[ ],...
    'gR',[ ],'Centroid',[ ]);

for i = 1:numel(bioInfo)
    % 荧光强度减去每一张图像的荧光背景信号
    dataCollect.intsfGFP(end+1:end+size(bioInfo(i).intsfGFP,1),:) = bioInfo(i).intsfGFP - bioInfo(i).BG_sfGFP;
    dataCollect.intmScarletI(end+1:end+size(bioInfo(i).intmScarletI,1),:) = bioInfo(i).intmScarletI - bioInfo(i).BG_mScarletI;
    dataCollect.intVenus(end+1:end+size(bioInfo(i).intVenus,1),:) = bioInfo(i).intVenus - bioInfo(i).BG_Venus;
    dataCollect.intPVD(end+1:end+size(bioInfo(i).intPVD,1),:) = bioInfo(i).intPVD - bioInfo(i).BG_PVD;
    dataCollect.intCyOFP(end+1:end+size(bioInfo(i).intCyOFP,1),:) = bioInfo(i).intCyOFP - bioInfo(i).BG_CyOFP;
    dataCollect.intTDsmURFP(end+1:end+size(bioInfo(i).intTDsmURFP,1),:) = bioInfo(i).intTDsmURFP - bioInfo(i).BG_TDsmURFP;
    dataCollect.MajorAxisLength(end+1:end+size(bioInfo(i).MajorAxisLength,1),:) = bioInfo(i).MajorAxisLength;
    dataCollect.gR(end+1:end+size(bioInfo(i).gR,1),:) = bioInfo(i).gR;
    dataCollect.Centroid(end+1:end+size(bioInfo(i).Centroid,1),:) = bioInfo(i).Centroid;
end

% 删除为空的结构体dataCollect的field
fields = fieldnames(dataCollect);%结构体的fields of struct
for iFds = 1 :numel(fields)%iFds is short for Fields of struct
    if isempty(dataCollect.(fields{iFds}))
        dataCollect = rmfield(dataCollect,fields{iFds});
    end
end

fields = fieldnames(dataCollect);
% case for empty samples
if ~isfield(dataCollect,'Centroid')|| isempty(fields)
    return
end

% gR筛选
gRTF = ~ismissing(dataCollect.gR(:,1));% 找到非nan数据
TF = gRTF;
% posTF = (dataCollect.Centroid(:,1)>= xPosLimit(1)&dataCollect.Centroid(:,1)<= xPosLimit(2));
% if strcmp(fixedChannel,'BF1')
%     %     intTF = dataCollect.(fields{1})> intThreshold;%光强的筛选
%     intTF = true(size(posTF));
% else
%     intTF = dataCollect.(strcat('int',fixedChannel))> intThreshold;
% 
% end
% TF = posTF & gRTF & intTF;
% 
% for iFds = 1 : numel(fields)
%     if strlength(fields{iFds})<3% case for string less than 3 letters
%         continue
%     end
%     if strcmp(fields{iFds}(1:3),'int') %光强减去背景100
%         dataCollect.(fields{iFds}) =  dataCollect.(fields{iFds})-100;
%         % eachChIntTF = dataCollect.(fields{iFds})> 0;% 对所有荧光通道进行筛选 保证为正值
%         % TF = TF & eachChIntTF;
%     end
% end
% 
for iFds = 1 : numel(fields)
    dataCollect.(fields{iFds}) = dataCollect.(fields{iFds})(TF,:);
end

dataCollect = rmfield(dataCollect,'Centroid');

save(strcat(dirSave,'\dataCollect.mat'),'dataCollect','-v7.3');
%plot function
basicInfoPlot(dataCollect,scale,dirSave)
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
subplot1 = subplot(3,1,1,'Parent',f1);
for iFds = 1:channelNum %除去MajorAxisLength和gR fields,-2
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

subplot2 = subplot(3,1,2,'Parent',f1);
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

subplot3 = subplot(3,1,3,'Parent',f1);
meanValue = mean(dataCollect.gR(:,1));
stdValue = std(dataCollect.gR(:,1));
dataName = strcat('gR',32, num2str(meanValue,'%.3f'),'±',num2str(stdValue,'%.3f'),'min-1');
histogram(dataCollect.gR(:,1),'DisplayName',dataName,'FaceColor','k','EdgeColor','none')

set(subplot3,'FontSize',12);
box(subplot3,'off');
legend(subplot3,'show');
xlabel({'gR (min-1)'});
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
