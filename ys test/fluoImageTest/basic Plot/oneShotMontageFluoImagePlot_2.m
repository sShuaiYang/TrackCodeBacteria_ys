function [dataCollect] = oneShotMontageFluoImagePlot_2(dirFile,bioInfo,fixedChannel)
% 用于对经过[bioInfo] = cellFluoIntensityCorrection(bioInfo)光强进行校正后的作图
% bioInfo 中有bioInfo结构体中有类似intsfGFP_c的field
% 经过串色后的荧光强度作图 intsfGFP_c 已经减过背景 已经串光校正
% allCorrectedFluoName = {'intsfGFP_c','intVenus_c',...
%     'intCyOFP_c','intmScarletI_c','intTDsmURFP_c'};% c represents for corrected
% shuai Yang 2020.06.09
% update 将每个通道的光强进行筛选要求>0 去掉, 2021/10/13
% 只要求fixedChannel的光强>intThreshold,2021/10/13
% 取消了位置的筛选和光强的筛选 2021/10/28

allFluoName_c = {'intsfGFP','intVenus','intCyOFP',...
    'intmScarletI','intTDsmURFP'};%校正的荧光通道
% scale 获得
try
    load([dirFile,filesep,'mip.mat']);
    scale = mip.Calib.scale ; 
catch
    scale = 0.065;%100x pixel to um %default
end

% intThreshold = 5;%已经减过背景 default 5 最小为0  intensity Threshold
% 
% if strcmp(fixedChannel,'BF1')
%     %     xPosLimit = [300,1800];
%     %BF的mask，荧光图像的视野范围小
%     xPosLimit = [0,2048];%2020.12.30 视野照明范围变大了 不用再截了
% else
%     xPosLimit = [0,2048];
% end

dirSave = strcat(dirFile,'\result_basic2');
mkdir(dirSave);
dataCollect = struct('intsfGFP',[ ],'intmScarletI',[ ],'intVenus',[ ],...
    'intPVD',[ ],'intCyOFP',[ ],'intTDsmURFP',[ ],'MajorAxisLength',[ ],...
    'MinorAxisLength',[ ],'Centroid',[ ]);

for i = 1:numel(bioInfo)
    dataCollect.intsfGFP(end+1:end+size(bioInfo(i).intsfGFP,1),:) = bioInfo(i).intsfGFP_c ;
    dataCollect.intmScarletI(end+1:end+size(bioInfo(i).intmScarletI,1),:) = bioInfo(i).intmScarletI_c;
    dataCollect.intVenus(end+1:end+size(bioInfo(i).intVenus,1),:) = bioInfo(i).intVenus_c;
    dataCollect.intPVD(end+1:end+size(bioInfo(i).intPVD,1),:) = bioInfo(i).intPVD - bioInfo(i).BG_PVD;%PVD无校正的数据 直接减背景
    dataCollect.intCyOFP(end+1:end+size(bioInfo(i).intCyOFP,1),:) = bioInfo(i).intCyOFP_c;
    dataCollect.intTDsmURFP(end+1:end+size(bioInfo(i).intTDsmURFP,1),:) = bioInfo(i).intTDsmURFP_c;
    dataCollect.MajorAxisLength(end+1:end+size(bioInfo(i).MajorAxisLength,1),:) = bioInfo(i).MajorAxisLength;
    dataCollect.MinorAxisLength(end+1:end+size(bioInfo(i).MinorAxisLength,1),:) = bioInfo(i).MinorAxisLength;
    dataCollect.Centroid(end+1:end+size(bioInfo(i).Centroid,1),:) = bioInfo(i).Centroid;
end

% 删除为空的结构体dataCollect的field
fields = fieldnames(dataCollect);%结构体的fields of struct
for iFds = 1 :numel(fields)%iFds is short for iFields of struct
    if isempty(dataCollect.(fields{iFds}))
        dataCollect = rmfield(dataCollect,fields{iFds});
    end
end
fields = fieldnames(dataCollect);

% case for empty samples 
if isempty(fields)
    save(strcat(dirSave,'\dataCollect.mat'),'dataCollect','-v7.3');
    save(strcat(dirSave,'\bioInfo.mat'),'bioInfo','-v7.3');
    return
end

% posTF = (dataCollect.Centroid(:,1)>= xPosLimit(1)&dataCollect.Centroid(:,1)<= xPosLimit(2));
% if strcmp(fixedChannel,'BF1')
%     %     intTF = dataCollect.(fields{1})> intThreshold;%光强的筛选
%     intTF = true(size(posTF));
% else
%     intTF = dataCollect.(strcat('int',fixedChannel))> intThreshold;
%     
% end
% TF = posTF & intTF;
% 
% for iFds = 1 : numel(fields)
%     if strlength(fields{iFds})<3% case for string less than 3 letters
%         continue
%     end
%     if strcmp(fields{iFds}(1:3),'int')
%         if ~ismember(fields{iFds},allFluoName_c)
%             %没有校正的荧光通道光强减去100，其他在cellFluoIntensityCorrection已经减过
%             dataCollect.(fields{iFds}) =  dataCollect.(fields{iFds})-100;
%         end
%         % 对所有荧光通道进行筛选 保证为正值
%         %         eachChIntTF = dataCollect.(fields{iFds})> 0;
%         %         TF = TF & eachChIntTF;
%     end
% end
% 
% for iFds = 1 : numel(fields)
%     dataCollect.(fields{iFds}) = dataCollect.(fields{iFds})(TF,:);
% end

dataCollect = rmfield(dataCollect,'Centroid');
    
save(strcat(dirSave,'\dataCollect.mat'),'dataCollect','-v7.3');
save(strcat(dirSave,'\bioInfo.mat'),'bioInfo','-v7.3');
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
subplot1 = subplot(2,1,1,'Parent',f1);
for iFds = 1:channelNum
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
xlabel(subplot1,{'intensity (a.u.)'});
ylabel(subplot1,{'cell Num (a.u.)'});

subplot2 = subplot(2,1,2,'Parent',f1);
meanValue = mean(dataCollect.MajorAxisLength)*scale;
stdValue = std(dataCollect.MajorAxisLength)*scale;
dataName = strcat('cellLength',32, num2str(meanValue,'%.1f'),'±',num2str(stdValue,'%.1f'),'um');
histogram(dataCollect.MajorAxisLength*scale,'DisplayName',dataName,'FaceColor','k','EdgeColor','none')

hold(subplot2,'on')
meanValue = mean(dataCollect.MinorAxisLength)*scale;
stdValue = std(dataCollect.MinorAxisLength)*scale;
dataName = strcat('cellWidth',32, num2str(meanValue,'%.1f'),'±',num2str(stdValue,'%.1f'),'um');
histogram(dataCollect.MinorAxisLength*scale,'DisplayName',dataName,'FaceColor',"auto",'EdgeColor','none')

set(subplot2,'FontSize',12);
box(subplot2,'off');
legend(subplot2,'show');
xlabel(subplot2,{'cell Size (um)'});
ylabel(subplot2,{'cell Num (a.u.)'});
xlim(subplot2,[0,8])
hold(subplot2,'off')
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
                scatter(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),20,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.05);
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
                scatter(dataCollect.(fields{iCh}),dataCollect.(fields{jCh}),20,[0,0.45,0.74],'filled','MarkerFaceAlpha',0.05)
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