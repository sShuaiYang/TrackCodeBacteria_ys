function timeLapseMultiFieldFluoImagePlot_singleField(dirFile,fixedChannel)
%timeLapseMontageFluoImagePlot 用于将单个视野的荧光强度变化作图
% shuai Yang 2020.08.04
% 未经过串光校正的荧光强度作图
% update 将每个通道的光强进行筛选要求>0 去掉,2021/10/13
% 只要求fixedChannel的光强>intThreshold,2021/10/13
% 取消了位置的筛选和光强的筛选 2021/10/28

% scale =0.065;%100x pixel to um
% intThreshold = 105;%没有减背景 default105 最小为100  intensity Threshold
% if strcmp(fixedChannel,'BF1')
%     %     xPosLimit = [244,1860];
%     %BF的mask，荧光图像的视野范围小
%     xPosLimit = [0,2048];%2020.12.30 视野照明范围变大了 不用再截了
% else
%     xPosLimit = [0,2048];
% end

dirSave =strcat(dirFile,'\result_basic');
if ~isfolder(dirSave)
    mkdir(dirSave);
end

disp('Single field fluo intensity change plot')

fieldList = dir([dirFile,filesep,'field*']);
fieldList = fieldListClean (fieldList);%只保留field的文件和两个系统文件
for iField= 1:(length(fieldList))
    if ~strcmp(fieldList(iField).name(1:5),'field')
        continue
    end
    disp(fieldList(iField).name);
    dirField = strcat(dirFile,'\',fieldList(iField).name);


    if isfile(strcat(dirField,'\','bioInfo.mat'))
        load([dirField,'\bioInfo.mat']);
    end
    dataCollect_tps = cell(1,numel(bioInfo));

    for iTime = 1 : numel(bioInfo) % 每个时间点的数据获得
        dataCollect = struct('intsfGFP',[ ],'intmScarletI',[ ],'intVenus',[ ],...
            'intPVD',[ ],'intCyOFP',[ ],'intTDsmURFP',[ ],'MajorAxisLength',[ ],...
            'Centroid',[ ]);
        % here dataCollect represents for timepoint datacollect

        dataCollect.intsfGFP(end+1:end+size(bioInfo(iTime).intsfGFP,1),:) = bioInfo(iTime).intsfGFP - bioInfo(iTime).BG_sfGFP;
        dataCollect.intmScarletI(end+1:end+size(bioInfo(iTime).intmScarletI,1),:) = bioInfo(iTime).intmScarletI - bioInfo(iTime).BG_mScarletI;
        dataCollect.intVenus(end+1:end+size(bioInfo(iTime).intVenus,1),:) = bioInfo(iTime).intVenus - bioInfo(iTime).BG_Venus;
        dataCollect.intPVD(end+1:end+size(bioInfo(iTime).intPVD,1),:) = bioInfo(iTime).intPVD - bioInfo(iTime).BG_PVD;
        dataCollect.intCyOFP(end+1:end+size(bioInfo(iTime).intCyOFP,1),:) = bioInfo(iTime).intCyOFP - bioInfo(iTime).BG_CyOFP;
        dataCollect.intTDsmURFP(end+1:end+size(bioInfo(iTime).intTDsmURFP,1),:) = bioInfo(iTime).intTDsmURFP -bioInfo(iTime).BG_TDsmURFP;
        dataCollect.MajorAxisLength(end+1:end+size(bioInfo(iTime).MajorAxisLength,1),:) = bioInfo(iTime).MajorAxisLength;
        dataCollect.Centroid(end+1:end+size(bioInfo(iTime).Centroid,1),:) = bioInfo(iTime).Centroid;

        % 通过Centroid判断是否有数据点 如果为空说明此时间点没有细菌的数据
        if isempty(dataCollect.Centroid)
            continue
        end

        % 删除为空的结构体dataCollect的field
        fields = fieldnames(dataCollect);%结构体的fields of struct
        for iFds = 1 : numel(fields)%iFds is short for Fields of struct
            if isempty(dataCollect.(fields{iFds}))
                dataCollect = rmfield(dataCollect,fields{iFds});
            end
        end

%         fields = fieldnames(dataCollect);
        %通过'Centroid'判断dataCollect是否为空
        if isfield(dataCollect,'Centroid')
            %             posTF = (dataCollect.Centroid(:,1)>= xPosLimit(1)&dataCollect.Centroid(:,1)<= xPosLimit(2));
            %             if strcmp(fixedChannel,'BF1')
            %                 % intTF = dataCollect.(fields{1})> intThreshold;%光强的筛选
            %                 intTF = true(size(posTF));
            %             else
            %                 intTF = dataCollect.(strcat('int',fixedChannel))> intThreshold;
            %             end
            %             TF = posTF & intTF;
            %
            %             for iFds = 1 : numel(fields)
            %                 if strlength(fields{iFds})<3% case for string less than 3 letters
            %                     continue
            %                 end
            %                 if strcmp(fields{iFds}(1:3),'int') %光强减去背景100
            %                     dataCollect.(fields{iFds}) =  dataCollect.(fields{iFds})-100;
            %                     % eachChIntTF = dataCollect.(fields{iFds})> 0;% 对所有荧光通道进行筛选 保证为正值
            %                     % TF = TF & eachChIntTF;
            %                 end
            %             end
            %             for iFds = 1 : numel(fields)
            %                 dataCollect.(fields{iFds}) = dataCollect.(fields{iFds})(TF,:);
            %             end

            dataCollect = rmfield(dataCollect,'Centroid');
        end
        dataCollect_tps{iTime} = dataCollect;
    end

    intensityChangePlot_singleField(dataCollect_tps,dirField,dirSave);
end

end
%%
function  intensityChangePlot_singleField(dataCollect_tps,dirField,dirSave)
%case for no cells in this timepoint 有可能某个时间点数据为空
for iTime = 1:numel(dataCollect_tps)
    if ~isempty(dataCollect_tps{iTime})
        fields = fieldnames(dataCollect_tps{iTime});
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
if channelNum == 0
    return
end

load([dirField,'\Tracking\frameInfo.mat']);
frameInfo = frameInfo(1:numel(dataCollect_tps),:);
time_lag = abs(etime (frameInfo(1,1:6),frameInfo(:,1:6))/60);

fname = strcat('intenisytChange','_',dirField(end-8:end));
f1 = figure('Name',fname);
for iFds = 1:channelNum
    fluoInt = zeros(numel(dataCollect_tps),2);
    for iTime = 1:numel(dataCollect_tps)
        try
            fluoInt(iTime,1) = mean(dataCollect_tps{iTime}.(fields{iFds}));
            fluoInt(iTime,2) = std(dataCollect_tps{iTime}.(fields{iFds}));
        catch
            msg = ['No cell was dectected in frame ',num2str(iTime),';',...
                ' Specify the ',fields{iFds}(4:end),' values as NaN.'];
            warning(msg);
            fluoInt(iTime,1) = NaN;
            fluoInt(iTime,2) = NaN;
        end
    end
    subplot1 = subplot(channelNum,1,iFds,'Parent',f1);
    %     errorbar(1:iTime,fluoInt(:,1),fluoInt(:,2),'o')
    %         line(1:iTime,fluoInt(:,1),'Color',[0,0.45,0.74],'LineWidth',2,'Marker','o' );
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
function fieldList = fieldListClean (fieldList)
% 只保留field 文件和两个系统文件
templogic = false (1, numel (fieldList));
for iField = 1: numel (fieldList)

    if(isequal(fieldList(iField).name,'.')||... % 系统自带的两个隐文件夹
            isequal(fieldList(iField).name,'..')||...
            strcmp(fieldList(iField).name(1:5),'field'))
        templogic(iField) = true;
    end

end
fieldList = fieldList(templogic);
end