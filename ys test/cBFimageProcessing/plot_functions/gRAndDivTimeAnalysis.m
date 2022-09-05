function [growthDataAll,gRDiv] = gRAndDivTimeAnalysis(growthDataAll,dirSave)
% gRAndDivTimeAnalysis
% % dirAll='\\192.168.1.12\e\2019-12-23 PAO1_IP32_100x ys';
scale=growthDataAll{1, 1}.scale;
% dirSave=strcat(dirAll,'\growthAnalysis');
% mkdir(dirSave);
gRDiv=[];
%第1列 field
%第2列 代数
%第3列 出生菌长
%第4列 分裂时菌长
%第5列 出生时面积
%第6列 分裂时面积
%第7列 生长率
%第8列 division time

n=1;
for iField=1:numel(growthDataAll)
    if isempty(growthDataAll{iField}) % 判断此视野的数据是否为空%有的视野数据没有处理
        continue
    else
        if numel(growthDataAll{iField}.growthInfo)< 2 % 判断此视野的数据是否为空
            continue
        end
    end
    for iGenr=2:numel(growthDataAll{iField}.growthInfo)
        for iCell=1:numel(growthDataAll{iField}.growthInfo{iGenr})
            if growthDataAll{iField}.growthInfo{iGenr}(iCell).is2Node==1
                gRDiv(n,1)=growthDataAll{iField}.growthInfo{iGenr}(iCell).fieldTag;
                gRDiv(n,2)=growthDataAll{iField}.growthInfo{iGenr}(iCell).genrTag;
                gRDiv(n,3)=growthDataAll{iField}.growthInfo{iGenr}(iCell).MajorAxisLength(1);
                gRDiv(n,4)=growthDataAll{iField}.growthInfo{iGenr}(iCell).MajorAxisLength(end);
                gRDiv(n,5)=growthDataAll{iField}.growthInfo{iGenr}(iCell).cellArea(1);
                gRDiv(n,6)=growthDataAll{iField}.growthInfo{iGenr}(iCell).cellArea(end);
                gRDiv(n,7)=growthDataAll{iField}.growthInfo{iGenr}(iCell).gR;
                gRDiv(n,8)=growthDataAll{iField}.growthInfo{iGenr}(iCell).divTime;
                n=n+1;
            end
            
        end
        
    end
    
end
save([dirSave,'\gRDiv.mat'],'gRDiv');

dirFigSave = strcat(dirSave,'\gRAndDivTime');
mkdir(dirFigSave);
gRAndDivTimePlot(gRDiv,scale,dirFigSave);

end





