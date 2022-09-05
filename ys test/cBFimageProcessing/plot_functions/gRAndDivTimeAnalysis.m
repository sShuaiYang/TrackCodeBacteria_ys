function [growthDataAll,gRDiv] = gRAndDivTimeAnalysis(growthDataAll,dirSave)
% gRAndDivTimeAnalysis
% % dirAll='\\192.168.1.12\e\2019-12-23 PAO1_IP32_100x ys';
scale=growthDataAll{1, 1}.scale;
% dirSave=strcat(dirAll,'\growthAnalysis');
% mkdir(dirSave);
gRDiv=[];
%��1�� field
%��2�� ����
%��3�� ��������
%��4�� ����ʱ����
%��5�� ����ʱ���
%��6�� ����ʱ���
%��7�� ������
%��8�� division time

n=1;
for iField=1:numel(growthDataAll)
    if isempty(growthDataAll{iField}) % �жϴ���Ұ�������Ƿ�Ϊ��%�е���Ұ����û�д���
        continue
    else
        if numel(growthDataAll{iField}.growthInfo)< 2 % �жϴ���Ұ�������Ƿ�Ϊ��
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





