function [growthDataAll,gRPos_class] = gRAndDistanceAnalysis(growthDataAll,dirSave)
% dirSave='\\192.168.1.12\e\2019-12-23 PAO1_IP32_100x ys';
% load([dirSaveFile,'\mip.mat']);
% % scale=mip.Calib.scale;
% 
% scale=0.0650;
% xyPosition=mip.fieldTag.xyPosition;
% fieldPosition=xyPosition-xyPosition(1,:);
% 
% growthDataAll{1}.scale=scale;
% growthDataAll{1}.xyPosition=xyPosition;
% growthDataAll{1}.fieldPosition=fieldPosition;
% % 
% for iField=1:numel(growthDataAll)
%     fieldTag=growthDataAll{iField}.growthInfo{1}.fieldTag;
%     for iGenr=1:numel(growthDataAll{iField}.growthInfo)
%         for iCell=1:numel(growthDataAll{iField}.growthInfo{iGenr})
%             growthDataAll{iField}.growthInfo{iGenr}(iCell).absCentroid=...
%                 growthDataAll{iField}.growthInfo{iGenr}(iCell).Centroid*scale+fieldPosition(fieldTag,:);
%         end
%         
%     end
%     
% end

k=1;
for iField=1:numel(growthDataAll)
    if isempty(growthDataAll{iField}) % 判断此视野的数据是否为空%有的视野数据没有处理
        continue
    else
        if numel(growthDataAll{iField}.growthInfo)== 0 % 判断此视野的数据是否为空
            continue
        end
    end
    fieldTag=growthDataAll{iField}.growthInfo{1}.fieldTag;
    for iGenr=1:numel(growthDataAll{iField}.growthInfo)
        for iCell=1:numel(growthDataAll{iField}.growthInfo{iGenr})
            gRPos(k,1)=fieldTag;%field
            gRPos(k,2)=iGenr;%generation
            gRPos(k,3)=growthDataAll{iField}.growthInfo{iGenr}(iCell).gR;%growth rate
            gRPos(k,4:5)=growthDataAll{iField}.growthInfo{iGenr}(iCell).absCentroid(1,:);% position
            k=k+1;
        end
        
    end
    
end


gRPos_class=cell(1,max(gRPos(:,2)));
% parfor i=1:max(gRPos(:,2))% 可以用parfor 如果内存不会爆的话
for i=1:max(gRPos(:,2))
    templogic=gRPos(:,2)==i;
    gRPos_class{i}.gRPos=gRPos(templogic,:);
    if size(gRPos_class{i}.gRPos,1)>=2
        cellNum=size(gRPos_class{i}.gRPos,1);
        D= pdist(gRPos_class{i}.gRPos(:,4:5));
        gRPos_class{i}.D=D;
        gRPos_class{i}.Z=squareform(D);
        plotData=zeros((cellNum*(cellNum-1))/2,4);
        n=1;
        for ix=1:size(gRPos_class{i}.gRPos,1)-1
            for iy=ix+1:size(gRPos_class{i}.gRPos,1)
                plotData(n,1)=gRPos_class{i}.gRPos(ix,3);%g(r1)
                plotData(n,2)=gRPos_class{i}.gRPos(iy,3);%g(r2)
                plotData(n,3)=plotData(n,1)*plotData(n,2);%g(r1)*g(r2)
                plotData(n,4)=gRPos_class{i}.Z(ix,iy);%dist(r1,r2)
                n=n+1;
            end
        end
        gRPos_class{i}.plotData=plotData;
    end
end
save([dirSave,'\gRPos.mat'],'gRPos','gRPos_class','-v7.3');
dirFigSave = strcat(dirSave,'\gRAndDistance');
mkdir(dirFigSave);
gRAndDistancePlot(gRPos_class,dirFigSave);
end
