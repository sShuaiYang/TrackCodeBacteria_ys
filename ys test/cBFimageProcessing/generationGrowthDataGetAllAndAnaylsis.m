function [growthDataAll] = generationGrowthDataGetAllAndAnaylsis()
% generationGrowthDataGetAllAndAnaylsis
% allCellsGenerationGrowthRateGetFromOneAgroseExperiment
dirAll='F:\2020-05-09-PAO1_IP31_100x 15agar chip Temp30 ys';
nameList = dir(dirAll);
% fieldNum=length(nameList)-2;
fieldNum = 64;
growthDataAll = cell(1,fieldNum);
tempFile = cell(1,fieldNum);
for iField = 1:fieldNum%可以用parfor 如果内存够用
    dirField=[dirAll,'\',nameList(iField+2).name];
    %     temp=load([dirFile,'\','bacInfo&tree']);
    if isfile(strcat(dirField,'\growthResult2','\','growthData.mat'))
        tempFile{iField} = load([dirField,'\growthResult2','\','growthData']);
        growthDataAll{iField} = tempFile{iField}.growthData;
    end
    %     [growthData,~]=generationGrowthRateCal_agrose_new(temp.bioTree,dirFile)
%     growthDataAll{iField} = tempFile{iField}.growthData;
end
clear tempFile

gRGetAll_Num=numel(growthDataAll{1}.gRGet);
gR_kinNum=numel(growthDataAll{1}.gR_kin);
for iField=1:numel(growthDataAll)
    if isempty(growthDataAll{iField})
        continue
    end
    if gRGetAll_Num<numel(growthDataAll{iField}.gRGet)
        gRGetAll_Num=numel(growthDataAll{iField}.gRGet);
    end
    if gR_kinNum<numel(growthDataAll{iField}.gR_kin)
        gR_kinNum=numel(growthDataAll{iField}.gR_kin);
    end    
end

gRGetAll=cell(1,gRGetAll_Num);
gR_kinAll=cell(1,gR_kinNum);
for iField=1:numel(growthDataAll)
    if isempty(growthDataAll{iField})
        continue
    end
    for j=1:numel(growthDataAll{iField}.gRGet )
        gRGetAll{j}=[gRGetAll{j},growthDataAll{iField}.gRGet{j}];
    end
    for j=1:numel(growthDataAll{iField}.gR_kin)
        gR_kinAll{j}=[gR_kinAll{j};growthDataAll{iField}.gR_kin{j}];
    end
end


% 
load([dirAll,'\mip.mat']);
scale=mip.Calib.scale;
% scale=0.0650;
xyPosition=mip.fieldTag.xyPosition;
fieldPosition=xyPosition-xyPosition(1,:);% micromanager code V3.0 code %unit um
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
            growthDataAll{iField}.growthInfo{iGenr}(iCell).absCentroid=...
                growthDataAll{iField}.growthInfo{iGenr}(iCell).Centroid*scale+fieldPosition(fieldTag,:);
        end
    end
end
growthDataAll{1}.scale=scale;
growthDataAll{1}.fieldPosition=fieldPosition;
growthDataAll{1}.xyPosition=xyPosition;

growthDataAll{1}.gRGetAll=gRGetAll;
growthDataAll{1}.gR_kinAll=gR_kinAll;
dirSave=strcat(dirAll,'\growthAnalysis');
mkdir(dirSave);
[growthDataAll] = allGenerationgRPlot(growthDataAll,dirSave);
% save([dirAll,'\growthDataAll.mat'],'growthDataAll','-v7.3');
growthDataAllSave(growthDataAll,dirSave);
end
%%
function growthDataAllSave(growthDataAll,dirFile)

dirSave=strcat(dirFile,'\growthDataAll');
mkdir(dirSave);

fieldNum=numel(growthDataAll);
saveNum=sqrt(fieldNum);
for i=0:saveNum
    
    if i==0
        growthDataSub=growthDataAll(1);
        save([dirSave,'\',num2str( i,'%04.f')],'growthDataSub','-v7.3');
    else
        if i==1
            growthDataSub=growthDataAll(2:saveNum*i);
            save([dirSave,'\',num2str( i,'%04.f')],'growthDataSub','-v7.3');
        else
            growthDataSub=growthDataAll(saveNum*(i-1)+1:saveNum*i);
            save([dirSave,'\',num2str( i,'%04.f')],'growthDataSub','-v7.3')
        end
    end
    
end

end
%% growthDataAllLoad
function  growthDataAll=growthDataAllLoad()
dirFile='E:\2019-12-23 PAO1_IP32_100x ys\growthDataAll';
nameList=dir(dirFile);
for i=1:numel(nameList)-2
    load( [dirFile,'\',nameList(i+2).name]);
    if i==1
        growthDataAll=growthDataSub;
    else
        growthDataAll(end+1:end+numel(growthDataSub))=growthDataSub;
    end       
end

end

