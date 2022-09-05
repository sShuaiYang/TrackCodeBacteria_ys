function [ proteomeData ] = proteomeAnalysis()
dirFile  = 'C:\Users\XJY\Desktop\2020-08-质谱-色谱（LC-MS数据分析）-ys\cxy20200708-LC-MS数据分析';
data = readtable(strcat(dirFile,'\proteinGroups.txt'));
proteomeData =  table2struct(data);
% proteomeData filter
proteomeData =  proteomeDataFilter(proteomeData);
% protein gene name get
[proteomeData,geneSymbols] = proteomeDataGeneNameGet(proteomeData);
% geneOrderedLocusNamesGet
% proteomeData = geneOrderedLocusNamesGet(proteomeData);
save(strcat(dirFile,'\','proteomeData.mat'),'proteomeData');
writetable(struct2table(proteomeData),strcat(dirFile,'\','proteomeData.txt'),'Delimiter',' ') ;
writetable(struct2table(proteomeData),strcat(dirFile,'\','proteomeData.xlsx')) ;
writetable(struct2table(proteomeData),strcat(dirFile,'\','proteomeData.csv'));
save(strcat(dirFile,'\','geneSymbols.mat'),'geneSymbols');

% find LFQIntensity 进行处理
fields = fieldnames(proteomeData);
fldsTF = false(numel(fields),1);
for iFlds = 1:numel(fields)
    if strlength(fields{iFlds})>12 && strcmp(fields{iFlds}(1:12),'LFQIntensity')
        fldsTF(iFlds) = true;
    end
end
fields = fields(fldsTF);
%LFQInt数据获取
LFQInt = NaN(size(proteomeData,1),numel(fields));
for iProtein = 1:size(proteomeData,1)
    for iLabel = 1:numel(fields)
        if proteomeData(iProtein).(fields{iLabel}) ~= 0
            LFQInt(iProtein,iLabel) = proteomeData(iProtein).(fields{iLabel});
        end
    end
end
% 对LFQInt进行归一化处理
LFQInt_norm1 = LFQInt;% 用中位数median进行归一 Median normlizaiton
LFQInt_norm2 = LFQInt;% 用sum(xi)进行归一 Global normlizaiton
LFQInt_norm3 = LFQInt;% 用zscore进行归一 Z-score normlizaiton
LFQInt_norm4 = LFQInt;%  Max-min normlizaiton
for iSamp = 1: size(LFQInt,2)%sample 样品
    valueTF = ~isnan(LFQInt(:,iSamp));
    value_med = median(LFQInt(valueTF,iSamp));
    value_sum = sum(LFQInt(valueTF,iSamp));
    value_mean = mean(LFQInt(valueTF,iSamp));
    value_std = std(LFQInt(valueTF,iSamp));
    value_max = max(LFQInt(valueTF,iSamp));% 可以用类似的表达代替：max(A,[],'omitnan')
    value_min = min(LFQInt(valueTF,iSamp));
    LFQInt_norm1(:,iSamp) = LFQInt(:,iSamp)/value_med;
    LFQInt_norm2(:,iSamp) = LFQInt(:,iSamp)/value_sum;
    LFQInt_norm3(:,iSamp) = (LFQInt(:,iSamp)-value_mean)/value_std;
    LFQInt_norm4(:,iSamp) = (LFQInt(:,iSamp)-value_min)/(value_max-value_min);
end
save(strcat(dirFile,'\','LFQInt.mat'),'LFQInt','LFQInt_norm1','LFQInt_norm2','LFQInt_norm3','LFQInt_norm4');
label = {'No normalizaiton','Median normalizaiton','Global normalizaiton',...
    'Z-score normalizaiton','Max-min normalizaiton'};

smapleNum = size(LFQInt,2);
groupNum = 2;
volcanoDataGetAndPlot(LFQInt,proteomeData,fields,smapleNum,groupNum,dirFile,label{1});
volcanoDataGetAndPlot(LFQInt_norm1,proteomeData,fields,smapleNum,groupNum,dirFile,label{2});
volcanoDataGetAndPlot(LFQInt_norm2,proteomeData,fields,smapleNum,groupNum,dirFile,label{3});
% volcanoDataGetAndPlot(LFQInt_norm3,proteomeData,fields,smapleNum,groupNum,dirFile,label{4});
volcanoDataGetAndPlot(LFQInt_norm4,proteomeData,fields,smapleNum,groupNum,dirFile,label{5});
end


function volcanoDataGetAndPlot(A,proteomeData,fields,smapleNum, groupNum, dirFile, label)
group_Idx = false(groupNum,smapleNum);
for iGroup = 1:groupNum
    posIdx1 = (iGroup-1)*smapleNum/groupNum +1;
    posIdx2= iGroup*smapleNum/groupNum;
    group_Idx(iGroup,posIdx1:posIdx2) = true;
end
% groupNum = 4;
% group_Idx = [...
%     1,1,1,0,0,0,0,0,0,0,0,0;...
%     0,0,0,1,1,1,0,0,0,0,0,0;...
%     0,0,0,0,0,0,1,1,1,0,0,0;...
%     0,0,0,0,0,0,0,0,0,1,1,1;];
% group_Idx = logical(group_Idx);

dirSave = strcat(dirFile,'\',label);
if ~isfolder(dirSave)
    mkdir(dirSave)
end
% A = LFQInt_norm3;
groupData = struct();
for iGroup = 1:groupNum
    groupData(iGroup).data = A(:,group_Idx(iGroup,:));
    groupData(iGroup).meanValue = mean(groupData(iGroup).data, 2, 'omitnan');
end


for iGroup = 1:groupNum
    for iTest = 1:groupNum
        [~,p] = ttest2( groupData(iGroup).data', groupData(iTest).data');
        groupData(iTest).(strcat('pValue_x',num2str(iGroup))) = p';
    end
end
save(strcat(dirSave,'\','groupData.mat'),'groupData')

for iGroup = 1:groupNum-1
    for iTest = iGroup+1:groupNum
        fname = strcat('volcanoPlot_group',num2str(iGroup),'vs',num2str(iTest));
        f1 = figure('Name',fname);
        xData = log2(groupData(iGroup).meanValue./groupData(iTest).meanValue);
        yData = -log10( groupData(iGroup).(strcat('pValue_x',num2str(iTest))));
        line(xData,yData,'LineStyle','none','MarkerSize',4,'Marker','o')
        xlabel(['log2(group',num2str(iGroup),'/',num2str(iTest),')'])
        ylabel('-log10(pValue)')
        saveas(f1,strcat(dirSave,'\',fname,'.fig'));
        
        plotData.oriData1 = groupData(iGroup).data;
        plotData.oriData2 = groupData(iTest).data;
        plotData.oriData1_mean = groupData(iGroup).meanValue;
        plotData.oriData2_mean = groupData(iTest).meanValue;
        plotData.xData = xData;
        plotData.yData = yData;
        save(strcat(dirSave,'\','plotData_group',num2str(iGroup),'vs',num2str(iTest),'.mat'),'plotData');
    end
end
proteomeData_New = proteomeData;
for iProtein = 1:size(proteomeData,1)
    for iSamp = 1: smapleNum
        proteomeData_New(iProtein).(strcat('Normlized',fields{iSamp})) = ...
            A(iProtein,iSamp);
    end
end
save(strcat(dirSave,'\','proteomeData_New.mat'),'proteomeData_New')
proteomeData_New = struct2table(proteomeData_New);
writetable(proteomeData_New,strcat(dirSave,'\','proteomeData_New.txt'),'Delimiter',' ') ;
writetable(proteomeData_New,strcat(dirSave,'\','proteomeData_New.xlsx')) ;
writetable(proteomeData_New,strcat(dirSave,'\','proteomeData_New.csv')) ;
close all
end



