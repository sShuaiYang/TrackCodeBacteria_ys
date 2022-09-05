function reCalculateGrowthDataMultiField()
%用于根据bioTree重新对多视野板底实验进行计算
%用函数generationGrowthRateCal_agrose_new重新对各个视野的数据进行计算

dirFile='\\192.168.1.14\e\2019-12-26-deltpslpelfliC_IP31 ys';
fieldList=dir(dirFile);
parfor iField = 1:length(fieldList)-2
    if ~(strcmp(fieldList(iField+2).name(1:5),'field') )
        continue
    end
    dirField=[dirFile,'\',fieldList(iField+2).name];
    if ~isfile(strcat(dirField,'\bacInfo&tree.mat'))
        continue
    end
    temp = load (strcat(dirField,'\bacInfo&tree.mat'));
    bioTree = temp.bioTree;
    
    
%     frameInfoNew =zeros(numel(bioTree),size(frameInfo,2));
%     frameInfoNew(1:2:numel(bioTree)-1,:)=frameInfo(1:numel(bioTree)/2,:);
%     frameInfoNew(2:2:numel(bioTree),:)=frameInfo(1:numel(bioTree)/2,:);
%     bioTreeInitialTime=frameInfoNew(1,1:6);
%     bioTreeTimer =zeros(1,numel(bioTree));
%     for i=1:numel(bioTree)
%         bioTreeTimer(i)=etime(frameInfoNew(i,1:6),bioTreeInitialTime)/60;   %min
%     end

    bioTree{1}.bioTreeTimer=bioTreeTimer;
    [~,~]=generationGrowthRateCal_agrose_new(bioTree,dirField);
end
end

