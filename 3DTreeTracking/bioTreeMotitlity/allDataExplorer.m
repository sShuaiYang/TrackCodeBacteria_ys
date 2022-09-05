function allDataExplorer(allData,dirSave)
% dirSave=uigetdir();
cd(dirSave)
for iData=1:size(allData.dataInfo,2)
    mkdir(dirSave,'Graphic')
    saveGraphic=strcat(dirSave,'\Graphic');
% iData=4;
    if isempty(allData.dataInfo{iData}.branchIndex)
        rootInfo=allData.dataInfo{iData}.rootInfo;
        folderName=strcat(num2str(allData.dataList{iData}(1)),'-',num2str(allData.dataList{iData}(2)),'__root',num2str(rootInfo(1)),'_',num2str(rootInfo(2)));
        mkdir(folderName);
        saveFile1=strcat(dirSave,'\',folderName);
        saveGraphic=strcat(saveGraphic,'\',strcat('Graphic___',folderName));
        timeResult=timeSeries(allData.dataInfo{iData});
        histResult=histVelocity(allData.dataInfo{iData}.velocityData);
        plotResultsAndSave(saveGraphic,saveFile1,timeResult,histResult,allData.dataList{iData},iData);
        close all;
    elseif  ~isempty(allData.dataInfo{iData}.branchIndex)
        branchIndex=allData.dataInfo{iData}.branchIndex;
        folderName=strcat(num2str(allData.dataList{iData}(1)),'-',num2str(allData.dataList{iData}(2)),'_branch',num2str(branchIndex));
        mkdir(folderName);
        cd(folderName);
        if allData.dataInfo{iData}.isNode==false
            rootInfo=allData.dataInfo{iData}.rootInfo;
            folderName1=strcat(num2str(allData.dataList{iData}(1)),'-',num2str(allData.dataList{iData}(2)),'_root',num2str(rootInfo(1)),'_',num2str(rootInfo(2)));
            mkdir(folderName1);
            saveFile2=strcat(dirSave,'\',folderName,'\',folderName1);
            saveGraphic=strcat(saveGraphic,'\',strcat('Graphic___',folderName));
            timeResult=timeSeries(allData.dataInfo{iData});
            histResult=histVelocity(allData.dataInfo{iData}.velocityData);
            plotResultsAndSave(saveGraphic,saveFile2,timeResult,histResult,allData.dataList{iData},iData);
            close all;
        elseif allData.dataInfo{iData}.isNode==true
            nodeInfo=allData.dataInfo{iData}.nodeInfo;
            folderName1=strcat(num2str(allData.dataList{iData}(1)),'-',num2str(allData.dataList{iData}(2)),'_node',num2str(nodeInfo(1)),'_',num2str(nodeInfo(2)),'_',num2str(nodeInfo(3)));
            mkdir(folderName1);
            saveFile2=strcat(dirSave,'\',folderName,'\',folderName1);
            saveGraphic=strcat(saveGraphic,'\',strcat('Graphic___',folderName));
            timeResult=timeSeries(allData.dataInfo{iData});
            histResult=histVelocity(allData.dataInfo{iData}.velocityData);
            plotResultsAndSave(saveGraphic,saveFile2,timeResult,histResult,allData.dataList{iData},iData);
            close all;
        end
        cd('..');
    end
end
end