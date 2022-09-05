function getBioTreeAndMakeDemo(stepFrame)
disp('backGrond correction')
dirAll=superBatchConvert();
% dirAll='D:\Original deta  NiLei\2012-09-12 bacteria fast blank 30degree\5678';
nameList=dir(dirAll);
cd(dirAll);
for i=1:(length(nameList)-2)
    disp(strcat('AllJob',num2str(i)));
    disp(nameList(i+2).name);
    dirFile=strcat(dirAll,'\',nameList(i+2).name);
    bioTree=batchTreeTrackingJob(dirFile);
    [bioTree,allDataOri]=runaLL(bioTree);
    save(strcat(dirFile,'\allDataOri'),'allDataOri');
    startFrame=1;
    endFrame=size(bioTree,2);
    motilityDemoMaker(allDataOri,startFrame,stepFrame,endFrame,dirFile);
    [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri);
    disp('rockFrame is below: ');disp(rockFrame);
    correctionList=autoGetCorrectionList(bioTree,allDataOri);
    allData=dataCorrection(allDataOri,correctionList,[],[],dirFile,[]); 
%     load(strcat(dirFile,'\allData'))
    dirResultSave=strcat(dirFile,'\allResult');
    mkdir(dirResultSave);
    allDataExplorer(allData,dirResultSave);
end
end





% % 不做背景校正，从生成bioTree开始
% function getBioTreeAndMakeDemo(stepFrame)
% % disp('backGrond correction')
% % dirAll=superBatchConvert();
% dirAll='总文件夹名';
% nameList=dir(dirAll);
% cd(dirAll);
% for i=1:(length(nameList)-2)
%     disp(strcat('AllJob',num2str(i)));
%     disp(nameList(i+2).name);
%     dirFile=strcat(dirAll,'\',nameList(i+2).name);
%     bioTree=batchTreeTrackingJob(dirFile);
%     [bioTree,allDataOri]=runaLL(bioTree);
%     save(strcat(dirFile,'\allDataOri'),'allDataOri');
%     startFrame=1;
%     endFrame=size(bioTree,2);
%     motilityDemoMaker(allDataOri,startFrame,stepFrame,endFrame,dirFile);
%     [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri);
%     disp('rockFrame is below: ');disp(rockFrame);
%     correctionList=autoGetCorrectionList(bioTree,allDataOri);
%     allData=dataCorrection(allDataOri,correctionList,[],[],dirFile,rockFrame); 
% %     load(strcat(dirFile,'\allData'))
%     dirResultSave=strcat(dirFile,'\allResult');
%     mkdir(dirResultSave);
%     allDataExplorer(allData,dirResultSave);
% end
% end
% 
% 
% 
% % 不做bioTree，从runAll开始
% function getBioTreeAndMakeDemo(stepFrame)
% % disp('backGrond correction')
% % dirAll=superBatchConvert();
% dirAll='总文件夹名';
% nameList=dir(dirAll);
% cd(dirAll);
% for i=1:(length(nameList)-2)
%     disp(strcat('AllJob',num2str(i)));
%     disp(nameList(i+2).name);
%     dirFile=strcat(dirAll,'\',nameList(i+2).name);
%     bioTree=load(strcat(dirFile,'\bioTreeResult\allTree\bioTreeStack'));
%     bioTree=bioTree.bioTree;
% %     bioTree=batchTreeTrackingJob(dirFile);
%     [bioTree,allDataOri]=runaLL(bioTree);
%     save(strcat(dirFile,'\allDataOri'),'allDataOri');
%     startFrame=1;
%     endFrame=size(bioTree,2);
%     motilityDemoMaker(allDataOri,startFrame,stepFrame,endFrame,dirFile);
%     [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri);
%     disp('rockFrame is below: ');disp(rockFrame);
%     correctionList=autoGetCorrectionList(bioTree,allDataOri);
%     allData=dataCorrection(allDataOri,correctionList,[],[],dirFile,rockFrame); 
% %     load(strcat(dirFile,'\allData'))
%     dirResultSave=strcat(dirFile,'\allResult');
%     mkdir(dirResultSave);
%     allDataExplorer(allData,dirResultSave);
% end
% end
% 
% 
% % 从demo生成开始
% function getBioTreeAndMakeDemo(stepFrame)
% % disp('backGrond correction')
% % dirAll=superBatchConvert();
% dirAll='总文件夹名';
% nameList=dir(dirAll);
% cd(dirAll);
% for i=1:(length(nameList)-2)
%     disp(strcat('AllJob',num2str(i)));
%     disp(nameList(i+2).name);
%     dirFile=strcat(dirAll,'\',nameList(i+2).name);
%     bioTree=load(strcat(dirFile,'\bioTreeResult\allTree\bioTreeStack'));
%     bioTree=bioTree.bioTree;
%     %     bioTree=batchTreeTrackingJob(dirFile);
%     %     [bioTree,allDataOri]=runaLL(bioTree);
%     %     save(strcat(dirFile,'\allDataOri'),'allDataOri');
%     allDataOri=load(strcat(dirFile,'\allDataOri'));
%     allDataOri=allDataOri.allDataOri;
%     startFrame=1;
%     endFrame=size(bioTree,2);
%     motilityDemoMaker(allDataOri,startFrame,stepFrame,endFrame,dirFile);
%     [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri);
%     disp('rockFrame is below: ');disp(rockFrame);
%     correctionList=autoGetCorrectionList(bioTree,allDataOri);
%     allData=dataCorrection(allDataOri,correctionList,[],[],dirFile,rockFrame);
%     %     load(strcat(dirFile,'\allData'))
%     dirResultSave=strcat(dirFile,'\allResult');
%     mkdir(dirResultSave);
%     allDataExplorer(allData,dirResultSave);
% end
% end
% 
% 
% 
% 
% % 从做correctionList开始
% function getBioTreeAndMakeDemo(stepFrame)
% % disp('backGrond correction')
% % dirAll=superBatchConvert();
% dirAll='总文件夹名';
% nameList=dir(dirAll);
% cd(dirAll);
% for i=1:(length(nameList)-2)
%     disp(strcat('AllJob',num2str(i)));
%     disp(nameList(i+2).name);
%     dirFile=strcat(dirAll,'\',nameList(i+2).name);
%     bioTree=load(strcat(dirFile,'\bioTreeResult\allTree\bioTreeStack'));
%     bioTree=bioTree.bioTree;
%     %     bioTree=batchTreeTrackingJob(dirFile);
%     %     [bioTree,allDataOri]=runaLL(bioTree);
%     %     save(strcat(dirFile,'\allDataOri'),'allDataOri');
%     allDataOri=load(strcat(dirFile,'\allDataOri'));
%     allDataOri=allDataOri.allDataOri;
% %     startFrame=1;
% %     endFrame=size(bioTree,2);
% %     motilityDemoMaker(allDataOri,startFrame,stepFrame,endFrame,dirFile);
%     [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri);
%     disp('rockFrame is below: ');disp(rockFrame);
%     correctionList=autoGetCorrectionList(bioTree,allDataOri);
%     allData=dataCorrection(allDataOri,correctionList,[],[],dirFile,rockFrame);
%     %     load(strcat(dirFile,'\allData'))
%     dirResultSave=strcat(dirFile,'\allResult');
%     mkdir(dirResultSave);
%     allDataExplorer(allData,dirResultSave);
% end
% end
% 
% 
% 
% % 直接根据allData做图
% function getBioTreeAndMakeDemo(stepFrame)
% % disp('backGrond correction')
% % dirAll=superBatchConvert();
% dirAll='总文件夹名';
% nameList=dir(dirAll);
% cd(dirAll);
% for i=1:(length(nameList)-2)
%     disp(strcat('AllJob',num2str(i)));
%     disp(nameList(i+2).name);
%     dirFile=strcat(dirAll,'\',nameList(i+2).name);
% %     bioTree=load(strcat(dirFile,'\bioTreeResult\allTree\bioTreeStack'));
% %     bioTree=bioTree.bioTree;
% %     bioTree=batchTreeTrackingJob(dirFile);
% %     [bioTree,allDataOri]=runaLL(bioTree);
% %     save(strcat(dirFile,'\allDataOri'),'allDataOri');
% %     allDataOri=load(strcat(dirFile,'\allDataOri'));
% %     allDataOri=allDataOri.allDataOri;
% %     startFrame=1;
% %     endFrame=size(bioTree,2);
% %     motilityDemoMaker(allDataOri,startFrame,stepFrame,endFrame,dirFile);
% %     [rockFrame,allDataOri]=findRockFrameAndDelete(allDataOri);
% %     disp('rockFrame is below: ');disp(rockFrame);
% %     correctionList=autoGetCorrectionList(bioTree,allDataOri);
% %     allData=dataCorrection(allDataOri,correctionList,[],[],dirFile,rockFrame);
%     allData=load(strcat(dirFile,'\allData'));
%     allData=allData.allData;
%     dirResultSave=strcat(dirFile,'\allResult');
%     mkdir(dirResultSave);
%     allDataExplorer(allData,dirResultSave);
% end
% end