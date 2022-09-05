function timeLapseMultiFieldAbsCentroidGetAndIntCorrection(dirFile)
disp('Cell abs centroid Get and fluo intensity correction')
fieldList = dir([dirFile,filesep,'field*']);
for iField= 1:(length(fieldList))
    %     t0 = clock;
    if ~strcmp(fieldList(iField).name(1:5),'field')
        continue
    end
    disp(fieldList(iField).name);

    dirField = strcat(dirFile,'\',fieldList(iField).name);
    load([dirField,'\bioInfo.mat']);

    %%细菌绝对位置的获得 unit um
    %     [bioInfo] = absCentroidGetInMultiFields(bioInfo,dirFile);
    %荧光光强校正
    [bioInfo] = cellFluoIntensityCorrection(bioInfo,dirFile);
    save (strcat(dirField,'\bioInfo.mat'),'bioInfo','-v7.3');
    clear bioInfo
    %     t1 = clock;
    %     disp(['共耗时:',num2str(etime(t1,t0)),'秒']);
end
end