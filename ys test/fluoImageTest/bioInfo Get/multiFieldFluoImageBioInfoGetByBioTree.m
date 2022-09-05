function multiFieldFluoImageBioInfoGetByBioTree(dirFile)
disp('Each field bioInfo get by bioTree');
fieldList =  dir([dirFile,filesep,'field*']);
for iField = 1:length(fieldList)
    t0=clock;
    if ~strcmp(fieldList(iField).name(1:5),'field')
        continue
    end
    disp(fieldList(iField).name);
    disp(compose("Start time %d-%d-%d %.2d:%.2d:%02.0f.",t0));
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    load([dirField,'\bioTreeResult\finalTree\bioTree.mat']);
    
    [fluo2TrackingIdx] = findTrackingImagesForUnsynFluo(dirField);
    bioTree{1}.fluo2TrackingIdx = fluo2TrackingIdx;
    disp('adding Fluo Info to BioTree');
    [bioTree] = addingFluoInfo2BioTree( bioTree,dirField);
    disp('adding Growth Rate to BioTree');
    [bioTree] = addingGrowthRate2BioTree(bioTree);
    disp('adding Cell Lineage to BioTree');
    [bioTree] = addingCellLineage2BioTree( bioTree,dirField);
    disp('bioInfo Get By BioTree');
    [bioInfo] = bioInfoGetByBioTree(bioTree,fluo2TrackingIdx,dirField);

    disp('Save bioInfo and BioTree');
    save (strcat(dirField,'\bioTree.mat'),'bioTree','-v7.3');
    save (strcat(dirField,'\bioInfo.mat'),'bioInfo','-v7.3');

    clear bioTree bioInfo
    t1 = clock;
    disp(compose("End time %d-%d-%d %.2d:%.2d:%02.0f.",t1))
    disp(['共耗时:',num2str(etime(t1,t0)),'秒']);
end
end