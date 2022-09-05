function BatchParticleSimulation()
dirFile=uigetdir();
cd(dirFile);
% pesistanceLength=10;
pesistanceLength=logspace(1,2,8);  %[0,180]  0 means line, 180 means stochastic move
for i=1:numel(pesistanceLength)
    dirResultFile=strcat(dirFile,'\',num2str(pesistanceLength(i)));
    mkdir(dirResultFile);
    cd(dirResultFile)
%     load('bestBioTree\bioTree_7.mat');
%     [bioTree,branchList,~,~]=myBiograph_new2(bioTree);
%     bacteriaFrameInfo=getEachBacteriaInFrame(bioTree);
%     save('bacteriaFrameInfo','bacteriaFrameInfo');
    stochasticMovingSimulation(30,[2160,2560],4500,pesistanceLength(i),dirResultFile);
end
end
function bacteriaFrameInfo=getEachBacteriaInFrame(bioTree)
for iframe=1:size(bioTree,2)
    bacteriaFrameInfo{iframe}.centroidInfo=[];
    bacteriaFrameInfo{iframe}.lengthInfo=[];
    bacteriaFrameInfo{iframe}.bacteriaInfo=[];
end
for iframe=1:size(bioTree,2)
    disp(iframe)
    if ~isempty(bioTree{iframe}.root)
        if ~isempty(bioTree{iframe}.root)
            for iRoot=1:size(bioTree{iframe}.root,2)
                traceInfo=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList;
                if isempty(traceInfo)
                    continue
                end
                for iTrace=1:size(traceInfo,2)
                    yCo=ceil((bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{iTrace}-1)/2160);
                    xCo=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{iTrace}-2160*(yCo-1);
                    bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;[yCo,xCo]];
                    bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo=[bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo;[iframe,iRoot,0,iTrace,bioTree{iframe}.root{iRoot}.branchIndex]];
                end
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2);
                traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList;
                for iTrace=1:size(traceInfo,2)
                    yCo=ceil((bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}-1)/2160);
                    xCo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}-2160*(yCo-1);
                    bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;[yCo,xCo]];
                    bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo=[bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo;[iframe,iNode,iOut,iTrace,bioTree{iframe}.node{iNode}.branchIndex]];
                end
            end
        end
    end
end
end