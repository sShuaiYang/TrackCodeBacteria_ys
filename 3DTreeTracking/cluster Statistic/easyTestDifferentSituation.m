function easyTestDifferentSituation(bioTree,bacteriaFrameInfo,dirFile)
sizeRatio=[0.5,1.5,2.5,4];
for i=1:numel(sizeRatio)
    clusterTree=accumulateNetwork(bacteriaFrameInfo,bioTree,200,sizeRatio(i));
    dirResult=strcat(dirFile,'\',num2str(sizeRatio(i)));
    mkdir(dirResult)
    cd(dirResult)
    save('clusterTree.mat','clusterTree','bacteriaFrameInfo');
    testClusterTree(clusterTree,dirResult,200);
end
end