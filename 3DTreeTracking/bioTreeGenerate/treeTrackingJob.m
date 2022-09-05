function [demoMovie,bioTree]=treeTrackingJob(imageStack,maskImages,firstFrame,endFrame,frameShift)
% clc;
disp('1.generate CC Array, using multicore computing');
tic;treeBranch=makeBranch(maskImages,firstFrame,endFrame);toc; %generate CC Array, using multicore computing
disp('2.generate the bioTree');
xSize=treeBranch(1).ImageSize(1);
ySize=treeBranch(1).ImageSize(2);
tic;bioTree=makeBioTree(treeBranch,firstFrame,endFrame,frameShift,xSize,ySize);toc; %generate the bioTree, and one by one connect.
% disp('3.measure the object in the tree');
% tic;bioTree=bioTreeMeasure(bioTree,frameShift,xSize,ySize);toc;%measure the object in the Tree
% disp('4.refind the leaf and root in the Tree');
% tic;bioTree=leafRootRefineBioTree(bioTree,frameShift);toc; % further connected the leaf in the ith frame and the root in the i+1 frame, in this function you need set a distence threshold for the link
% disp('5.refind the Node');
% tic;bioTree=bioTreeSizeReduction(bioTree,frameShift);
%     bioTree=nodeEdgeReduction(bioTree,frameShift);
% toc;
% bioTree=nodeRefineBioTree(bioTree,frameShift);
% bioTree=nodeConfusionRefine(bioTree,frameShift); % further Refind the Node in bioTree

% disp('5.Node measurement');
% tic;bioTree=nodeMeasurement(bioTree);tic;
% disp('6.basic counting');
% tic;result=basicSearching(bioTree);tic;
% plotResult(result);
demoMovie=1;
% disp('7.make demo Movie');
% tic;demoMovie=treeDemoMaker(imageStack,bioTree);toc;
end