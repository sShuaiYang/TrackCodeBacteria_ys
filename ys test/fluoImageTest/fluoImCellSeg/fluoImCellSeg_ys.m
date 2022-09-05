 function [processedImages] = fluoImCellSeg_ys(fluoImages)
%% 利用edge图像进行边缘检测
% watershed transform进行segementaion
% shuai Yang 2020.10.11
% fluoImSeg_edge_Test.mlx;segCellTest_Watershedtransform.mlx;
% segCellTest_splitKinkyCells.mlx;segCellTest_splitFatCells.mlx
% default 100x P.aeruginosa paramters; uint pixel
% objective  = '100x';default
p.minAreaThreshold = 50;
p.maxAreaThreshold = 600;
p.minCellLength = 10;%[20 30]
p.maxCellLength = 64;%~65
p.maxCellWidth = 15;
p.minSolidity = 0.75;%[0.75 0.85]
p.maxSegLineNum = 4;
p.minCellThickness = sqrt(2);
p.H = 1.6;%imextendedmin Extended-minima transform-minima transform, nonnegative scalar
p.maxH = 2;
p.conservativeH = 2.2;% used for split long cell；case for overseg

%%
% get mask
% &&& not suitable for  weak fluorescence &&&%
% &&& not suitable for  strong fluo loci in cells &&& %
processedImages = false(size(fluoImages));

parfor iFrame = 1:size(fluoImages,3)
    fluoIm = fluoImages(:,:, iFrame);
    [fluoImMask] = fluoImCellMask_ys(fluoIm);
    %     [fluoImMask] = fluoImCellMask_edge_w60x(fluoIm);%water 60x
    processedImages(:,:,iFrame) = fluoImMask;
end


%%
% % mask segmentation
% parfor iFrame = 1:size(processedImages,3)% 可以不用par
% % for  iFrame = 1:size(processedImages,3)
%     bw = processedImages(:,:, iFrame);
%     segbw = fluoCellMaskSeg(bw,p);% one image
%     processedImages(:,:, iFrame) = segbw;
% end
%%
% delete wrong size region based on cell morphology
% [processedImages] = removeWrongMaskRegions_cellMorphology(processedImages);
% imshowlabel_ys(processedImages(:,:,end));
end

