function [processedImages,mergedImages] = phCImProcessing_supperSegger(phaseImages)
% 利用supperSegger对P.aeruginosa的相差图像（phase contrast Images）进行分割
% Shuai Yang 2020.09.24
% convert images type to 8bit
convertedImages = uint8 (zeros(size(phaseImages)));
parfor iFrame = 1:size(phaseImages,3)
    I0 = phaseImages(:,:, iFrame);
    I0 = rescale(double(I0))*255;
    I0 = uint8(I0);
    convertedImages (:,:,iFrame) = I0;
end
phaseImages = convertedImages;

im_sz = [size(phaseImages,1),size(phaseImages,2),...
    size(phaseImages,3)];% image or imageStack size

processedImages = false(im_sz);
mergedImages = zeros(im_sz(1),im_sz(2),3,im_sz(3));
mergedImages = uint8(mergedImages);

CONST = loadConstants ('100XPa',0, 0);
% not verbose state
CONST.parallel.verbose = 0;
parfor iFrame = 1:size(phaseImages,3)%非必需用parfor
% for iFrame = 1:size(phaseImages,3)
    dataname = 'PhCIm';
    phaseIm = phaseImages(:,:,iFrame);
    segFile = CONST.seg.segFun( phaseIm, CONST, '', dataname, []);
    cell_mask = logical(segFile.mask_cell);
    cell_mask = imfill( cell_mask,'holes');
    cell_mask = imopen( cell_mask,strel('disk',2));
    processedImages(:,:,iFrame) = cell_mask;
    mergedImages(:,:,:,iFrame) = showSegDataMergePhase( segFile);
    close all
end
%% 
% delete wrong size region based on cell morphology
[processedImages] = removeWrongMaskRegions_cellMorphology(processedImages);
%%
% if numel(im_sz) == 2
%    figure;
%    imshow(mergedImages) 
% end
% imshowlabel_ys(processedImages(:,:,end));
end
%%
function [mergedImage] = showSegDataMergePhase(segFile)
%将segData 与phase image Merge
% Shuai Yang 2020.09.24
cell_mask = logical(segFile.mask_cell);
cell_mask = imfill( cell_mask,'holes');
outline = bwmorph(cell_mask,'remove');
mergedImage = imoverlay(segFile.phase,outline,[1 0 0]);
end