function [processedImages] = removeWrongMaskRegions_cellMorphology(processedImages)
%Shuai Yang 2020.12.11
% delete wrong size region based on cell morphology
scale = 0.065;%100x pixel to um
areaThreshold = 80;
% BW = bwareafilt(BW,[minAreaThreshold maxAreaThreshold]);直接筛选 Shuai Yang
% 2021.07.15
% BW2 = bwpropfilt(BW,'EulerNumber',[1 1]);
parfor  iFrame = 1:size(processedImages,3)
    bw = processedImages(:,:, iFrame);
    bw = bwareaopen(bw,areaThreshold); %remove the objecvies which area is smaller than threshold
    L = bwlabel(bw);
    stats = regionprops(bw,'MajorAxisLength','MinorAxisLength');
    cellWid = [stats.MinorAxisLength]*scale;
    cellLen = [stats.MajorAxisLength]*scale;
    fwrong = find((cellLen./ cellWid) < 1.2 | cellWid/2 > 0.7 | cellLen < 0.5);
    % L>1.2*W & (W/2)<0.7 & L >0.5;right region
    if ~isempty(fwrong)
        for iRegion = 1:numel(fwrong)
            L(L==fwrong(iRegion))= 0;
        end
    end
    processedImages(:,:, iFrame) = logical(L);
end
end

