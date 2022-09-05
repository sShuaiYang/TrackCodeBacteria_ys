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
                    if iframe+iTrace-1>numel(bioTree)
                        continue
                    end
                    bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;bioTree{iframe}.root{iRoot}.traceInfo.measurment{iTrace}(1).Centroid];
% yCo=ceil((bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{iTrace}-1)/2160);
% xCo=bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList{iTrace}-2160*yCo;
% bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;[yCo,xCo]];
                    bacteriaFrameInfo{iframe+iTrace-1}.lengthInfo=[bacteriaFrameInfo{iframe+iTrace-1}.lengthInfo;bioTree{iframe}.root{iRoot}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
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
                    if iframe+iTrace-1>numel(bioTree)
                        continue
                    end
                    bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).Centroid];
% yCo=ceil((bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}-1)/2160);
% xCo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList{iTrace}-2160*yCo;
% bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo=[bacteriaFrameInfo{iframe+iTrace-1}.centroidInfo;[yCo,xCo]];
                    bacteriaFrameInfo{iframe+iTrace-1}.lengthInfo=[bacteriaFrameInfo{iframe+iTrace-1}.lengthInfo;bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment{iTrace}(1).MajorAxisLength];
                    bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo=[bacteriaFrameInfo{iframe+iTrace-1}.bacteriaInfo;[iframe,iNode,iOut,iTrace,bioTree{iframe}.node{iNode}.branchIndex]];
                end
            end
        end
    end
end
end
% function bacteriaFrameInfo=getEachBacteriaInFrame(bioTree)
% for iframe=size(bioTree,2)
%     disp(iframe)
%     bacteriaFrameInfo{iframe}.centroidInfo=[];
%     bacteriaFrameInfo{iframe}.lengthInfo=[];
%     bacteriaFrameInfo{iframe}.bacteriaInfo=[];
%     for smallFrame=1:iframe
%         if ~isempty(bioTree{smallFrame}.root)
%             if ~isempty(bioTree{smallFrame}.root)
%                 for iRoot=1:size(bioTree{smallFrame}.root,2)
%                     traceInfo=bioTree{smallFrame}.root{iRoot}.traceInfo.pixelIdxList;
%                     if isempty(traceInfo)
%                         continue
%                     end
%                     if size(traceInfo,2)>=iframe+1-smallFrame
%                         branchIndex=bioTree{smallFrame}.root{iRoot}.branchIndex;
%                         bacteriaFrameInfo{iframe}.centroidInfo=[bacteriaFrameInfo{iframe}.centroidInfo;bioTree{smallFrame}.root{iRoot}.traceInfo.measurment{iframe+1-smallFrame}(1).Centroid];
%                         bacteriaFrameInfo{iframe}.lengthInfo=[bacteriaFrameInfo{iframe}.lengthInfo;bioTree{smallFrame}.root{iRoot}.traceInfo.measurment{iframe+1-smallFrame}(1).MajorAxisLength];
%                         bacteriaFrameInfo{iframe}.bacteriaInfo=[bacteriaFrameInfo{iframe}.bacteriaInfo;[smallFrame,iRoot,0,iframe+1-smallFrame,bioTree{smallFrame}.root{iRoot}.branchIndex]];
%                     end
%                 end
%             end
%         end
%         if ~isempty(bioTree{smallFrame}.node)
%             for iNode=1:size(bioTree{smallFrame}.node,2)
%                 for iOut=1:size(bioTree{smallFrame}.node{iNode}.Out,2);
%                     traceInfo=bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList;
%                     if size(traceInfo,2)>=iframe+1-smallFrame
%                         branchIndex=bioTree{smallFrame}.node{iNode}.branchIndex;
%                         bacteriaFrameInfo{iframe}.centroidInfo=[bacteriaFrameInfo{iframe}.centroidInfo;bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.measurment{iframe+1-smallFrame}(1).Centroid];
%                         bacteriaFrameInfo{iframe}.lengthInfo=[bacteriaFrameInfo{iframe}.lengthInfo;bioTree{smallFrame}.node{iNode}.Out{iOut}.traceInfo.measurment{iframe+1-smallFrame}(1).MajorAxisLength];
%                         bacteriaFrameInfo{iframe}.bacteriaInfo=[bacteriaFrameInfo{iframe}.bacteriaInfo;[smallFrame,iNode,iOut,iframe+1-smallFrame,bioTree{smallFrame}.node{iNode}.branchIndex]];
%                     end
%                 end
%             end
%         end
%     end
% end
% end