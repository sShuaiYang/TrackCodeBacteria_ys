function [bioInfo] = bioInfoGetByBioTree(bioTree,fluo2TrackingIdx,dirField)
%通过bioTree生成bioInfo 结构体
% shuai Yang 2020.05.21

bioInfo = struct('fieldIdx',[],'frameIdx',[],'trackingIdx',[],'intsfGFP',[ ],'intmScarletI',[ ],...
    'intVenus',[ ],'intPVD',[ ],'intCyOFP',[ ],'intTDsmURFP',[ ], ...
    'BG_sfGFP',[ ],'BG_mScarletI',[ ],'BG_Venus',[ ],'BG_PVD',[ ],'BG_CyOFP',[ ], ...
    'BG_TDsmURFP',[ ],'Centroid',[ ],'MajorAxisLength',[ ], 'MinorAxisLength',[ ], ...
    'FilledArea',[ ],'PixelIdxList',{ },'PixelList',{ },'gR',[],'cellLineage',{});

fluoIntNames = {'intsfGFP','intmScarletI','intVenus','intPVD','intCyOFP','intTDsmURFP'};
cellLineage = struct('birthField',[],'birthFrame',[],'birthRoot',[],'label',[]);

bioTreeFrame = fluo2TrackingIdx;
for  iPic = 1:numel(bioTreeFrame)
    k = 1;
    bioInfo(iPic).fieldIdx = str2double(dirField(end-3:end));
    bioInfo(iPic).frameIdx = iPic;
    bioInfo(iPic).trackingIdx = bioTreeFrame(iPic);
    for iFrame = 1:bioTreeFrame(iPic)
        for iRoot = 1:numel(bioTree{iFrame}.root)
            traceInfo = bioTree{iFrame}.root{iRoot}.traceInfo;
            if numel(traceInfo.pixelIdxList)+iFrame>bioTreeFrame(iPic)
                for  iChannel = 1:numel(fluoIntNames)
                    bioInfo(iPic).(fluoIntNames{iChannel})(k,:)= ...
                        bioTree{iFrame}.root{iRoot}.(fluoIntNames{iChannel})(bioTreeFrame(iPic)-iFrame+1);
                end
                bioInfo(iPic).Centroid(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.Centroid;
                bioInfo(iPic).MajorAxisLength(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.MajorAxisLength;
                bioInfo(iPic).MinorAxisLength(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.MinorAxisLength;
                bioInfo(iPic).FilledArea(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.FilledArea;
                bioInfo(iPic).PixelIdxList{k} = traceInfo.pixelIdxList{bioTreeFrame(iPic)-iFrame+1};
                bioInfo(iPic).PixelList{k} = pixelIdxList2XY(bioInfo(iPic).PixelIdxList{k},[2048,2048]);
                bioInfo(iPic).gR(k,:) = bioTree{iFrame}.root{iRoot}.gR;
                if isfield(bioTree{iFrame}.root{iRoot},'cellLineage')
                    bioInfo(iPic).cellLineage(k) = bioTree{iFrame}.root{iRoot}.cellLineage;
                else
                    bioInfo(iPic).cellLineage(k) = cellLineage;
                end
                k = k+1;
            end
            
        end
        
        for iNode = 1:numel(bioTree{iFrame}.node)
            for iOut = 1:numel(bioTree{iFrame}.node{iNode}.Out)
                traceInfo = bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo;
                if numel(traceInfo.pixelIdxList)+iFrame>bioTreeFrame(iPic)
                    for  iChannel = 1:numel(fluoIntNames)
                        bioInfo(iPic).(fluoIntNames{iChannel})(k,:)= ...
                            bioTree{iFrame}.node{iNode}.Out{iOut}.(fluoIntNames{iChannel})(bioTreeFrame(iPic)-iFrame+1);
                    end
                    bioInfo(iPic).Centroid(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.Centroid;
                    bioInfo(iPic).MajorAxisLength(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.MajorAxisLength;
                    bioInfo(iPic).MinorAxisLength(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.MinorAxisLength;
                    bioInfo(iPic).FilledArea(k,:) = traceInfo.measurment{bioTreeFrame(iPic)-iFrame+1}.FilledArea;
                    bioInfo(iPic).PixelIdxList{k} = traceInfo.pixelIdxList{bioTreeFrame(iPic)-iFrame+1};
                    bioInfo(iPic).PixelList{k} = pixelIdxList2XY(bioInfo(iPic).PixelIdxList{k},[2048,2048]);
                    bioInfo(iPic).gR(k,:) = bioTree{iFrame}.node{iNode}.Out{iOut}.gR;
                    if isfield(bioTree{iFrame}.node{iNode}.Out{iOut},'cellLineage')
                        bioInfo(iPic).cellLineage(k) = bioTree{iFrame}.node{iNode}.Out{iOut}.cellLineage;
                    else
                        bioInfo(iPic).cellLineage(k) = cellLineage;
                    end
                    k = k+1;
                end
            end
        end
    end
    
    % 如果等于NaN 赋值为空[]
    if k>1% 此fluo field至少有一个菌
        for  iChannel = 1:numel(fluoIntNames)
            if isnan(bioInfo(iPic).(fluoIntNames{iChannel})(1))
                bioInfo(iPic).(fluoIntNames{iChannel}) = [];
            end
        end
    end
    
end


end
% 此函数 现在可以用ind2sub替换
function PixelList = pixelIdxList2XY(PixelIdxList,imSize)
image = false (imSize);
image(PixelIdxList) = true;
stats = regionprops(image,'PixelIdxList','PixelList');
PixelList = stats.PixelList;
end
