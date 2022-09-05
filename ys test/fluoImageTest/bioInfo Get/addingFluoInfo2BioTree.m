function [ bioTree ] = addingFluoInfo2BioTree( bioTree,dirField)
%适用于bioTree由tracking生成，但拍摄的荧光图像与tracking不同步
%核心函数 bioTreeAllInfoGet
%shuai yang 2020.05.21
load([dirField,'\Tracking\frameInfo.mat'])
frameInfo_tracking = frameInfo;
trackingNum = numel(bioTree);
bioTreeInitialTime = frameInfo_tracking(1,1:6);
bioTreeTimer = zeros(1,numel(bioTree));
for i = 1:trackingNum
    bioTreeTimer(i)=etime(frameInfo_tracking(i,1:6),bioTreeInitialTime)/60;   %min
end
bioTree{1}.bioTreeInitialTime = bioTreeInitialTime;
bioTree{1}.bioTreeTimer = bioTreeTimer;

fluoChannels = {'sfGFP','mScarletI','Venus','PVD','CyOFP','TDsmURFP'};
fluoIntNames = {'intsfGFP','intmScarletI','intVenus','intPVD','intCyOFP','intTDsmURFP'};

bioTree = bioTreeAllChannel(bioTree,fluoIntNames);
for iChannel = 1:numel(fluoChannels)
    dirImage=[dirField,'\',fluoChannels{iChannel}];
    try
        load([dirImage,'\frameInfo.mat'])
    catch ME
        continue
    end
    
%     imageList = dir(dirImage);
    imageList = dir([dirImage,filesep,'*.tif']);
    frameInfo_fluo = frameInfo;
    % 寻找tracking的时间窗口内，所拍摄的荧光图像的数目
    fluo2TrackingIdx = zeros (1,size (frameInfo_fluo,1));
    for iFluo = 1 :size (frameInfo_fluo,1)
        time_lag = abs(etime (frameInfo_fluo(iFluo,1:6),frameInfo_tracking(:,1:6))/60);
        %min 与frameInfo_tracking拍摄的时间间隔绝对值
        trckingIdx = find(time_lag==min(time_lag));
        fluo2TrackingIdx(iFluo) = trckingIdx(end);
        %use trckingIdx(end),for the case of trckingIdx 不唯一 ;case:1,1,2,2...
    end
    fluo2TrackingIdx = fluo2TrackingIdx(fluo2TrackingIdx <= trackingNum);
    [~,checkFrame] = find(diff(fluo2TrackingIdx)==0);
    fluo2TrackingIdx(checkFrame+1)=0;% case for:两张荧光同时对应一张tracking
    bioTreeFrame = fluo2TrackingIdx;
    
    % 寻找所有root和node中traceInfo中包含fluo2TrackingIdx的fluo图像的frame 
    for iPic = 1:numel(bioTreeFrame)
        if bioTreeFrame(iPic)==0
            continue
        end
        imFluo = import_tiff_stack( strcat(dirImage,'\',imageList(iPic).name) );
        
        for iFrame = 1:bioTreeFrame(iPic)
            for iRoot = 1:numel(bioTree{iFrame}.root)
                traceInfo = bioTree{iFrame}.root{iRoot}.traceInfo;
                bioTree{iFrame}.root{iRoot}.linkTimer = bioTree{1}.bioTreeTimer(iFrame:iFrame+numel(traceInfo.pixelIdxList)-1)';
                if numel(traceInfo.pixelIdxList)+iFrame>bioTreeFrame(iPic)
                    tracePixel = traceInfo.pixelIdxList{bioTreeFrame(iPic)-iFrame+1};
                    %找到与fluoImage对应的那一个frame，进行赋值，没有的仍然是NaN
                    switch iChannel
                        case 1
                            bioTree{iFrame}.root{iRoot}.intsfGFP(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                        case 2
                            bioTree{iFrame}.root{iRoot}.intmScarletI(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                        case 3
                            bioTree{iFrame}.root{iRoot}.intVenus(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                        case 4
                            bioTree{iFrame}.root{iRoot}.intPVD(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                        case 5
                            bioTree{iFrame}.root{iRoot}.intCyOFP(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                        case 6
                            bioTree{iFrame}.root{iRoot}.intTDsmURFP(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                    end
                end
            end
            for iNode=1:numel(bioTree{iFrame}.node)
                for iOut=1:numel(bioTree{iFrame}.node{iNode}.Out)
                    traceInfo=bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo;
                    bioTree{iFrame}.node{iNode}.Out{iOut}.linkTimer=bioTree{1}.bioTreeTimer(iFrame:iFrame+numel(traceInfo.pixelIdxList)-1)';
                    if numel(traceInfo.pixelIdxList)+iFrame>bioTreeFrame(iPic)
                        tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iFrame+1};
                        switch iChannel
                            case 1
                                bioTree{iFrame}.node{iNode}.Out{iOut}.intsfGFP(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                            case 2
                                bioTree{iFrame}.node{iNode}.Out{iOut}.intmScarletI(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                            case 3
                                bioTree{iFrame}.node{iNode}.Out{iOut}.intVenus(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                            case 4
                                bioTree{iFrame}.node{iNode}.Out{iOut}.intPVD(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                            case 5
                                bioTree{iFrame}.node{iNode}.Out{iOut}.intCyOFP(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                            case 6
                                bioTree{iFrame}.node{iNode}.Out{iOut}.intTDsmURFP(bioTreeFrame(iPic)-iFrame+1) = mean(imFluo(tracePixel));
                        end
                    end
                end
            end
        end
    end
end


end
function bioTree=bioTreeAllChannel(bioTree,fluoIntNames)
for iFrame=1:numel(bioTree)
    for iRoot=1:numel(bioTree{iFrame}.root)
        for iChannel = 1:numel(fluoIntNames)
            bioTree{iFrame}.root{iRoot}.(fluoIntNames{iChannel}) =...
                ones(numel(bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList),1)*nan;
        end
    end
    for iNode=1:numel(bioTree{iFrame}.node)
        for iOut=1:numel(bioTree{iFrame}.node{iNode}.Out)
            for iChannel = 1:numel(fluoIntNames)
                bioTree{iFrame}.node{iNode}.Out{iOut}.(fluoIntNames{iChannel}) =...
                    ones(numel(bioTree{iFrame}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList),1)*nan;
            end
        end
    end
end
end