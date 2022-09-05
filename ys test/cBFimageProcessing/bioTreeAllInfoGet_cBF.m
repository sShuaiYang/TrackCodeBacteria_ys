function [ bioTree ] = bioTreeAllInfoGet_cBF( bioTree,dirFile)
%BIOTREEALLINFOGET Summary of this function goes here
%   Detailed explanation goes here
load([dirFile,'\Tracking\frameInfo.mat']);
bioTreeInitialTime=frameInfo(1,1:6);
frameInfoNew(1:2:numel(bioTree)-1,:)=frameInfo(1:numel(bioTree)/2,:);
frameInfoNew(2:2:numel(bioTree),:)=frameInfo(1:numel(bioTree)/2,:);
% frameInfoNew=frameInfo;
% bestPositionAccumulation=bioTree{1}.bestPositionAccumulation;
for i=1:numel(bioTree)
    bioTreeTimer(i)=etime(frameInfoNew(i,1:6),bioTreeInitialTime)/60;   %min
end
bioTree{1}.bioTreeInitialTime=bioTreeInitialTime;
bioTree{1}.bioTreeTimer=bioTreeTimer;

% load([dirFile,'\Tracking\frameInfo.mat'])

bioTreeFrame=[];
nameList=dir([dirFile,'\Tracking']);
for i=1:size(frameInfo,1)
    bioTreeTimer=etime(frameInfo(i,1:6),bioTree{1}.bioTreeInitialTime)/60;
    diffTime=abs(bioTree{1}.bioTreeTimer-bioTreeTimer);
    [~,index]=find(diffTime==min(diffTime));
    bioTreeFrame(i)=index(1);
end
[~,checkFrame]=find(diff(bioTreeFrame)==0);
bioTreeFrame(checkFrame+1)=0;
for iPic=1:numel(bioTreeFrame)
    if bioTreeFrame(iPic)==0
        continue
    end
    %         image=load([dirFluoFile,'\image',fluoChannel{iChannel},num2str(iPic,'%05.f'),'.mat']);
%     tempimage=load([dirFile,'\Tracking,'\image','Tacking',num2str(iPic,'%05.f'),'.mat']);
%     image=tempimage.imageTracking;
    
    for iframe=1:bioTreeFrame(iPic)
        for iRoot=1:numel(bioTree{iframe}.root)
            traceInfo=bioTree{iframe}.root{iRoot}.traceInfo;
            bioTree{iframe}.root{iRoot}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
            if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
                tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
            end
        end
        for iNode=1:numel(bioTree{iframe}.node)
            for iOut=1:numel(bioTree{iframe}.node{iNode}.Out)
                traceInfo=bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo;
                bioTree{iframe}.node{iNode}.Out{iOut}.linkTimer=bioTree{1}.bioTreeTimer(iframe:iframe+numel(traceInfo.pixelIdxList)-1)';
                if numel(traceInfo.pixelIdxList)+iframe>bioTreeFrame(iPic)
                    tracePixel=traceInfo.pixelIdxList{bioTreeFrame(iPic)-iframe+1};
                end
            end
        end
    end
end

end