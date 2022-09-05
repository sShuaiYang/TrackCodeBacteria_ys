function [bioInfo] = absCentroidGetInMultiFields(bioInfo,dirFile)
%  multifield 中细菌绝对位置的获得 unit um
%需要mip文件 shuai Yang 2020.05.19
disp('Multi field cells abs centroid get ')
load([dirFile,'\mip.mat']);
scale = mip.Calib.scale;
% scale=0.0650;
xyPosition = mip.fieldTag.xyPosition*0.1;% unit um
% micromanager code V4.0 code %unit 0.1 um
% micromanager code V3.0 code %unit 1 um
fieldPosition = xyPosition-xyPosition(1,:);
if isempty(bioInfo)
    return
end
for iInfo = 1:numel(bioInfo)
    if  isempty(bioInfo(iInfo).Centroid)
        continue
    end
    fieldTag =  bioInfo(iInfo).fieldIdx;
    bioInfo(iInfo).absCentroid = bioInfo(iInfo).Centroid *scale + fieldPosition(fieldTag,:);
end
% save(strcat(dirFile,'\bioInfo.mat'),'bioInfo','-v7.3');
end