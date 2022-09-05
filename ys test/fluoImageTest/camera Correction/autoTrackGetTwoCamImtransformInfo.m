function [transformInfo] = autoTrackGetTwoCamImtransformInfo(Im_cam1,Im_cam2,dirFile)
% 通过追踪 两张图像有overlap 生成biotree  来获取校正的点
% Im_cam1 相机1拍摄；Im_cam2 相机2拍摄
% Shuai Yang 2020.05.27
% camera 1 拍摄的为fixed channels，camera2 拍摄的为moving channels

%y轴镜像对称，2;x轴镜像对称1；
Im_cam2 = flip(Im_cam2,2);%图像镜像对称 Micromanger2; MATLAB 1

disp('Auto get camera transformInfo by two Cam images')
allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道

%初始化
for iChannel = 1:numel(allChannels)
    transformInfo.(allChannels{iChannel}) = NaN;
end

%'PhC','Fluo','Binary'
[fixedImageMask,ImType_cam1] = imageTypeGetAndSeg(Im_cam1); 
disp(['Im_Cam1 is ',ImType_cam1, 'Image'])
[movingImageMask,ImType_cam2] = imageTypeGetAndSeg(Im_cam2);
disp(['Im_Cam2 is ',ImType_cam2, 'Image'])

[selectedMovingPoints,selectedFixedPoints]=...
    autoTrackGetCalibrationPoints_ys(fixedImageMask,movingImageMask);
tform = getTransformInfoAndCheckResult(selectedMovingPoints,...
    selectedFixedPoints,Im_cam2,Im_cam1);

for iChannel = 1:numel(allChannels)
    if ismember(allChannels{iChannel},channels_cam1)
        transformInfo.(allChannels{iChannel}) = 0;
    end
    if ismember(allChannels{iChannel},channels_cam2)
        transformInfo.(allChannels{iChannel}) = tform;
    end
end
save(strcat(dirFile,'\transformInfo.mat'),'transformInfo');
end

%%
function [selectedMovingPoints,selectedFixedPoints] = autoTrackGetCalibrationPoints_ys(fixedImage,movingImage)

centerMarkFixedImage = false(size(fixedImage));
centerMarkMovingImage = false(size(movingImage));

fp = [];
mp = [];

% 通过两张图像tracking 获得细菌的质心来找到对应点
trackingImages = cat(3,fixedImage,movingImage);

[ ~ , bioTree] = bacteriaTracking( trackingImages(:,:,1) , trackingImages(:,:,2) );
bioTree = bioTreeMeasure(bioTree,0,size(trackingImages,1),size(trackingImages,2));

n = 1;
for iroot=1:size(bioTree{1}.root,2)
    is2Node=bioTree{1}.root{iroot}.is2Node;
    leafInfo=bioTree{1}.root{iroot}.leafInfo;
    if is2Node == 0 && leafInfo(1)==2
        pixelIdxList1 = bioTree{1}.root{iroot}.traceInfo.pixelIdxList{1};
        pixelIdxList2 = bioTree{1}.root{iroot}.traceInfo.pixelIdxList{2};
        if numel(pixelIdxList1)>1500 || numel(pixelIdxList1)<100
            continue
        end
        if abs(numel(pixelIdxList1)-numel(pixelIdxList2))>20
            continue
        end
        centerMarkFixedImage(pixelIdxList1) = true;
        centerMarkMovingImage(pixelIdxList2) = true;
        if size(bioTree{1}.root{iroot}.traceInfo.measurment{1},1)==1 && ...
                size(bioTree{1}.root{iroot}.traceInfo.measurment{2},1)==1
            Centroid1 = bioTree{1}.root{iroot}.traceInfo.measurment{1}.Centroid;
            Centroid2 = bioTree{1}.root{iroot}.traceInfo.measurment{2}.Centroid;
            if pdist2(Centroid1,Centroid2)<= 10
                fp(n,1:2) = Centroid1;
                mp(n,1:2) = Centroid2;
                n = n+1;
            end
        end
    end
    
end

% selectedMovingPoints mp;selectedFixedPoints fp；
selectedFixedPoints = fp;
selectedMovingPoints = mp;

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'projective');
Rfixed = imref2d(size(centerMarkMovingImage));%fixed image size
cmmIregistered = imwarp(centerMarkMovingImage,tform,'OutputView',Rfixed);%centermarked moving image registered
% C1 = imfuse(centerMarkFixedImage,centerMarkMovingImage,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
% C2 = imfuse(centerMarkFixedImage,cmmIregistered,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
C1 = cat(3,centerMarkFixedImage,centerMarkMovingImage,centerMarkFixedImage*0);
C2 = cat(3,centerMarkFixedImage,cmmIregistered,centerMarkFixedImage*0);
figure, 
imshowpair(C1,C2,'montage')
title('befor correction (left)- afer correction (right)')
% C1 = imoverlay(false(size(centerMarkFixedImage)),centerMarkFixedImage,[1,0,0]);
% C2 = imoverlay(false(size(centerMarkMovingImage)),centerMarkMovingImage,[0,1,0]);
% figure, imshowpair(C1,C2,'montage');
end

%% bacteriaTracking
function [ CCbranch , bioTree] = bacteriaTracking( preImage , maskImage )
%BACTERIATRACKING Summary of this function goes here
%   Detailed explanation goes here
connectMask = false(size(maskImage,1),size(maskImage,2),2);
connectMask(:,:,1) = preImage;
connectMask(:,:,2) = maskImage;
CCbranch = bwconncomp(connectMask(:,:,:),26);
bioTree = twoFrameConnect(CCbranch,1);
bioTree=bioTreeMeasure(bioTree,0,CCbranch.ImageSize(1),CCbranch.ImageSize(2));

end

function bioTree=twoFrameConnect(CC,startFrame) %initial conenct i and i+1 frame as well as intial the data structure
xSize=CC.ImageSize(1);
ySize=CC.ImageSize(2);
allSize=xSize*ySize;
endFrame=startFrame+1;
bioTree{startFrame}.leavies=[];
bioTree{startFrame}.node=[];
bioTree{startFrame}.root=[];
bioTree{endFrame}.root=[];
bioTree{endFrame}.node=[];
bioTree{endFrame}.leavies=[];
count1=1;count2=1;count3=1;count4=1;
for iObject=1:CC.NumObjects
    
    pixelFrame=fix(CC.PixelIdxList{iObject}./allSize);
    startFramePixel=CC.PixelIdxList{iObject}(pixelFrame==0);
    endFramePixel=CC.PixelIdxList{iObject}(pixelFrame==1);
    
    if isempty(startFramePixel)
        bioTree{1,endFrame}.root{count1}.rootPixelDetail=endFramePixel-allSize;
        bioTree{1,endFrame}.root{count1}.rootPixelNum=numel(endFramePixel-allSize);
        bioTree{1,endFrame}.root{count1}.is2Node=false;
        bioTree{1,endFrame}.root{count1}.leafInfo=[endFrame,count4];
        bioTree{1,endFrame}.root{count1}.nodeInfo=[];
        bioTree{1,endFrame}.root{count1}.traceInfo.pixelIdxList{1}=endFramePixel-allSize;
        bioTree{1,endFrame}.leavies{count4}.is2Node=false;
        bioTree{1,endFrame}.leavies{count4}.rootInfo=[endFrame,count1];
        bioTree{1,endFrame}.leavies{count4}.nodeInfo=[];
        bioTree{1,endFrame}.leavies{count4}.leaviesPixelDetail=endFramePixel-allSize;
         bioTree{1,endFrame}.leavies{count4}.leaviesPixelNum=numel(endFramePixel-allSize);
        count1=count1+1;
        count4=count4+1;
        continue;
    end
    
    if isempty(endFramePixel)
        bioTree{1,startFrame}.root{count3}.is2Node=false;
        bioTree{1,startFrame}.root{count3}.leafInfo=[startFrame,count2];
        bioTree{1,startFrame}.root{count3}.nodeInfo=[];
        bioTree{1,startFrame}.root{count3}.rootPixelDetail= startFramePixel;
        bioTree{1,startFrame}.root{count3}.rootPixelNum= numel( startFramePixel);
        bioTree{1,startFrame}.root{count3}.traceInfo.pixelIdxList{1}=startFramePixel;
        bioTree{1,startFrame}.leavies{count2}.is2Node=false;
        bioTree{1,startFrame}.leavies{count2}.nodeInfo=[];
        bioTree{1,startFrame}.leavies{count2}.rootInfo=[startFrame,count3];
        bioTree{1,startFrame}.leavies{count2}.leaviesPixelDetail=startFramePixel;
         bioTree{1,startFrame}.leavies{count2}.leaviesPixelNum=numel( startFramePixel);
        count2=count2+1;
        count3=count3+1;
        continue;
    end
    
    if (~isempty(startFramePixel))&&(~isempty(endFramePixel))
        bioTree{1,startFrame}.root{count3}.is2Node=false;
        bioTree{1,startFrame}.root{count3}.leafInfo=[endFrame,count4];
        bioTree{1,startFrame}.root{count3}.nodeInfo=[];
        bioTree{1,startFrame}.root{count3}.rootPixelDetail=startFramePixel;
        bioTree{1,startFrame}.root{count3}.rootPixelNum=numel(startFramePixel);
        bioTree{1,startFrame}.root{count3}.traceInfo.pixelIdxList{1}=startFramePixel;
        bioTree{1,startFrame}.root{count3}.traceInfo.pixelIdxList{2}=endFramePixel-allSize;
        bioTree{1,endFrame}.leavies{count4}.is2Node=false;
        bioTree{1,endFrame}.leavies{count4}.rootInfo=[startFrame,count3];
        bioTree{1,endFrame}.leavies{count4}.nodeInfo=[];
        bioTree{1,endFrame}.leavies{count4}.leaviesPixelDetail=endFramePixel-allSize;
        bioTree{1,endFrame}.leavies{count4}.leaviesPixelNum=numel(endFramePixel-allSize);
        count3=count3+1;
        count4=count4+1;
        continue;
    end    
end
end

