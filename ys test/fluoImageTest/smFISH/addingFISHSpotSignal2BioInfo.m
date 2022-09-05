function [bioInfo] = addingFISHSpotSignal2BioInfo(dirFile)
%   用于smFIHS中荧光亮点的光强值获取
%   参考smFISHSpotsSignalTest.mlx
%   目前不具有普适性 mScarletI 通道为smFISH通道
%   Shuai Yang 2022/5/18
%   给通过oneShot protocal获得bioInfo中添加 int_FISHspts

xPosLimit = [60,2000];
yPosLimit = [60,2000];

BG_FISHInt = 120; % FISH 背景荧光光强 backGround % 可自定义
thre_FISHInt = 140; %Fish spot 荧光阈值 threshold % 可自定义

load([dirFile,filesep,'result_basic/bioInfo.mat']);
parfor iIm = 1:numel(bioInfo)

    if isempty(bioInfo(iIm).Centroid)
        continue
    end

    bioInfo(iIm).BG_FISHInt = BG_FISHInt;
    bioInfo(iIm).thre_FISHInt = thre_FISHInt;

    fieldIdx = bioInfo(iIm).fieldIdx;
    frameIdx = bioInfo(iIm).frameIdx;

    dirField = [dirFile,'\field',num2str(fieldIdx,'%04d')];
    % FISH stain
    im_FISH = imread([dirField,'\mScarletI\imagemScarletI',num2str(frameIdx,'%05d'),'.tif']);
    im_FISH = imgaussfilt(im_FISH,1);% sigma =1 size ：2*ceil(2*sigma)+1

    % segmentation mask
    dirMask = [dirField,'\segmentation'];
    im_mask = imread([dirMask,'\imagesfGFP00001.tif']);

    posTF1 = (bioInfo(iIm).Centroid(:,1)>= xPosLimit(1)&bioInfo(iIm).Centroid(:,1)<= xPosLimit(2));
    posTF2 = (bioInfo(iIm).Centroid(:,2)>= yPosLimit(1)&bioInfo(iIm).Centroid(:,2)<= yPosLimit(2));
    posTF = posTF1 & posTF2;
    cell_PixelIdxList = bioInfo(iIm).PixelIdxList(posTF);
    bioInfo(iIm).posTF = posTF;
    celluar = smFISHSpotsSignalGet(im_FISH,cell_PixelIdxList,BG_FISHInt,thre_FISHInt);
    bioInfo(iIm).int_FISHspts = celluar.int_spts;
   
end

save([dirFile,filesep,'result_basic/bioInfo.mat'],'bioInfo');

end

function celluar = smFISHSpotsSignalGet(im_FISH,cell_PixelIdxList,BG_FISHInt,thre_FISHInt)
% 参考smFISHSpotsSignalTest.mlx
% Shuai Yang 2022/5/16

imsz = size(im_FISH);
celluar = {};
celluar.int_spts = [];
se = strel('disk',5);
for iCell = 1:numel(cell_PixelIdxList)
    cellMask = false(imsz);
    cellMask(cell_PixelIdxList{iCell}) = true;
    cellMask = imdilate(cellMask,se);
    stats = regionprops(cellMask,'BoundingBox');
    BB = stats.BoundingBox; %cell BoundingBox
    % 取出局部图像
    x_leftup = floor(BB(2)-3);% 左上角x y坐标 向上拓展3pixel
    y_leftup = floor(BB(1)-3);
    x_RightDw = x_leftup + BB(4)+7;% 右下角x y坐标 % 向下拓展3pixel
    y_RightDw = y_leftup + BB(3)+7; 
    cellMask_local = cellMask(x_leftup:x_RightDw,y_leftup:y_RightDw);
    cellFISH_local = im_FISH(x_leftup:x_RightDw,y_leftup:y_RightDw);
    % imshowMerge_ys(cellMask_local,cellFISH_local);
    spots = imregionalmax(cellFISH_local,8);
    spots = spots&cellMask_local;
    

    CC = bwconncomp(spots);% 找出spots中的连通区域 多个连通只留一个
    for iCC = 1:CC.NumObjects 
        if numel(CC.PixelIdxList{iCC})>1
            cc_FISHInts = double(cellFISH_local(CC.PixelIdxList{iCC}));
            k = find(cc_FISHInts == max(cc_FISHInts));
            CC.PixelIdxList{iCC} = CC.PixelIdxList{iCC}(k(1));
        end
    end
    spots_array = [CC.PixelIdxList{1:CC.NumObjects}]';
    % new spots generation
    spots = false(size(spots));
    spots(spots_array) = true;


    % 对每个spots 进行分析
    spts_Int = [];
    for j = numel(spots_array):-1:1
        spot_mask = false(size(cellMask_local));
        spot_mask(spots_array(j)) = true ;
        spot_neighbor = imdilate(spot_mask,strel('disk', 4));
        % 通过imdilalte 向外扩展3 pixel
        % imshowMerge_ys(spot_neighbor,cellFISH_local);

        neighbor_FISHInts = double(cellFISH_local(spot_neighbor));% 获取FISH intensity
        spot_FISHInt = double(cellFISH_local(spot_mask));% spot 中心点FISH intensity

        % 如果 spot 的中心点小于阈值不要
        if spot_FISHInt < thre_FISHInt
            spots_array(j) = [];% 倒序等于[]，相当于直接删除此数值
            continue
        end
        oneSpt_Int = mean(neighbor_FISHInts)-BG_FISHInt;
        spts_Int = [spts_Int,oneSpt_Int];
    end
    celluar.int_spts{iCell} = spts_Int;

end

end 