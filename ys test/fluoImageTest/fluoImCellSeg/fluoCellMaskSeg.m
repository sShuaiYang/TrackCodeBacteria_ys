function outbw = fluoCellMaskSeg(bw,p)
% shuai Yang 2020.10.12
% bw: cell masks of one fluo image 
bw1 = bw;
L1 = bwlabel(bw1);
% imshowlabel_ys(L,I0);
% imshowlabel_ys(L1);
%% 1. remove small regions
% based on Area 或者 majorAxisLength
% disp('remove small regions');

% r1 = regionprops(bw1,'Area');
% fsmall = find([r1.Area] < p.minAreaThreshold);
r1 = regionprops(bw1,'Area','MajorAxisLength');
fsmall = [find([r1.Area] < p.minAreaThreshold),...
find([r1.MajorAxisLength] < p.minCellLength)];
fsmall = unique(fsmall,'sorted');

if ~isempty(fsmall)
    for iRegion = 1:numel(fsmall)
        L1(L1==fsmall(iRegion))= 0;
    end    
end

%% 2. big regions split, large cells judged by Area
% disp('split big regions ');
bw2 = logical(L1);
L2 = bwlabel(bw2);
r2 = regionprops(bw2,'Area','Centroid');
% 先找出面积比较大的区域 然后用watershed进行分割
%根据length进行判断分割结果
% segCellTest_Watershedtransform.mlx
fbig = find([r2.Area] > p.maxAreaThreshold);
goodones = setdiff(1:max(max(L2)), fbig);
goodbw = ismember(L2, goodones);%没有fbig剩下的bw图像
% imshow(goodbw)
% bigbw = ismember(L2, fbig);
% imshow(bigbw)
% rbig = regionprops(bigbw,'PixelIdxList','MajoraxisLength');
[outcell] = doSegWatershed(L2,fbig,p);
% [outcell] = removeRegionsBySolidity(outcell,p);
bw3 = or(goodbw,outcell);

%% 3.kinky split, kinky cells judged by solidity
% disp('split kinky cells ');
L3 = bwlabel(bw3);
% imshowlabel_ys(L3);
r3 = regionprops(bw3,'solidity');
fkinky = find(([r3.Solidity] < p.minSolidity));
goodones = setdiff(1:max(max(L3)), fkinky);
goodbw = ismember(L3, goodones);%没有fkinky剩下的bw图像
% imshow(goodbw)
% kinkybw = ismember(L3, fkinky);
% rkinky = regionprops(kinkybw,'ConvexImage','ConvexHull','PixelIdxList');
[outcell] = doSegWatershed(L3,fkinky,p);
% [outcell] = removeRegionsBySolidity(outcell,p);
bw4 = or(goodbw,outcell);

%% 4.fat split, kinky cells judged by width
% disp('split fat cells ');
L4 = bwlabel(bw4);
% imshowlabel_ys(L4);
r4 = regionprops(bw4,'MinorAxisLength');
ffat= find(([r4.MinorAxisLength] > p.maxCellWidth));
goodones = setdiff(1:max(max(L4)), ffat);
goodbw = ismember(L4, goodones);%没有ffat剩下的bw图像
% imshow(goodbw)
% fatbw = ismember(L4, ffat);
% figure,imshow(fatbw)
% rfat = regionprops(fatbw,'MinorAxisLength','PixelIdxList');
[outcell] = doSegWatershed(L4,ffat,p);
% [outcell] = removeRegionsBySolidity(outcell,p);
bw5 = or(goodbw,outcell);

%% 5.long split, long cells judged by length
% disp('split long cells ');
L5 = bwlabel(bw5);
% imshowlabel_ys(L5);
r5 = regionprops(bw5,'MajorAxisLength');
flong= find(([r5.MajorAxisLength] > p.maxCellLength));
goodones = setdiff(1:max(max(L5)), flong);
goodbw = ismember(L5, goodones);%没有ffat剩下的bw图像
% imshow(goodbw)
% longbw = ismember(L5, flong);
% figure,imshow(longbw)
% rlong = regionprops(longbw,'MajorAxisLength','PixelIdxList');
[outcell] = doSegWatershed(L5,flong,p,p.conservativeH);
% [outcell] = removeRegionsBySolidity(outcell,p);
bw6 = or(goodbw,outcell);

%% 6. cell thinkness segmentation 
%如果细菌连通区域具有<=2pixel的薄弱点 认为可以分开
L6 = bwlabel(bw6);
r6 = regionprops(bw6,'MajorAxisLength');
flong= find(([r6.MajorAxisLength] > p.maxCellLength));
goodones = setdiff(1:max(max(L6)), flong);
goodbw = ismember(L6, goodones);%没有ffat剩下的bw图像
% imshow(goodbw)
% longbw = ismember(L6, flong);
% figure,imshow(longbw)
% rlong = regionprops(longbw,'MajorAxisLength','PixelIdxList');
[outcell] = doSegThickness(L6,flong,p);
% [outcell] = removeRegionsBySolidity(outcell,p);
outbw = or(goodbw,outcell);
%%
[outbw] = removeRegionsBySolidity(outbw,p);
[outbw] = removeRegionsByEllipseFitting(outbw);
end
%% Shuai Yang 2021.07.08
function [outbw] = removeRegionsBySolidity(bw,p)
%根据Solidity判断是否有细菌形态
%Solidity
%凸包中区域内像素所占的比例，以标量形式返回。计算为 Area/ConvexArea

L = bwlabel(bw);
r = regionprops(bw,'solidity');
fkinky = find([r.Solidity] < p.minSolidity);

if ~isempty(fkinky)
    for iRegion = 1:numel(fkinky)
        L(L == fkinky(iRegion)) = 0;
    end    
end
outbw = logical(L);
end