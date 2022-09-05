function [outbw] = removeRegionsByEllipseFitting(bw)
%根据 椭圆内细菌区域的面积/椭圆的面积为阈值进行判断
%小于阈值 说明椭圆拟合不好 不认为是一个细菌 删除
%参考 maskRegionsEllipseFit_subCell.mlx
%maskRegionsEllipseFit.mlx
%凸包中区域内像素所占的比例，以标量形式返回。计算为 Area/ConvexArea
% Shuai Yang 2021.07.20

ellInPixelsPropThreshold = 0.7;
L = bwlabel(bw);
s = regionprops(bw, 'Area');
mask_ell = false(size(bw));
% s_Ells = struct([]);
% tic;
for iCell = 1:numel(s)
    Lcell= +(L == iCell); % + converts logical to double
    Lcell(Lcell == 1)= iCell;
    % disp(num2str(iCell));
    % extract subimages
    [fx, fy] = find(Lcell);
    extra= 8;
    xmin = max(min(fx) - extra, 1);
    xmax = min(max(fx) + extra, size(Lcell,1));
    ymin = max(min(fy) - extra, 1);
    ymax = min(max(fy) + extra, size(Lcell,2));
    subcell = Lcell(xmin:xmax, ymin:ymax);
    
    s_subcell = regionprops(logical(subcell),'MajorAxisLength',...
        'MinorAxisLength', 'Orientation','PixelIdxList','Centroid');
    
    subCell_Ell = false(size(subcell));% each cell ellipse
    
    % ellipse fitting
    angle = s_subcell.Orientation; % 角度
    angle = angle/180*pi;
    r = 0:0.01:2*pi; % 转一圈
    
    a = s_subcell.MajorAxisLength/2;%semi_MajorAxis
    b = s_subcell.MinorAxisLength/2;%semi_minorAxis
    p = [(a*cos(r))' (b*sin(r))']; % 未旋转
    alpha = [cos(angle) -sin(angle)
        sin(angle) cos(angle)];
    p1 = p*alpha; % 旋转后
    p1(:,1) =  p1(:,1) + s_subcell.Centroid(1);
    p1(:,2) =  p1(:,2) + s_subcell.Centroid(2);
    xIdx = round(p1(:,1));
    yIdx = round(p1(:,2));
    
    % 如果超出索引 说明细菌椭圆拟合的piexel数值处于边缘 不要
    if max(xIdx) > size(subcell,2) || max(yIdx) > size(subcell,1) ...
            || any(min(xIdx(:)) < 1) ||  any(min(yIdx(:)) < 1)
        %         mask_ell(xmin:xmax, ymin:ymax) = ...
        %             or(mask_ell(xmin:xmax, ymin:ymax),subCell_Ell);
        continue
    else
        ind = sub2ind(size(subcell),yIdx(:),xIdx(:));
    end

    subCell_Ell(ind) = true;
    subCell_EllMask = imfill(subCell_Ell,'holes');
    
    s_subCellEll = regionprops(subCell_EllMask, 'Area');    
%     s_Ells(iCell).Area = s_subCellEll.Area;
    
    ellIn = and(subCell_EllMask,subcell);
    %椭圆内cell region 面积占比
    ellOut = and(~subCell_EllMask,subcell);
    %椭圆外cell region面积占比

    s(iCell).ellArea = s_subCellEll.Area;
    s(iCell).ellInPixels = sum(ellIn(:));
    s(iCell).ellOutPixels = sum(ellOut(:));
    s(iCell).ellInPixelsProportion = s(iCell).ellInPixels/s(iCell).ellArea;
    
    mask_ell(xmin:xmax, ymin:ymax) = or(mask_ell(xmin:xmax, ymin:ymax),subCell_Ell);
    
end
% toc;
% figure,imshowpair(bw,mask_ell)
% ell = patch(p1(:,1),p1(:,2),[255 225 225]/255,'EdgeColor',[1 0 0]);

% case for  empty image
% case for all cells are edge pixels, not subjected to ellipseFitting
if numel(s) == 0 || ~isfield(s,'ellInPixelsProportion')                         
    outbw = bw;
    return
end

% 椭圆内细菌区域的面积/椭圆的面积
fwrong = find([s.ellInPixelsProportion] < ellInPixelsPropThreshold);
% wrongbw = ismember(L, fwrong);
goodones = setdiff(1:max(max(L)), fwrong);
goodbw = ismember(L, goodones);%没有fwrong剩下的bw图像
outbw = goodbw;
end
