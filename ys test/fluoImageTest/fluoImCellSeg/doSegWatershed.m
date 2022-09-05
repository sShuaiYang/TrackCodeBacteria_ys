function [outcell] = doSegWatershed(L,fIdx,p,H)
% shuai Yang 2020.10.12
% L = bwlabel(bw);
% fIdx 为bwlabel标记后的bw图像，要分割区域的label index
% H 可以不输入
outcell = false(size(L));
if nargin == 3
    H = p.H;
end
% disp(['need seg cell numbers ',num2str(numel(fIdx))]);
for iCell = 1:numel(fIdx)
    Lcell= +(L == fIdx(iCell)); % + converts logical to double
    Lcell(Lcell == 1)= fIdx(iCell);
    % disp(num2str(iCell));
    % extract subimages
    [fx, fy] = find(Lcell);
    extra= 5;
    xmin = max(min(fx) - extra, 1);
    xmax = min(max(fx) + extra, size(Lcell,1));
    ymin = max(min(fy) - extra, 1);
    ymax = min(max(fy) + extra, size(Lcell,2));
    subcell = Lcell(xmin:xmax, ymin:ymax);
    
    % if the area larger than p.maxAreaThreshold
    % delete the regions. Shuai Yang;2021.07.08
    s_subcell = regionprops(logical(subcell),'Area');
    if s_subcell.Area > p.maxAreaThreshold
        cutcell = false(size(subcell));
        outcell(xmin:xmax, ymin:ymax) = or(outcell(xmin:xmax, ymin:ymax),cutcell);
        continue
    end
    
    %     thin = bwmorph(subcell, 'thin', inf);
    [cutcell,Ld2] = segCellWatershedTransform(subcell,H);
    %     figure,imshowpair(subcell,cutcell,'montage');
    
    
    s_cutcell = regionprops(cutcell,'MajorAxisLength');
    while sum([s_cutcell.MajorAxisLength] < p.minCellLength) && H <= p.maxH
        %based on celllength judge oversegmentaion
        H = H + 0.1;
        [cutcell,Ld2] = segCellWatershedTransform(subcell,H);
        s_cutcell = regionprops(cutcell,'MajorAxisLength');
    end
    %     figure,imshowpair(subcell,cutcell,'montage');
    
    % ridgeLine = Ld2 == 0;% watershed ridge lines
    if ~isequal(logical(subcell),cutcell) % 判断是否被分割 如果没有不再进行分割恢复
        %否则会陷入死循环
        if numel(s_cutcell) > p.maxSegLineNum
            % 如果分割线过多 认为不对 也直接delete
            cutcell = false(size(cutcell));
        else
            cutline = (Ld2 == 0)&subcell;%cut cell lines
            % 如果升高H还会有oversegementation， 人工删除细胞分割线
            while sum([s_cutcell.MajorAxisLength] < p.minCellLength)
                [cutcell,cutline] = overSegRecover(cutcell,cutline);
                s_cutcell = regionprops(cutcell,'MajorAxisLength');
            end
        end
    end
    %     figure,imshowpair(subcell,cutcell,'montage');
    
    % 将cutcell结果复原
    outcell(xmin:xmax, ymin:ymax) = or(outcell(xmin:xmax, ymin:ymax),cutcell);
end
% figure,imshow(outcell)
% bw3 = or(goodbw,outcell);
end
%%
function [cutcell,Ld2] = segCellWatershedTransform(subcell,H)
% shuai Yang 2020.10.12
subcell = logical(subcell);
D = -bwdist(~subcell);
mask = imextendedmin(D,H);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
cutcell = subcell;
cutcell(Ld2 == 0) = 0;

end
%%
function [cutcell,cutline] = overSegRecover(cutcell,cutline)
% shuai Yang 2020.10.12
% 找出被分割的最小的区域 每次添加一次分割线进行overlap复原
L_cutcell = bwlabel(cutcell);

s_cutcell = regionprops(cutcell,'MajorAxisLength','PixelIdxList');
fminRegionIdx = find([s_cutcell.MajorAxisLength] == min([s_cutcell.MajorAxisLength]),1);
s_cutline = regionprops(cutline,'PixelIdxList');

oversegRegion = L_cutcell == fminRegionIdx;

for i = 1:numel(s_cutline)
    overlapRegion = oversegRegion;
    overlapRegion(s_cutline(i).PixelIdxList) = 1;
    s_OLRegion = regionprops(overlapRegion,'PixelIdxList');
    
    if numel(s_OLRegion) == 1
        cutline(s_cutline(i).PixelIdxList) = 0;
        cutcell = cutcell|overlapRegion;
        return
    end
    
end
end

