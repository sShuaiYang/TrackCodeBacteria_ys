function [outcell] = doSegThickness(L,fIdx,p)
% shuai Yang 2020.10.12
% L = bwlabel(bw);
% fIdx 为bwlabel标记后的bw图像，要分割区域的label index
% 根据thickness 的值进行分割 薄弱连接点认为可以分割
outcell = false(size(L));
k = 0;
for iCell = 1:numel(fIdx)
    Lcell= +(L == fIdx(iCell)); % + converts logical to double
    Lcell(Lcell == 1)= fIdx(iCell);
    
    % extract subimages
    [fx, fy] = find(Lcell);
    %此处边缘有限制 对于靠近太靠近边缘的细菌 extra以图像边缘为准
    %无法有extra的预留 
    extra= 7;
    xmin = max(min(fx) - extra, 1);
    xmax = min(max(fx) + extra, size(Lcell,1));
    ymin = max(min(fy) - extra, 1);
    ymax = min(max(fy) + extra, size(Lcell,2));
    subcell = Lcell(xmin:xmax, ymin:ymax);
    %没有extra预留的 细菌不处理; case for 转置切割会出错
    if (max(fx) + extra)  > size(Lcell,1) || ...
            (max(fy) + extra)  > size(Lcell,2)||...
            xmin <= 1 || ymin <= 1
        continue
    end
  
    subcell = logical(subcell);
    % cell thickness calculate
    D = -bwdist(~subcell);
    bw_thin = bwmorph(subcell, 'thin', inf);
    %     bw_perm = bwperim(subcell);
    %     figure,imshowpair(bw_thin,D);
    %     figure,imshowpair(bw_perm,subcell);
    s_thin = regionprops(bw_thin,'PixelIdxList','PixelList');
    
    cellThickness = D([s_thin.PixelIdxList]);
    cellThickness = abs(cellThickness);% actually half of cell thickness +1
    
    fposIdx = find((cellThickness <= p.minCellThickness));%找到<2pixel宽度的点 8连通 包括对角线
    % 薄弱连接点
    
    thinLength = sum(bw_thin(:));%获得核心骨架的长度
    
    cutcell = subcell;
    % 将连接<2pixel的点 进行人为切断分割
    segline  = false(size(subcell));
   
    markPosIdx = 0;% 需要判断薄弱连接点之间的距离 不能小于minCellLength
    for iPts = 1:numel(fposIdx)
        
        if min([fposIdx(iPts)-markPosIdx,thinLength-fposIdx(iPts)]) < p.minCellLength
            continue
        end
        markPosIdx = fposIdx(iPts);
        [row,col] = ind2sub(size(subcell),s_thin.PixelIdxList(fposIdx(iPts)));
        extraPt = 3;%2对应 5*5box；3 7*7box ; 2*extraPt+1
        segline(row-extraPt:row+extraPt,col-extraPt:col+extraPt) ...
            = bw_thin(row-extraPt:row+extraPt,col-extraPt:col+extraPt)';%7*7box转置进行切割 目的：垂直
        %         figure,imshowpair(segline,bw_thin)
        %     segline = imdilate(segline,ones(2));
        segline = bwmorph(segline,'diag');%Uses diagonal fill to eliminate 8-connectivity
        cutcell = and(cutcell,~segline);
        %         figure,imshowpair(cutcell,segline);
    end
    
    CC = bwconncomp(cutcell);
    if CC.NumObjects == 1 %equal to 1;means no seg or segmentaion failed
        cutcell = subcell; % case for cutline too short to break cells
    end
    
    if ~isequal(cutcell,subcell)
        k = k+1;
    end
    
    %     figure,imshow(cutcell)
    outcell(xmin:xmax, ymin:ymax) = or(outcell(xmin:xmax, ymin:ymax),cutcell);
end
% disp(['Thickness Breaking up long cells(', num2str(k),').']);
end