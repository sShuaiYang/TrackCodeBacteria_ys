function [pedigree] = genealogicalTreeGeneration_Delta(res)
% Shuai Yang
% 2022/01/8
% 用于从delta 程序处理的result来获得谱系信息
% genealogicalTreeGeneration from delta results
pedigree = struct('frameIdx',[ ],'rootIdx',[ ],'tree',{ },'leafNum',[ ]);

lineage = res{1, 1}.lineage;
iRoot = 0;% root细菌的数目
for iCell = 1:numel(lineage)

    if lineage{iCell}.mother == 0 %0 if no mother detected (eg first timepoint) 寻找root
        initial_label = 1;
        iRoot = iRoot+1;
        pedigree(iRoot).rootIdx = iRoot;

        frames = lineage{iCell}.frames;
        frameIdx = frames(1);
        pedigree(iRoot).frameIdx = frameIdx;
        leafCount = 1;
        labelNum = 1;



        daughters = lineage{iCell}.daughters ;
        [~,divPts,daughter_id] = find(daughters); % 获取daughter cell 信息 divison points

        if isempty(daughter_id) % case for没有发生分裂事件
            pedigree(iRoot).tree(leafCount).mother = lineage{iCell}.mother;
            pedigree(iRoot).tree(leafCount).label = initial_label; % root的label定义为1
            pedigree(iRoot).tree(leafCount).linkInfo = [frameIdx,frames(end)];
            pedigree(iRoot).tree(leafCount).cell_id = iCell;
            pedigree(iRoot).tree(leafCount).timer = frameIdx :frames(end);
            leafCount = leafCount+1;
            continue
        else
            daughter_labels = zeros(size(daughter_id));
            for iDiv = 1 :numel(divPts)
                if iDiv == 1
                    label = initial_label;
                end
                pedigree(iRoot).tree(leafCount).mother = lineage{iCell}.mother;
                pedigree(iRoot).tree(leafCount).label = label;
                pedigree(iRoot).tree(leafCount).cell_id = iCell;
                pedigree(iRoot).tree(leafCount).linkInfo = [frameIdx,frames(divPts(iDiv)-1)];

                pedigree(iRoot).tree(leafCount).timer = frameIdx:frames(divPts(iDiv))-1;

                frameIdx = divPts(iDiv);
                leafCount = leafCount+1;
                label = label*2+1;
                daughter_labels(iDiv) = label-1;
            end
            pedigree(iRoot).tree(leafCount).linkInfo = [frameIdx,frames(end)];
            pedigree(iRoot).tree(leafCount).mother = lineage{iCell}.mother;
            pedigree(iRoot).tree(leafCount).label = label;
            pedigree(iRoot).tree(leafCount).cell_id = iCell;
            pedigree(iRoot).tree(leafCount).timer = (frameIdx:frames(end));
            leafCount = leafCount+1;
        end

        tree = pedigree(iRoot).tree;
        [tree,leafCount] = generateTreeFrDaughterCellID(daughter_id,daughter_labels,tree,lineage,leafCount);
        pedigree(iRoot).tree = tree;
        pedigree(iRoot).leafNum = leafCount -1;
    end
end
end


function [tree,leafCount] = generateTreeFrDaughterCellID(daughter_id,daughter_labels,tree,lineage,leafCount)

for iDaughter = 1:numel(daughter_id)

    initial_label = daughter_labels(iDaughter);
    cell_id = daughter_id(iDaughter);
    daughters = lineage{cell_id}.daughters;
    frames = lineage{cell_id}.frames;
    frameIdx = frames(1);
    [~,divPts,next_daughter_id] = find(daughters);

    if isempty(next_daughter_id) % case for没有发生分裂事件
        tree(leafCount).mother = lineage{cell_id}.mother;
        tree(leafCount).label = initial_label;
        tree(leafCount).linkInfo = [frameIdx,frames(end)];
        tree(leafCount).cell_id = cell_id;
        tree(leafCount).timer = frameIdx :frames(end);
        leafCount = leafCount+1;
        continue
    else
        next_daughter_labels = zeros(size(next_daughter_id));
        for iDiv = 1 :numel(divPts)
            if iDiv == 1
                label = initial_label;
            end
            tree(leafCount).mother = lineage{cell_id}.mother;
            tree(leafCount).label = label;
            tree(leafCount).cell_id = cell_id;
            tree(leafCount).linkInfo = [frameIdx,frames(divPts(iDiv))-1];
            tree(leafCount).timer = frameIdx:frames(divPts(iDiv))-1;

            frameIdx = divPts(iDiv);
            leafCount = leafCount+1;
            label = label*2+1;
            next_daughter_labels(iDiv) = label-1;
        end
        tree(leafCount).linkInfo = [frameIdx,frames(end)];
        tree(leafCount).mother = lineage{cell_id}.mother;
        tree(leafCount).label = label;
        tree(leafCount).cell_id = cell_id;
        tree(leafCount).timer = (frameIdx:frames(end));
        leafCount = leafCount+1;
    end
    [tree,leafCount] = generateTreeFrDaughterCellID(next_daughter_id,next_daughter_labels,tree,lineage,leafCount);
end
end
