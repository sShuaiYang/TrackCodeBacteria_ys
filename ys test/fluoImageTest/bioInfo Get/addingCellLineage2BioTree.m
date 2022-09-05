function [bioTree] = addingCellLineage2BioTree( bioTree,dirField)
%将细菌的lineage信息添加到bioTree中
% cellLineage 结构
fieldIdx = str2double(dirField(end-3:end));
for iFrame = 1:size(bioTree,2)
    if isempty(bioTree{iFrame})
        continue
    end
    
    for iRoot = 1:numel(bioTree{iFrame}.root)
        cellLineage.birthField =fieldIdx;
        cellLineage.birthFrame = iFrame;
        cellLineage.birthRoot = iRoot;
        cellLineage.label = '1';
        bioTree{iFrame}.root{iRoot}.cellLineage = cellLineage;
        
        if bioTree{iFrame}.root{iRoot}.is2Node==1
            agePos = zeros(3,2);
            %agePos细菌分裂前的最后位置，第1,2,3列 0端,1端，质心位置
            traceInfo = bioTree{iFrame}.root{iRoot}.traceInfo;
            if isempty(traceInfo.measurment)
                continue
            end
            agePos(1,:) = traceInfo.measurment{end}.p1Position;%0 端 root 人为认定p1为0端，p2为1端
            agePos(2,:) = traceInfo.measurment{end}.p2Position;
            agePos(3,:) = traceInfo.measurment{end}.Centroid;
            nodeInfo = bioTree{iFrame}.root{iRoot}.nodeInfo();
            [bioTree] = generateCellLineageFrNode(bioTree,cellLineage,agePos,nodeInfo);
        end
        
        
    end
    
end
end
%%
function [bioTree] = generateCellLineageFrNode(bioTree,motherCellLineage,agePos,nodeInfo)

%用try尝试，判断node是否存在
try
    bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.In{nodeInfo(1,3)};
catch
    disp('error node')
    return
end

if ~isempty(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)&&...
        numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)>=nodeInfo(1,3)*2&&...
        numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out)==...
        numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.In)*2
    %确保有Out以及1个In 2个Out %判断node为1分2的node
    
    for iOut = nodeInfo(1,3)*2-1:nodeInfo(1,3)*2 %确认1In对应2Out后，找到两个Out位置
        if mod(iOut,2)==1 %iOut 第一个为奇数 第二个为偶数
            traceInfo1 = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.traceInfo;
            pos_daughter1 = zeros(3,2);%获取刚分裂时head tail centroid位置
            pos_daughter1(1,:) = traceInfo1.measurment{1}.p1Position;%
            pos_daughter1(2,:) = traceInfo1.measurment{1}.p2Position;
            pos_daughter1(3,:) = traceInfo1.measurment{1}.Centroid;
        else
            traceInfo2 = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.traceInfo;
            pos_daughter2 = zeros(3,2);
            pos_daughter2(1,:) = traceInfo2.measurment{1}.p1Position;%
            pos_daughter2(2,:) = traceInfo2.measurment{1}.p2Position;
            pos_daughter2(3,:) = traceInfo2.measurment{1}.Centroid;
        end
    end
    
    [lineage_daugher1,lineage_daugher2,agePos_daughter1,agePos_daughter2] =...
        judgeDivisionDaughterCellPoleAge(agePos,motherCellLineage,pos_daughter1,pos_daughter2);
    
    bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut-1}.cellLineage = lineage_daugher1;
    bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.cellLineage = lineage_daugher2;
    % daughter cell 1
    if bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut-1}.is2Node==1
        nextNodeInfo=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut-1}.nodeInfo;
        if isequal(pos_daughter1,agePos_daughter1)% 如果相等0,1就是p1和p2,否则调换
            agePos = zeros(3,2);%获下一次分裂前0,1端位置
            agePos(1,:) = traceInfo1.measurment{end}.p1Position;%
            agePos(2,:) = traceInfo1.measurment{end}.p2Position;
            agePos(3,:) = traceInfo1.measurment{end}.Centroid;
        else
            agePos = zeros(3,2);
            agePos(2,:) = traceInfo1.measurment{end}.p1Position;%
            agePos(1,:) = traceInfo1.measurment{end}.p2Position;
            agePos(3,:) = traceInfo1.measurment{end}.Centroid;
        end
        [bioTree] = generateCellLineageFrNode(bioTree,lineage_daugher1,agePos,nextNodeInfo);
    end
     % daughter cell 2
    if bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.is2Node==1
        nextNodeInfo=bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.nodeInfo;
        if isequal(pos_daughter2,agePos_daughter2)% 如果相等0,1就是p1和p2,否则调换
            agePos = zeros(3,2);%获下一次分裂前0,1端位置
            agePos(1,:) = traceInfo2.measurment{end}.p1Position;%
            agePos(2,:) = traceInfo2.measurment{end}.p2Position;
            agePos(3,:) = traceInfo2.measurment{end}.Centroid;
        else
            agePos = zeros(3,2);
            agePos(2,:) = traceInfo2.measurment{end}.p1Position;%
            agePos(1,:) = traceInfo2.measurment{end}.p2Position;
            agePos(3,:) = traceInfo2.measurment{end}.Centroid;
        end
        [bioTree] = generateCellLineageFrNode(bioTree,lineage_daugher2,agePos,nextNodeInfo);
    end
    
end

end
function [lineage_daugher1,lineage_daugher2,agePos_daughter1,agePos_daughter2] =...
    judgeDivisionDaughterCellPoleAge(agePos_mother,lineage_mother,pos_daughter1,pos_daughter2)
lineage_daugher1 = lineage_mother;
lineage_daugher2 = lineage_mother;
p0_mother = agePos_mother(1,:);%0 端位置 mother cell
p1_mother = agePos_mother(2,:);%1 端位置 mother cell
pC_mother = agePos_mother(3,:);%centroid 位置 mother cell
%比较两个daughter cell质心与mothercell p0和p1的距离来判断其来自于mother的 0 或1 端
D1 = pdist ([pos_daughter1(3,:); p0_mother; p1_mother]);% 取第一列与第二列的值
D2 = pdist ([pos_daughter2(3,:); p0_mother; p1_mother]);
if D1(1)  < D1(2)
    lineage_daugher1.label = strcat(lineage_mother.label,'0');
else
    lineage_daugher1.label = strcat(lineage_mother.label,'1');
end
if D2(1)  < D2(2)
    lineage_daugher2.label = strcat(lineage_mother.label,'0');
else
    lineage_daugher2.label = strcat(lineage_mother.label,'1');
end
if strcmp(lineage_daugher1.label,lineage_daugher2.label)
    if D1(1) < D2(1)
        lineage_daugher1.label = strcat(lineage_mother.label,'0');
        lineage_daugher2.label = strcat(lineage_mother.label,'1');
    else
        lineage_daugher1.label = strcat(lineage_mother.label,'1');
        lineage_daugher2.label = strcat(lineage_mother.label,'0');
    end
end
%比较两个daughter cell的两端与mothercell质心的距离来判断daughter cell 的0端和1端
D1 = pdist ([pC_mother; pos_daughter1(1:2,:)]);
D2 = pdist ([pC_mother; pos_daughter2(1:2,:)]);
if D1(1) < D1(2)
    agePos_daughter1 = pos_daughter1;%第一个点就是0端，离mother的质心近
else
    agePos_daughter1 = [pos_daughter1(2,:); pos_daughter1(1,:); pos_daughter1(3,:)];% 1，2行调换位置
end
if D2(1) < D2(2)
    agePos_daughter2 = pos_daughter2;%第一个点就是0端，离mother的质心近
else
    agePos_daughter2 = [pos_daughter2(2,:); pos_daughter2(1,:); pos_daughter2(3,:)];% 1，2行调换位置
end
end