function [pedigreeAll] = genealogicalTreeGeneration(bioTree,dirSave)
% shuai Yang 2020.05.10
%label为i的两个子代的label为2i和2i+1，如1的子代为2和3，2的子代为4和5...

% dirSave = 'C:\Users\XJY\Desktop\genealogicalTreeTest1';

pedigreeAll = struct('frameIdx',[ ],'rootIdx',[ ],'tree',{ },'leafNum',[ ]);
if isempty(bioTree)
    disp('bioTree is emmpty')
    return
end
if isempty(bioTree{1}.root)
    disp('no cells in the 1st frame')
    return
end

for iFrame = 1:size(bioTree,2) 
    if isempty(bioTree{iFrame})
        continue
    end 
    if isempty(bioTree{iFrame}.root)
        continue
    end 
    
    pedigree = struct('frameIdx',[ ],'rootIdx',[ ],'tree',{ },'leafNum',[ ]);
    for iRoot = 1:numel(bioTree{iFrame}.root)
        leafNum = 1;
        labelNum = 1;
        
        pedigree(iRoot).frameIdx = iFrame;
        pedigree(iRoot).rootIdx = iRoot;
        
        pedigree(iRoot).tree(leafNum).label = 1;%第一个root的label定义为1
        pedigree(iRoot).tree(leafNum).rootInfo = [iFrame,iRoot];
        pedigree(iRoot).tree(leafNum).traceInfo = bioTree{iFrame}.root{iRoot}.traceInfo;
        pedigree(iRoot).tree(leafNum).linkInfo = [iFrame,numel(bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList)+iFrame-1];
        pedigree(iRoot).tree(leafNum).is2Node = bioTree{iFrame}.root{iRoot}.is2Node;
        pedigree(iRoot).tree(leafNum).nodeInfo = [];
        pedigree(iRoot).tree(leafNum).leafInfo = [];
        pedigree(iRoot).tree(leafNum).timer = ...
            bioTree{1, 1}.bioTreeTimer(iFrame:numel(bioTree{iFrame}.root{iRoot}.traceInfo.pixelIdxList)+iFrame-1);
       
        % 可以将想要加入的信息加入到pedigree中  Shuai Yang 2022/3/30
        pedigree(iRoot).tree(leafNum).gR = bioTree{iFrame}.root{iRoot}.gR;
        pedigree(iRoot).tree(leafNum).intsfGFP = bioTree{iFrame}.root{iRoot}.intsfGFP;

        if bioTree{iFrame}.root{iRoot}.is2Node==1
            nodeInfo = bioTree{iFrame}.root{iRoot}.nodeInfo();
            pedigree(iRoot).tree(leafNum).nodeInfo = nodeInfo;
            [pedigree,leafNum] = generateTreeInfoFrNode(bioTree,pedigree,iRoot,leafNum,labelNum,nodeInfo);
        else
             pedigree(iRoot).tree(leafNum).leafInfo = bioTree{iFrame}.root{iRoot}.leafInfo;
        end


        pedigree(iRoot).leafNum = leafNum;
    end
    pedigreeAll(end+1:end+numel(pedigree)) = pedigree;
end
pedigreeAll = pedigreeTidy(pedigreeAll);

if isempty (pedigreeAll)
    return
end
[pedigreeAll] = genealogicalTreePlot_3D(pedigreeAll);
[pedigreeAll] = genealogicalTreePlot_Line(pedigreeAll);
[pedigreeAll] = genealogicalTreePlot_Arc(pedigreeAll);
save(strcat(dirSave,'\pedigreeAll.mat'),'pedigreeAll','-v7.3');
% genealogicalTreePlot_MultiType(pedigreeAll,dirSave);
end
%%
function [pedigree,leafNum] = generateTreeInfoFrNode(bioTree,pedigree,iRoot,leafNum,labelNum,nodeInfo)

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

        leafNum = leafNum+1;
        pedigree(iRoot).tree(leafNum).traceInfo = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.traceInfo;
        pedigree(iRoot).tree(leafNum).linkInfo=[nodeInfo(1,1),...
            nodeInfo(1,1)+numel(bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.traceInfo.pixelIdxList)-1];
        pedigree(iRoot).tree(leafNum).timer = bioTree{1, 1}.bioTreeTimer(pedigree(iRoot).tree(leafNum).linkInfo(1)...
            :pedigree(iRoot).tree(leafNum).linkInfo(2));
        pedigree(iRoot).tree(leafNum).is2Node = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.is2Node;
        if bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.is2Node==1
            pedigree(iRoot).tree(leafNum).nodeInfo = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.nodeInfo;
            pedigree(iRoot).tree(leafNum).leafInfo = [];            
        else
            pedigree(iRoot).tree(leafNum).nodeInfo = [];
            pedigree(iRoot).tree(leafNum).leafInfo = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.leafInfo;
        end

        % 可以将想要加入的信息加入到pedigree中  Shuai Yang 2022/3/30
        pedigree(iRoot).tree(leafNum).gR = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.gR;
        pedigree(iRoot).tree(leafNum).intsfGFP = bioTree{nodeInfo(1,1)}.node{nodeInfo(1,2)}.Out{iOut}.intsfGFP;
        
        if mod(iOut,2)==1 %iOut 第一个为奇数 第二个为偶数
            pedigree(iRoot).tree(leafNum).label = labelNum*2;
        else
            pedigree(iRoot).tree(leafNum).label = labelNum*2+1;
        end   
    end

    
    markleafNum=leafNum;
    
    if pedigree(iRoot).tree(markleafNum-1).is2Node==1
        nodeInfo=pedigree(iRoot).tree(markleafNum-1).nodeInfo;
        labelNum=pedigree(iRoot).tree(markleafNum-1).label;
        [pedigree,leafNum]=generateTreeInfoFrNode(bioTree,pedigree,iRoot,leafNum,labelNum,nodeInfo);
    end
    
    if pedigree(iRoot).tree(markleafNum).is2Node==1
        nodeInfo=pedigree(iRoot).tree(markleafNum).nodeInfo;
        labelNum=pedigree(iRoot).tree(markleafNum).label;
        [pedigree,leafNum]=generateTreeInfoFrNode(bioTree,pedigree,iRoot,leafNum,labelNum,nodeInfo);
    end
end
end

% pedigree filter
function pedigreeAll = pedigreeTidy(pedigreeAll)

% 1.将leaf数目少的pedigree进行删除
filterLogical = false (1,numel(pedigreeAll));

for i = 1:numel(pedigreeAll)
    if pedigreeAll(i).leafNum >=4
        filterLogical(i) = true;
    end
    
end
pedigreeAll = pedigreeAll(filterLogical);

%2.根据timer将分裂时间明显不对的pedigree删除
%因为leafNum>=4,此pedigree至少两次分裂，label为1,2,3的细菌定存在
%label等于1的细菌不计算，计算1→2,3的细菌分裂时间
% 如果首次分裂1→2和1→3的时间均很短，则认为错误

thresTime = 15;%min
filterLogical = false (1,numel(pedigreeAll));

for i = 1:numel(pedigreeAll) 
    %首次分裂，位置2和3刚好对应label为2和3的细菌%15min
    if range(pedigreeAll(i).tree(2).timer)>thresTime||...
            range(pedigreeAll(i).tree(3).timer)>thresTime
        filterLogical(i) = true;
    end
end
pedigreeAll = pedigreeAll(filterLogical);
end