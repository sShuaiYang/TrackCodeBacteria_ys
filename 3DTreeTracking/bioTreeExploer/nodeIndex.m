function [bioTree,nodeNum]=nodeIndex(bioTree)
nodeIndex=1;
for iframe=1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                bioTree{iframe}.root{iroot}.isCluster=true;
                bioTree{iframe}.root{iroot}.nodeIndex=nodeIndex;
                nodeIndex=nodeIndex+1;
            else
                bioTree{iframe}.root{iroot}.isCluster=false;
            end
        end
    end
    if ~isempty(bioTree{iframe}.leavies)
        for ileaf=1:size(bioTree{iframe}.leavies,2)
            if bioTree{iframe}.leavies{ileaf}.is2Node==true
                bioTree{iframe}.leavies{ileaf}.isCluster=true;
                bioTree{iframe}.leavies{ileaf}.nodeIndex=nodeIndex;
                nodeIndex=nodeIndex+1;
            else
                bioTree{iframe}.leavies{ileaf}.isCluster=false;
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            bioTree{iframe}.node{iNode}.isCluster=true;
            bioTree{iframe}.node{iNode}.nodeIndex=nodeIndex;
            nodeIndex=nodeIndex+1;
        end
    end
end
nodeNum=nodeIndex-1;
end