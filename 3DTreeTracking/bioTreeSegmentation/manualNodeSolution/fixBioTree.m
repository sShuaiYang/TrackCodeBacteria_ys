function bioTree=bioTreeDeleteLeaf(bioTree,leafInfo)
iframe=leafInfo(1);
bioTree{iframe}.leavies(leafInfo(2))=[];
for iLeaf=1:size(bioTree{iframe}.leavies,2)
    is2Node=bioTree{iframe}.leavies{iLeaf}.is2Node;
    if is2Node==1
        nodeInfo=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[iframe,iLeaf];
    end
    if is2Node==0
        rootInfo=bioTree{iframe}.leavies{iLeaf}.rootInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=[iframe,iLeaf];
    end
end
