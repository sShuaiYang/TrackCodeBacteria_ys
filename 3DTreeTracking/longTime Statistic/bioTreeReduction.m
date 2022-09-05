function bioTree=bioTreeReduction(bioTree)
for iframe=1:size(bioTree,2)
    rootNum=[];
    leafNum=[];
    for iRoot=1:size(bioTree{iframe}.root,2)
        if bioTree{iframe}.root{iRoot}.is2Node==0
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            if leafInfo(1)==iframe
                rootNum=[rootNum;iRoot];
                leafNum=[leafNum;leafInfo(2)];
            end
        end
    end
    bioTree{iframe}.root(rootNum)=[];
    for iRoot=1:size(bioTree{iframe}.root,2)
        if bioTree{iframe}.root{iRoot}.is2Node==0
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iRoot];
        end
        if bioTree{iframe}.root{iRoot}.is2Node==1
            nodeInfo=bioTree{iframe}.root{iRoot}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[iframe,iRoot];
        end
    end
    bioTree{iframe}.leavies(leafNum)=[];
    for iLeaf=1:size(bioTree{iframe}.leavies,2)
        if bioTree{iframe}.leavies{iLeaf}.is2Node==0
            rootInfo=bioTree{iframe}.leavies{iLeaf}.rootInfo;
            bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=[iframe,iLeaf];
        end
        if bioTree{iframe}.leavies{iLeaf}.is2Node==1
            nodeInfo=bioTree{iframe}.leavies{iLeaf}.nodeInfo;
            bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[iframe,iLeaf];
        end
    end
end