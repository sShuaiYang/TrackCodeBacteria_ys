function testBioTree(bioTree)
for iframe=1:size(bioTree,2)
    for iRoot=1:size(bioTree{iframe}.root,2)
        if bioTree{iframe}.root{iRoot}.is2Node==1
            nodeInfo=bioTree{iframe}.root{iRoot}.nodeInfo;
            if ~(size(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList,2)==nodeInfo(1)-iframe)
                disp(strcat('root','__',num2str(iframe),'__',num2str(iRoot)))
            end
        end
        if bioTree{iframe}.root{iRoot}.is2Node==0;
            leafInfo=bioTree{iframe}.root{iRoot}.leafInfo;
            if ~(size(bioTree{iframe}.root{iRoot}.traceInfo.pixelIdxList,2)==leafInfo(1)-iframe+1)
                disp(strcat('root','__',num2str(iframe),'__',num2str(iRoot)))
            end
        end
    end
    for iNode=1:size(bioTree{iframe}.node,2)
        for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==1
                nodeInfo=bioTree{iframe}.node{iNode}.Out{iOut}.nodeInfo;
                if ~(size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList,2)==nodeInfo(1)-iframe)
                    disp(strcat('node','__',num2str(iframe),'__',num2str(iNode),'_',num2str(iOut)))
                end
            end
            if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==0
                leafInfo=bioTree{iframe}.node{iNode}.Out{iOut}.leafInfo;
%                 disp(iframe)
%                 disp(iNode)
%                 disp(iOut)
                if ~(size(bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList,2)==leafInfo(1)-iframe+1)
                    disp(strcat('node','_',num2str(iframe),'_',num2str(iNode),'_',num2str(iOut)))
                end
            end
        end
    end
end
[bioTree,branchList,~,~]=myBiograph_new2(bioTree);
for i=1:size(bioTree,2)
    for iLeaf=1:size(bioTree{i}.leavies,2);
        pixel=bioTree{i}.leavies{iLeaf}.leaviesPixelDetail;
    end
end
end