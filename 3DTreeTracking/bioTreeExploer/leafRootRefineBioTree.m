function bioTree=leafRootRefineBioTree(bioTree,frameShift)
distThreshold=60;
for iframe=1+frameShift:size(bioTree,2)-1
%     disp(num2str(iframe));
    if  ~isempty(bioTree{iframe}.leavies)&&~isempty(bioTree{iframe+1}.root)
        [distLeafRoot,indexI]=distenceFinder(bioTree,iframe);
        [leafLink,rootLink]=findLinker(distLeafRoot,indexI,distThreshold);
        if ~isempty(leafLink)
            bioTree=refineConnecter(bioTree,iframe,leafLink,rootLink);
            bioTree{iframe}.leavies(leafLink)=[];
            bioTree{iframe+1}.root(rootLink)=[];
            if  ~isempty(bioTree{iframe}.leavies)
                bioTree=correctLeafIndex(bioTree,iframe);
            end
            if ~isempty(bioTree{iframe+1}.root)
                bioTree=correctRootIndex(bioTree,iframe+1);
            end
        end    
    end
end
end

function [distLeafRoot,indexI]=distenceFinder(bioTree,iframe)
leafXY=[];
rootXY=[];
for ileaf=1:size(bioTree{iframe}.leavies,2)
    pos=fastGetProps(bioTree{iframe}.leavies{ileaf}.leaviesPixelDetail,bioTree{1}.imageSize);
    leafXY=[leafXY;pos.Centroid];
end
for jroot=1:size(bioTree{iframe+1}.root,2)
    pos=fastGetProps(bioTree{iframe+1}.root{jroot}.rootPixelDetail,bioTree{1}.imageSize);
    rootXY=[rootXY;pos.Centroid];
end
[distLeafRoot,indexI]=pdist2(leafXY,rootXY,'euclidean','Smallest',1);
end
function [leafLink,rootLink]=findLinker(distLeafRoot,indexI,distThreshold)
leafLink=[];
rootLink=[];
indexI(distLeafRoot>distThreshold)=0;
for ilinker=1:size(indexI,2)
    if indexI(ilinker)>0
        numLinker=size(indexI(indexI==indexI(ilinker)),2);
        if numLinker==1
            leafLink=[leafLink,indexI(ilinker)];
            rootLink=[rootLink,ilinker];
        else
            rootIndexTemp=find(indexI==indexI(ilinker));
            [~,minIndex]=min(distLeafRoot(rootIndexTemp));
            rootIndex=rootIndexTemp(minIndex);
            leafLink=[leafLink,indexI(ilinker)];
            rootLink=[rootLink,rootIndex];
            indexI(indexI==indexI(ilinker))=0;
        end 
    end
end
end
function bioTree=refineConnecter(bioTree,iframe,leafLink,rootLink)
for ilink=1:size(leafLink,2)
    i=leafLink(ilink);
    j=rootLink(ilink);
    if  (bioTree{iframe}.leavies{i}.is2Node==false)&&(bioTree{iframe+1}.root{j}.is2Node==false)
        rootInTree1=bioTree{iframe}.leavies{i}.rootInfo;
        leafInTree2=bioTree{iframe+1}.root{j}.leafInfo;
        sizeTrace1=size(bioTree{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2);
        sizeTrace2=size(bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList,2);
        bioTree{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
%         bioTree{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.measurment(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.measurment(1:sizeTrace2);
        bioTree{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=leafInTree2;
        bioTree{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=rootInTree1;
        continue;
    end
    
    if  (bioTree{iframe}.leavies{i}.is2Node==false)&&(bioTree{iframe+1}.root{j}.is2Node==true)
        rootInTree1=bioTree{iframe}.leavies{i}.rootInfo;
        nodeInTree2=bioTree{iframe+1}.root{j}.nodeInfo;
        sizeTrace1=size(bioTree{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2);
        sizeTrace2=size(bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList,2);
        bioTree{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
%         bioTree{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.measurment(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.measurment(1:sizeTrace2);
        bioTree{rootInTree1(1)}.root{rootInTree1(2)}.is2Node=true;
        bioTree{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=[];
        bioTree{rootInTree1(1)}.root{rootInTree1(2)}.nodeInfo=nodeInTree2;
        bioTree{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=rootInTree1;
        continue;
    end
    
    if  (bioTree{iframe}.leavies{i}.is2Node==true)&&(bioTree{iframe+1}.root{j}.is2Node==false)
        nodeInTree1=bioTree{iframe}.leavies{i}.nodeInfo;
        leafInTree2=bioTree{iframe+1}.root{j}.leafInfo;
        sizeTrace1=size(bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2);
        sizeTrace2=size(bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList,2);
        bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
%         bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.measurment(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.measurment(1:sizeTrace2);
        bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=leafInTree2;
        bioTree{leafInTree2(1)}.leavies{leafInTree2(2)}.is2Node=true;
        bioTree{leafInTree2(1)}.leavies{leafInTree2(2)}.nodeInfo=nodeInTree1;
        bioTree{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=[];
        continue;
    end
    
    if (bioTree{iframe}.leavies{i}.is2Node==true)&&(bioTree{iframe+1}.root{j}.is2Node==true)
        nodeInTree1=bioTree{iframe}.leavies{i}.nodeInfo;
        nodeInTree2=bioTree{iframe+1}.root{j}.nodeInfo;
        sizeTrace1=size(bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2);
        sizeTrace2=size(bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList,2);
        bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
%         bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.measurment(sizeTrace1+1:sizeTrace1+sizeTrace2)=bioTree{iframe+1}.root{j}.traceInfo.measurment(1:sizeTrace2);
        bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.is2Node=true;
        bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=[];
        bioTree{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.nodeInfo=nodeInTree2;
        bioTree{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.isNode=true;
        bioTree{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=[];
        bioTree{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.nodeInfo=nodeInTree1;
        continue;
    end
end
end
function bioTree=correctLeafIndex(bioTree,iframe)

for ileaf=1:size(bioTree{iframe}.leavies,2)
    if bioTree{iframe}.leavies{ileaf}.is2Node==false
        rootInfo=bioTree{iframe}.leavies{ileaf}.rootInfo;
        bioTree{rootInfo(1)}.root{rootInfo(2)}.leafInfo=[iframe,ileaf];
        continue;
    end
    
    if bioTree{iframe}.leavies{ileaf}.is2Node==true
        nodeInfo=bioTree{iframe}.leavies{ileaf}.nodeInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[iframe,ileaf];
        continue;
    end
end

end
function bioTree=correctRootIndex(bioTree,iframe)

for iroot=1:size(bioTree{iframe}.root,2)
    if bioTree{iframe}.root{iroot}.is2Node==false
        leafInfo=bioTree{iframe}.root{iroot}.leafInfo;
        bioTree{leafInfo(1)}.leavies{leafInfo(2)}.rootInfo=[iframe,iroot];
        continue;
    end
    if bioTree{iframe}.root{iroot}.is2Node==true
        nodeInfo=bioTree{iframe}.root{iroot}.nodeInfo;
        bioTree{nodeInfo(1)}.node{nodeInfo(2)}.In{nodeInfo(3)}.rootInfo=[iframe,iroot];
        continue;
    end 
end
end