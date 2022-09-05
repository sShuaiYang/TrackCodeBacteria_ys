function bioTree1=bioTreeConnecter(bioTree1,bioTree2,connectPoint,xSize,ySize)
%link the leavies in i frame @bioTree1 and roots in the i+1 frame @bioTree2. (one to one link)
i=1;
%need input correct the iamge size in this function 
while i<=size(bioTree1{connectPoint}.leavies,2)
    ileaf=bioTree1{connectPoint}.leavies{i}.leaviesPixelDetail;
    for j=1:size(bioTree2{connectPoint}.root,2)
        jroot=bioTree2{connectPoint}.root{j}.rootPixelDetail;
        if isequal(ileaf,jroot)
            [bioTree1,bioTree2]=linearLinker(bioTree1,bioTree2,connectPoint,i,j);
            i=i-1;
            break;
        end
    end
    i=i+1;
end

if size(bioTree1{connectPoint}.leavies,2)>0
    nodeLeafInTree1=bioTree1{connectPoint}.leavies;
    nodeRootInTree2=bioTree2{connectPoint}.root;
    kNode=createNodeMask(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,xSize,ySize);
    
    for iNode=1:size(kNode,2)
        [nodeIn,nodeOut]=searchingNode(nodeLeafInTree1,nodeRootInTree2,kNode,iNode);
        [bioTree1,bioTree2]=nodeLinker(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,nodeIn,nodeOut,iNode);
    end
end

bioTree1=connectTwoTree(bioTree1,bioTree2,connectPoint);
end

function [bioTree1,bioTree2] =linearLinker(bioTree1,bioTree2,connectPoint,i,j)

if  (bioTree1{connectPoint}.leavies{i}.is2Node==false)&&(bioTree2{connectPoint}.root{j}.is2Node==false)
    rootInTree1=bioTree1{connectPoint}.leavies{i}.rootInfo;
    leafInTree2=bioTree2{connectPoint}.root{j}.leafInfo;
    sizeTrace1=size(bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=leafInTree2;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=rootInTree1;
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

if  (bioTree1{connectPoint}.leavies{i}.is2Node==false)&&(bioTree2{connectPoint}.root{j}.is2Node==true)
    rootInTree1=bioTree1{connectPoint}.leavies{i}.rootInfo;
    nodeInTree2=bioTree2{connectPoint}.root{j}.nodeInfo;
    sizeTrace1=size(bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.is2Node=true;
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=[];
    bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.nodeInfo=nodeInTree2;
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=rootInTree1;
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

if  (bioTree1{connectPoint}.leavies{i}.is2Node==true)&&(bioTree2{connectPoint}.root{j}.is2Node==false)
    nodeInTree1=bioTree1{connectPoint}.leavies{i}.nodeInfo;
    leafInTree2=bioTree2{connectPoint}.root{j}.leafInfo;
    sizeTrace1=size(bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=leafInTree2;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.is2Node=true;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.nodeInfo=nodeInTree1;
    bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=[];
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;
end

if (bioTree1{connectPoint}.leavies{i}.is2Node==true)&&(bioTree2{connectPoint}.root{j}.is2Node==true)
    nodeInTree1=bioTree1{connectPoint}.leavies{i}.nodeInfo;
    nodeInTree2=bioTree2{connectPoint}.root{j}.nodeInfo;
    sizeTrace1=size(bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2);
    sizeTrace2=size(bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList,2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.pixelIdxList(1:sizeTrace2);
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.is2Node=true;
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=[];
    bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.nodeInfo=nodeInTree2;
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.isNode=true;
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=[];
    bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.nodeInfo=nodeInTree1;
    bioTree1{connectPoint}.leavies(i)=[];
    bioTree2{connectPoint}.root(j)=[];
    return;  
end

end 
function kNode=createNodeMask(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,xSize,ySize) 
allSize=xSize*ySize;
nS1=false(xSize,ySize);
nS2=false(xSize,ySize);
nS3=false(xSize,ySize);

for ileaf=1:size(nodeLeafInTree1,2)
    pixelIdxList2=bioTree1{connectPoint}.leavies{ileaf}.leaviesPixelDetail;
    nS2(pixelIdxList2(pixelIdxList2>0))=1;
    clear pixelIdxList2;
    
    if nodeLeafInTree1{ileaf}.is2Node==false
        rootInTree1=nodeLeafInTree1{ileaf}.rootInfo;
        if size(bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList,2)>1
            pixelIdxList1=bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.pixelIdxList{end-1};
            nS1(pixelIdxList1)=1;
            clear pixelIdxList1;
        end
        continue;
    end
    if nodeLeafInTree1{ileaf}.is2Node==true
        nodeInTree1=nodeLeafInTree1{ileaf}.nodeInfo;
        if size(bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList,2)>1
            pixelIdxList1=bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.pixelIdxList{end-1};
            nS1(pixelIdxList1)=1;
            clear pixelIdxList1;
        end
        continue;
    end
end

for iroot=1:size(nodeRootInTree2,2)
    
    if size(bioTree2{connectPoint}.root{iroot}.traceInfo.pixelIdxList,2)>1
        pixelIdxList3=bioTree2{connectPoint}.root{iroot}.traceInfo.pixelIdxList{2};
        nS3(pixelIdxList3(pixelIdxList3>0))=1;
        clear pixelIdxList3;
    end
    
end

nodeStacks(:,:,1)=nS1;
nodeStacks(:,:,2)=nS2;
nodeStacks(:,:,3)=nS3;
CCnode=bwconncomp(nodeStacks,26);

for iObject=1:CCnode.NumObjects
    pixelFrame=fix(CCnode.PixelIdxList{iObject}./allSize);
    kNode{iObject}=CCnode.PixelIdxList{iObject}(pixelFrame==1)-allSize;
end
end 
function [nodeIn,nodeOut]=searchingNode(nodeLeafInTree1,nodeRootInTree2,kNode,iNode)
nodeIn=[];
nodeOut=[];
for iNodeIn=1:size(nodeLeafInTree1,2)
    if ismember(nodeLeafInTree1{iNodeIn}.leaviesPixelDetail,kNode{iNode})
        nodeIn=[nodeIn,iNodeIn];
    end
end

for iNodeOut=1:size(nodeRootInTree2,2)
    if ismember(nodeRootInTree2{iNodeOut}.rootPixelDetail,kNode{iNode})
        nodeOut=[nodeOut,iNodeOut];
    end
end

end
function [bioTree1,bioTree2]=nodeLinker(bioTree1,bioTree2,nodeLeafInTree1,nodeRootInTree2,connectPoint,nodeIn,nodeOut,iNode)

for iNodeIn=1:size(nodeIn,2)
    if  nodeLeafInTree1{nodeIn(iNodeIn)}.is2Node==false
        rootInTree1=nodeLeafInTree1{nodeIn(iNodeIn)}.rootInfo;
        bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.is2Node=true;
        bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.nodeInfo=[connectPoint,iNode,iNodeIn];
        bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.leafInfo=[];
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.isNode=false;
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.rootInfo=rootInTree1;
        clear rootInTree1;
        continue;
    end
    
    if  nodeLeafInTree1{nodeIn(iNodeIn)}.is2Node==true
        nodeInTree1=nodeLeafInTree1{nodeIn(iNodeIn)}.nodeInfo;
        bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.is2Node=true;
        bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.leafInfo=[];
        bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.nodeInfo=[connectPoint,iNode,iNodeIn];
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.isNode=true;
        bioTree2{connectPoint}.node{iNode}.In{iNodeIn}.nodeInfo=nodeInTree1;
        clear nodeInTree1;
        continue;
    end
end

for iNodeOut=1:size(nodeOut,2)
    if nodeRootInTree2{nodeOut(iNodeOut)}.is2Node==false
        leafInTree2=nodeRootInTree2{nodeOut(iNodeOut)}.leafInfo;
        bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.is2Node=true;
        bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.rootInfo=[];
        bioTree2{leafInTree2(1)}.leavies{leafInTree2(2)}.nodeInfo=[connectPoint,iNode,iNodeOut];
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.is2Node=false;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.leafInfo=leafInTree2;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.traceInfo=nodeRootInTree2{nodeOut(iNodeOut)}.traceInfo;
        clear leafInTree2;
        continue;
    end
    
    if nodeRootInTree2{nodeOut(iNodeOut)}.is2Node==true
        nodeInTree2=nodeRootInTree2{nodeOut(iNodeOut)}.nodeInfo;
        bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.isNode=true;
        bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.rootInfo=[];
        bioTree2{nodeInTree2(1)}.node{nodeInTree2(2)}.In{nodeInTree2(3)}.nodeInfo=[connectPoint,iNode,iNodeOut];
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.is2Node=true;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.nodeInfo=nodeInTree2;
        bioTree2{connectPoint}.node{iNode}.Out{iNodeOut}.traceInfo=nodeRootInTree2{nodeOut(iNodeOut)}.traceInfo;
        clear nodeInTree2;
        continue;
    end
    
end

end
function bioTree1=connectTwoTree(bioTree1,bioTree2,connectPoint)
bioTree1{connectPoint}.leavies=[];
bioTree2{connectPoint}.root=[];
bioTree1{connectPoint}.leavies=bioTree2{connectPoint}.leavies;
bioTree1{connectPoint}.node=bioTree2{connectPoint}.node;
bioTree1(connectPoint+1:1:size(bioTree2,2))=bioTree2(connectPoint+1:1:size(bioTree2,2));
end