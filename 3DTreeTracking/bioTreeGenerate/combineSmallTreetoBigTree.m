function bioTree=combineSmallTreetoBigTree(smallTreeNumber,dirBioTree) % this function cnc be use to connect the bioTree in diffrent time point
bioTree=bioTreeLoader(smallTreeNumber,dirBioTree);
% bioTree=leafRootRefineBioTree(bioTree,0);
% bioTree=bioTreeSizeReduction(bioTree,0);
% bioTree=nodeEdgeReduction(bioTree,0);
% allTreeSave(bioTree);
end
function bioTreeCombine=bioTreeLoader(smallTreeNumber,dirBioTree)
% dirBioTree=uigetdir();
cd(dirBioTree);
bioTreeCombine=[];
for i=1:smallTreeNumber
    disp(strcat('Load small Tree',num2str(i)));
    folderName=strcat('t',num2str(i));
    cd(folderName);
    nameList=dir();
    if size(nameList,1)==3
        fileName=strcat('bioTree',num2str(i));
        treeTemp=load(fileName);
        bioTree=treeTemp.bioTree;
        if i==1
            bioTreeCombine=bioTree;
            xSize=bioTree{1}.imageSize(1);
            ySize=bioTree{1}.imageSize(2);
        else
            connectPoint=size(bioTreeCombine,2);
            bioTreeCombine=bioTreeConnecterNew(bioTreeCombine,bioTree,connectPoint,xSize,ySize);
        end
        cd(dirBioTree);
    end
    if size(nameList,1)>=4
        treeIndex=size(nameList,1)-2;
        fileName1=strcat('bioTree',num2str(i),'-',num2str(treeIndex));
        treeTemp1=load(fileName1);
        bioTree=treeTemp1.bioTree;
        preTreeSize=size(treeTemp1.bioTreeStack,2);
        for j=2:size(nameList,1)-2
            treeIndex=size(nameList,1)-j-1;
            fileName1=strcat('bioTree',num2str(i),'-',num2str(treeIndex));
            treeTemp2=load(fileName1);
            bioTree(1:preTreeSize)=treeTemp2.bioTreeStack;
            preTreeSize=size(bioTree,2);
        end
        if i==1
            bioTreeCombine=bioTree;
            xSize=bioTree{1}.imageSize(1);
            ySize=bioTree{1}.imageSize(2);
        else
            connectPoint=size(bioTreeCombine,2);
            bioTreeCombine=bioTreeConnecterNew(bioTreeCombine,bioTree,connectPoint,xSize,ySize);
        end
        cd(dirBioTree);
    end
end
end
function bioTree1=bioTreeConnecterNew(bioTree1,bioTree2,connectPoint,xSize,ySize)
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
%     bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
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
%     bioTree1{rootInTree1(1)}.root{rootInTree1(2)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
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
%     bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
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
%     bioTree1{nodeInTree1(1)}.node{nodeInTree1(2)}.Out{nodeInTree1(3)}.traceInfo.measurment(sizeTrace1:sizeTrace1+sizeTrace2-1)=bioTree2{connectPoint}.root{j}.traceInfo.measurment(1:sizeTrace2);
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
    nS2(pixelIdxList2)=1;
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
        nS3(pixelIdxList3)=1;
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
function bioTree=bioTreeSizeReduction(bioTree,frameShift)
for iframe=frameShift+1:size(bioTree,2)
    if ~isempty(bioTree{iframe}.root)
        for iroot=1:size(bioTree{iframe}.root,2)
            if bioTree{iframe}.root{iroot}.is2Node==true
                bioTree{iframe}.root{iroot}.traceInfo.pixelIdxList(end)=[];
%                 bioTree{iframe}.root{iroot}.traceInfo.measurment(end)=[];
            end
        end
    end
    if ~isempty(bioTree{iframe}.node)
        for iNode=1:size(bioTree{iframe}.node,2)
            for iOut=1:size(bioTree{iframe}.node{iNode}.Out,2)
                if bioTree{iframe}.node{iNode}.Out{iOut}.is2Node==true
                    bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.pixelIdxList(end)=[];
%                     bioTree{iframe}.node{iNode}.Out{iOut}.traceInfo.measurment(end)=[];
                end
            end
        end
    end
end
end
function allTreeSave(bioTree)
mkdir('allTree');
saveFile=strcat('allTree\bioTeeStack');
bytesTemp=whos('bioTree');
bytes=bytesTemp.bytes;
treeNum=fix(bytes/2000000000)+1;
savePoint=findSavePoint(bioTree,treeNum);
startFrame=1;
for iTree=1:treeNum
    if treeNum== 1
        save(saveFile,'bioTree');
        return;
    end
    bioTreeStack=bioTree(startFrame:savePoint(iTree));
    saveFile1=strcat(saveFile,num2str(iTree));
    save(saveFile1,'bioTreeStack');
    startFrame=savePoint(iTree)+1;
    bioTreeStack=[];
end
end
function savePoint=findSavePoint(bioTree,treeNum)
savePoint=[];
if treeNum==1
    return;
end
startFrame=1;
for iframe=100:100:size(bioTree,2)
bioTreeTmep=bioTree(startFrame:iframe)  ;  
bytesTemp=whos('bioTreeTmep');
bytes=bytesTemp.bytes;
pointNum=fix(bytes/2000000000);
if pointNum==1
    savePoint=[savePoint,iframe-100];
    startFrame=iframe-99;
    if size(savePoint,2)==treeNum-1
        savePoint=[savePoint,size(bioTree,2)];
        return;
    end
end
end
end