iLeaf=8;
image=false(bioTree{1}.imageSize);
image(bioTree{end}.leavies{iLeaf}.leaviesPixelDetail)=1;
image=imfill(image,'holes');
pixelIdxList=[];

image1=MIJ.getCurrentImage;
image1=uint8(image1);
image1=im2bw(image1,0.1);
image1=bwmorph(image1,'remove');
pixelIdxList{1}=find(image1==1);

nodeInfo=bioTree{end}.leavies{iLeaf}.nodeInfo;
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.is2Node=1;
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.leafInfo=[];
newNode=[size(bioTree,2),size(bioTree{end}.node,2)+1];
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.nodeInfo=[newNode,1];
bioTree{nodeInfo(1)}.node{nodeInfo(2)}.Out{nodeInfo(3)}.traceInfo.pixelIdxList(end)=[];
bioTree{newNode(1)}.node{newNode(2)}.In{1}.nodeInfo=nodeInfo;
bioTree{newNode(1)}.node{newNode(2)}.In{1}.isNode=1;
bioTree{newNode(1)}.node{newNode(2)}.In{1}.leafInfo=[];
for i=1:size(pixelIdxList,2)
    if i==1
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.is2Node=0;
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.leafInfo=[size(bioTree,2),iLeaf];
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.nodeInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.Out{1}.traceInfo.pixelIdxList{1}=pixelIdxList{i};
        bioTree{end}.leavies{iLeaf}.is2Node=1;
        bioTree{end}.leavies{iLeaf}.nodeInfo=[newNode,1];
        bioTree{end}.leavies{iLeaf}.rootInfo=[];
        bioTree{end}.leavies{iLeaf}.leaviesPixelDetail=pixelIdxList{i};
    else
        leafNum=size(bioTree{end}.leavies,2)+1;
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.is2Node=0;
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.leafInfo=[size(bioTree,2),leafNum];
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.nodeInfo=[];
        bioTree{newNode(1)}.node{newNode(2)}.Out{i}.traceInfo.pixelIdxList{1}=pixelIdxList{i};
        bioTree{end}.leavies{leafNum}.is2Node=1;
        bioTree{end}.leavies{leafNum}.nodeInfo=[newNode,i];
        bioTree{end}.leavies{leafNum}.rootInfo=[];
        bioTree{end}.leavies{leafNum}.leaviesPixelDetail=pixelIdxList{i};
    end
end