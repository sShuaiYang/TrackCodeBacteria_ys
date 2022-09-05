% image1=imageStack;
% image2=imageStack;
% image3=imageStack;
% for i=1:27
% nodeList=newList(i,1:4);
% pixelIdxList=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList{nodeList(4)};
% image1(pixelIdxList)=255;image2(pixelIdxList)=0;image3(pixelIdxList)=0;
% end
% image=cat(3,image1,image2,image3);
% imshow(image)
% 

image1=imageStack;
image2=imageStack;
image3=imageStack;
bacteriaList=bacteriaFrameInfo{5695}.bacteriaInfo;
bacteriaBranch=bacteriaList(:,5);
for i=1:size(bacteriaBranch,1)
    if ismember(bacteriaBranch(i),aimBranch)
        if bacteriaList(i,3)==0
            rootInfo=bacteriaList(i,[1:2,4]);
            pixelIdxList=bioTree{rootInfo(1)}.root{rootInfo(2)}.traceInfo.pixelIdxList{rootInfo(3)};
        else
            nodeList=bacteriaList(i,1:4);
            pixelIdxList=bioTree{nodeList(1)}.node{nodeList(2)}.Out{nodeList(3)}.traceInfo.pixelIdxList{nodeList(4)};
        end
        image1(pixelIdxList)=255;image2(pixelIdxList)=0;image3(pixelIdxList)=0;
    end
end
image=cat(3,image1,image2,image3);
imshow(image)