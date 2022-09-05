function finalImage=getIndImage(image)
cc=bwconncomp(image);
finalImage=false(size(image));
for i=1:cc.NumObjects
    image1=false(size(image));
    image1(cc.PixelIdxList{i})=1;
    image2=bwmorph(image1,'endpoints');
    [xPix,yPix]=find(image2==1);
    xyPix=[xPix(1),yPix(1)];
    B=bwtraceboundary(image1,xyPix,'E');
    B(2:2:end,:)=[];
    for j=1:size(B,1)
        finalImage(B(j,1),B(j,2))=1;
    end
end
end