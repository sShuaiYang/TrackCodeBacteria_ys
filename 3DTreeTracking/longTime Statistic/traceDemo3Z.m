function traceDemo3Z(B1,B2,A,imageSize)
figure;hold on
for i=1
    xyIdx1=idx2XyIdx(B1{i},imageSize);
    xyIdx2=idx2XyIdx(B2{i},imageSize);
    xyIdx=[xyIdx1,xyIdx2];
    plot3(xyIdx(1,:),xyIdx(2,:),5*i*ones(size(xyIdx(1,:))),'lineStyle','none','Marker','.');hold on
end
for i=2
    xyIdx=idx2XyIdx(A{i},imageSize);
    plot3(xyIdx(1,:),xyIdx(2,:),5*(i+5)*ones(size(xyIdx(1,:))),'lineStyle','none','Marker','.','MarkerEdgeColor',[1,0,0]);hold on
end
end
function xyIdx=idx2XyIdx(pixelIdxList,pictureSize)
% this function is used to convert linearIdx to a small BWImage
xSize=pictureSize(1);
yresult=ceil(pixelIdxList/xSize);
xresult=pixelIdxList-(yresult-1)*xSize;
xyIdx=cat(1,xresult',yresult');
end