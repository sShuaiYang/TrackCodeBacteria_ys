function bestPosition=stepCorrectionForImageStack(imageStack)
% 计算出imageStack的偏移量，并对imageStack进行校正
bestPosition(1,:)=[size(imageStack,1),size(imageStack,2)];
parfor i=1:size(imageStack,3)-1
    correlationMatrix=caculateCrossCorrelationForImage(imageStack(:,:,i),imageStack(:,:,i+1),10);
%     correlationMatrix=xcorr2(imageStack(:,:,i)-100,imageStack(:,:,i+1)-100);
    [a,b]=find(correlationMatrix==max(max(correlationMatrix)));
    bestPosition(i+1,:)=[a(1),b(1)];
end
end
function  image=imageStackCorrectionWithBestPosition(image,bestPosition)
deltaPosition=bestPosition-512;
for i=1:size(deltaPosition,1)-1
    deltaPosition(i+1,:)=deltaPosition(i+1,:)+deltaPosition(i,:);
end
bestPosition1=deltaPosition+512;
fullSize=[size(imageStack,1)+size(imageStack,1)-1,size(imageStack,2)+size(imageStack,2)-1];
centroid=floor((fullSize+1)/2);
delta1=floor((size(imageStack,1)+1)/2)-1;
delta2=size(imageStack,1)-floor((size(imageStack,1)+1)/2);
parfor i=2:size(imageStack,3)
    fullImage=uint16(ones(fullSize)*backGround);
    fullImage(bestPosition1(i,1)-delta1:bestPosition1(i,1)+delta2,bestPosition1(i,2)-delta1:bestPosition1(i,2)+delta2)=imageStack(:,:,i);
    imageStack1(:,:,i)=fullImage(centroid(1)-delta1:centroid(1)+delta2,centroid(2)-delta1:centroid(2)+delta2);
end
end