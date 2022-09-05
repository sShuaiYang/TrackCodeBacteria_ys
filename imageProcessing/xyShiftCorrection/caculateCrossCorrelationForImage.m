function correlationMatrix=caculateCrossCorrelationForImage(image1,image2,step,backGround)
% for calculater the cross correlation of two image
% step is the searching range
% backGround should be calculater or set by user
if nargin==3
    imageIndex=image2(:);
    imageIndex=sort(imageIndex);
    backGround=mean(imageIndex(1:end/3));
end
fullSize=[size(image1,1)+size(image2,1)-1,size(image1,2)+size(image2,2)-1];
correlationMatrix=zeros(fullSize);
centroid=floor((fullSize+1)/2);
matrixSequence=(centroid(1)-step:centroid(2)+step)';
[x,y]=meshgrid(matrixSequence,matrixSequence);
x=x(:);
y=y(:);
delta1=floor((size(image2,1)+1)/2)-1;
delta2=size(image2,1)-floor((size(image2,1)+1)/2);
for i=1:numel(x)
    fullImage=double(ones(fullSize)*backGround);
    fullImage(x(i)-delta1:x(i)+delta2,y(i)-delta1:y(i)+delta2)=image2;
    image2New=fullImage(centroid(1)-delta1:centroid(1)+delta2,centroid(2)-delta1:centroid(2)+delta2);
    correlationMatrix(x(i),y(i))=sum(sum(double(image1).*image2New));
%     correlationMatrix(i)=sum(sum(double(image1).*image2New));
end
end


