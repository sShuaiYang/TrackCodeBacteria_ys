function [resultImage,x,y]=twoDimGaussian(positionInfo,sigma)
xMin=min(min(positionInfo(:,1)));
yMin=min(min(positionInfo(:,2)));
xMax=max(max(positionInfo(:,1)));
yMax=max(max(positionInfo(:,2)));
interval=0.01;
rowPos=xMin-2*sigma:interval:xMax+2*sigma;
rowPos=rowPos';
columnPos=yMin-2*sigma:interval:yMax+2*sigma;
resultImage=zeros(numel(rowPos),numel(columnPos));
for i=1:size(positionInfo,1)
    disp(i)
    rowGauss=1/(2*pi)^0.5/sigma*exp(-(rowPos-positionInfo(i,1)).^2/2/sigma^2);
    columnGauss=1/(2*pi)^0.5/sigma*exp(-(columnPos-positionInfo(i,2)).^2/2/sigma^2);
    resultImage=resultImage+rowGauss*columnGauss;
end
resultImage=resultImage/size(positionInfo,1);
x=rowPos';
y=columnPos;
end
