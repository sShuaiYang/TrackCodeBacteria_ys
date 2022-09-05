function result=twoDimGaussian(a,b)
% ¶þÎ¬¸ßË¹blur
xMin=1000;
yMin=1000;
xMax=30000;
yMax=30000;
sigma1=200;
interval1=100;
sigma2=200;
interval2=100;
rowPos=xMin:interval1:xMax;
rowPos=rowPos';
columnPos=yMin:interval2:yMax;
resultImage=zeros(numel(rowPos),numel(columnPos));
for i=1:numel(a)
    disp(i)
    rowGauss=1/(2*pi)^0.5/sigma1*exp(-(rowPos-a(i)).^2/2/sigma1^2);
    columnGauss=1/(2*pi)^0.5/sigma2*exp(-(columnPos-b(i)).^2/2/sigma2^2);
    resultImage=resultImage+rowGauss*columnGauss;
end
result.Image=resultImage;
result.numPoint=numel(a);
result.x=rowPos';
result.y=columnPos;
end