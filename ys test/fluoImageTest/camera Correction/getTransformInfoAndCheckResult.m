function transformInfo = getTransformInfoAndCheckResult(selectedMovingPoints,selectedFixedPoints,movingImage,fixedImage)
% 根据选择的点 获取两张图像的transformInfo 并展示校正后的图像
% Shuai Yang
% 2022/6/24

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'projective');
Rfixed = imref2d(size(movingImage));%fixed image size
mIregistered = imwarp(movingImage,tform,'OutputView',Rfixed);%moving image registered

transformInfo = tform;
I0 = fixedImage - 100;
im_PixelValues = sort(double(I0(:)),'descend'); 
I0 = uint8((rescale(double(I0),'InputMax',im_PixelValues(1000)))*255);
fixedImage = I0;

I0 = mIregistered - 100;
im_PixelValues = sort(double(I0(:)),'descend'); 
I0 = uint8((rescale(double(I0),'InputMax',im_PixelValues(1000)))*255);
mIregistered = I0;
% C = imfuse(fixedImage,mIregistered,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
C = cat(3,fixedImage,mIregistered,fixedImage*0);
figure, imshow(C,[])
end
