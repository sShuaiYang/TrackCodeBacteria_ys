function [outim1,outim2] = imshowMerge_ys(bwmask,images)
% 用于将mask图像和荧光/明场图像进行边缘和区域重叠
% 将mask的edge 和 region merge到原始图像上
% bwmask mask
% I 为原始图像 作为mask的背景 荧光或者明场 fluoIm or BF
% outim1: edge overlay; outim2: region overlay;
% Shuai Yang 2022/02/23

images = images-100;
mergeColor=[0,1,0];
imsz = size(images);
bw = false(size(bwmask));
I = uint8 (zeros(size(images)));
imNum = size(images,3);

for iIm = 1:imNum
    I0 = rescale(double(images(:,:,iIm)))*255;
    I(:,:,iIm) = uint8(I0);
    bw0 = logical(double(bwmask(:,:,iIm)));
    %     bw(:,:,iIm) = edge(bw0,'log',0);
    bw(:,:,iIm) = bwmorph(bw0,'remove');

end

outim1 = uint8(zeros(imsz(1),imsz(2),3,imNum));
outim2 = uint8(zeros(imsz(1),imsz(2),3,imNum));
for iIm = 1:imNum
    outim1(:,:,:,iIm) = imoverlay(I(:,:,iIm),bw(:,:,imNum),mergeColor);
    outim2(:,:,:,iIm) = labeloverlay(I(:,:,iIm),bwmask(:,:,imNum));
end

if imNum ==1
    figure,imshow(outim1)
    figure,imshow(outim2)
else
    implay(outim1)
end

end