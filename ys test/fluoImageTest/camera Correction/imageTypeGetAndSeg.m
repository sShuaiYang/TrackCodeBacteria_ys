%% 图像类型鉴定 %'PhC','Fluo','Binary'
function [imageMask,imType] = imageTypeGetAndSeg(myImage)
% 鉴定图像类型 并进行分割 获得mask
% myImage: single image
%'PhC','Fluo','Binary'
% Shuai Yang 2022/6/15

if numel(unique(double(myImage))) <= 2
    imType = 'Binary';
    imageMask = logical(myImage);
    return
end

imsz = size(myImage);
I1 = myImage;
I1 = I1(round(imsz(1)/2-imsz(1)/4):round(imsz(1)/2+imsz(1)/4), ...
    round(imsz(2)/2-imsz(2)/4):round(imsz(2)/2+imsz(2)/4));
I1 = double(I1);

% 通过背景判断 认为荧光图像的背景小于200
if min(I1(:)) < 200
    imType = 'Fluo';
    [imageMask] = fluoImCellMask_ys(myImage);
else
    imType = 'PhC';
    [imageMask,~] = phCImProcessing_supperSegger(myImage);
end

end