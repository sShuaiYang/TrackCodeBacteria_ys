function [mergedImages] = showSegMergeOriginalImage(originalImages,maskImages)
%将原始荧图像与mask图像进行轮廓重叠
% Shuai Yang 2021.06.21
% Shuai Yang 2021.08.16 update suitable for single image or imageStack
% 旧函数 改用imshowMerge_ys 和 imshowColorlabel_ys 函数
% 2022/6/8
imsz = size(originalImages);%image size
masksz = size(maskImages);%mask size
if ~isequal (imsz,masksz)
    disp('originalImage and masks are unqual in size');
    return
end
if numel(imsz) == 2 % 只有一张图像
    cell_mask = logical(maskImages);
    cell_mask = imfill( cell_mask,'holes');
    outline = bwmorph(cell_mask,'remove');
    I0 = rescale(double(originalImages))*255;
    I0 = uint8(I0);
    mergedImage = imoverlay(I0,outline,[1 0 0]);
    figure,imshow(mergedImage);
    mergedImages = mergedImage;
else
    imageType = 'uint8';
    mergedImages = zeros(imsz(1),imsz(2),3,imsz(3),imageType);
    for iImage = 1:size(originalImages,3)
        cell_mask = logical(maskImages(:,:,iImage));
        cell_mask = imfill( cell_mask,'holes');
        outline = bwmorph(cell_mask,'remove');
        I0 = rescale(double(originalImages(:,:,iImage)))*255;
        I0 = uint8(I0);
        mergedImage = imoverlay(I0,outline,[1 0 0]);
        % Save merged images using 8bitRGB Color，for movie demonstration
        %         dirSave = '';
        %         imwrite(mergedImage,[dirSave,filesep,num2str(iImage,'%03d'),'.tif']);
        mergedImages(:,:,:,iImage) = mergedImage;
    end
end
end