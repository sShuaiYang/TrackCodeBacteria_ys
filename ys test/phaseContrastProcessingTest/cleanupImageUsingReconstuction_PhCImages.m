function [IobrcbrStack] = cleanupImageUsingReconstuction_PhCImages(image)
% Shuai Yang 利用图像重构来cleanup 相差图像 ys 2020.04.17

imageType='uint8';
IobrcbrStack = zeros(size(image),imageType);
% tic;
parfor i=1:size(image,3)

    I = image(:,:,i);
    se = strel('disk',5);
    Io = imopen(I,se);
    % imshow(Io)
    % title('Opening Io')
    
    Ie = imerode(I,strel('disk',5));
    Iobr = imreconstruct(Ie,I);
    % imshow(Iobr)
    % title('Opening-by-Reconstruction Iobr')
    
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    IobrcbrStack(:,:,i) = imcomplement(Iobrcbr);
end
% toc;
end

