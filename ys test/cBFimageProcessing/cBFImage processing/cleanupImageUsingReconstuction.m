function [Iobr] = cleanupImageUsingReconstuction(image)

% image=import_tiff_stack('E:\2019-12-26 PAO1_deltpslpelfliC_IP32_100x ys\field0004\BF\imageBF04000.tif');
imageType='uint8';
Iobr = zeros(size(image),imageType);
% tic;
parfor i=1:size(image,3)
    I=image(:,:,i);
    se = strel('disk',2);
    Io = imopen(I,se);
    % imshow(I), title('Original Image')
    % imshow(Io),title('Opening I-Io')
    
    Ie = imerode(Io,strel('disk',3));
    % imshow(Ie), title('Eroded image I-Ie ')
    Iobr(:,:,i) = imreconstruct(Ie,Io);
    % imshow(Iobr),title('Opening-by-Reconstruction I-Ie-Iobr ')
end
% toc;
end

