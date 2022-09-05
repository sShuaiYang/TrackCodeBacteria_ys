function xyt_stack=import_tiff_stack(fname)

% returns an array with xy in first two dimensions, stack in thrid
% dimension, and multi [color] channels (if present) in fourth dimension
% NB - assumes all images in stack have same size, bit depth, and color ch.

info =imfinfo(fname);

[img_size1 img_size2]=size(imread(fname, 1, 'Info', info));
img_type=['uint' num2str(info(1).BitDepth)];

num_images = numel(info);

xyt_stack=zeros([img_size1 img_size2 num_images], img_type);
%preallocate array of correct type and size to improve speed

for k = 1:num_images
    A = imread(fname, k, 'Info', info);
    xyt_stack(:,:,k)=A;
end
