function meragedMask=createMeragedMask(mask1,mask2)
gaussianFilter=fspecial('gaussian', [3,3], 5);
parfor iframe=1:size(mask1,3)
    mask1(:,:,iframe)= mask1(:,:,iframe)+ mask2(:,:,iframe);
    mask1(:,:,iframe)=imfilter(mask1(:,:,iframe),gaussianFilter);
end
meragedMask=mask1;
end