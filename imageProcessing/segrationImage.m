function  afterProcessingImages=segrationImage(beforeProcessingImages,maskImages) % this function can overlay your mask and orignal images 
maskColor=[0,1,0];
imageType='uint8'; %here you can chenge your image type
afterProcessingImages=zeros(size(beforeProcessingImages,1),size(beforeProcessingImages,2),3,size(beforeProcessingImages,3),imageType);
parfor iframe=1:size(beforeProcessingImages,3)
%     afterProcessingImages(:,:,iframe)=immultiply((~bwperim(maskImages(:,:,iframe),4)),beforeProcessingImages(:,:,iframe));
%      afterProcessingImages(:,:,iframe)=immultiply(~maskImages(:,:,iframe),beforeProcessingImages(:,:,iframe));
afterProcessingImages(:,:,:,iframe)=imoverlay(beforeProcessingImages(:,:,iframe),maskImages(:,:,iframe),maskColor);
end
end