function [afterProcessingImages]=fluoImageProcessing_ystest1(beforeProcessingImages,thresh) 
%JZY旧版本的明场图像处理程序修改得到荧光图像处理code
%imageType='uint16'; %荧光图像拍摄的一般为16bit 后面转化为8bit处理

% %细菌较多时 直接减100 
% beforeProcessingImages = beforeProcessingImages-100;
 % 适合初始细菌较少时 提取背景  减去背景
 BG = imopen(beforeProcessingImages(:,:,1), strel('disk', 70));
parfor iframe=1:size(beforeProcessingImages,3)
    I0 = beforeProcessingImages(:,:,iframe);
    beforeProcessingImages(:,:,iframe) = imsubtract(I0, BG);
end

 % 将图像转化为8bit
 imageType='uint8';
 convertImageStack = zeros(size(beforeProcessingImages),imageType);
 parfor iframe = 1:size(beforeProcessingImages,3)
     I0 = beforeProcessingImages(:,:,iframe);
     convertImageStack(:,:,iframe) = uint8(double(I0)/double(max(I0(:)))*255);
 end
beforeProcessingImages = convertImageStack;
% gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter

% here you can change the process parameters
if nargin<2 || isempty(thresh)
    I = beforeProcessingImages (:,:,1);
    I = imgaussfilt(I, 2);% sigma = 2 Gaussian Blur Filter
    I = imfilter(I, edgeFilter);
    I = imgaussfilt(I, 2);
    T = graythresh(I);
    thresh = T *255;
end
areaThreshold=60;  % proper 60, overlap 60+40
imageProcessingInfo.grayThresh=thresh;
imageProcessingInfo.areaThreshold=areaThreshold;

afterProcessingImages = zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
parfor iframe=1:size(beforeProcessingImages,3)
    
    afterProcessingImages(:,:,iframe)=imgaussfilt(beforeProcessingImages(:,:,iframe),2);
    %     afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter);
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    %     afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imgaussfilt(afterProcessingImages(:,:,iframe),2);
    afterProcessingImages(:,:,iframe)=imbinarize(afterProcessingImages(:,:,iframe),thresh/(2^8-1));
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes');
    
end

afterProcessingImages=logical(afterProcessingImages);
end