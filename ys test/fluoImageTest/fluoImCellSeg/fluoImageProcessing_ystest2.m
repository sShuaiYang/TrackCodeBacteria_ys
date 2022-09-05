function [afterProcessingImages] = fluoImageProcessing_ystest2(beforeProcessingImages,thresh)
%在fluoImageProcessing_ystest程序中，如果一个图像中有极个别荧光特别强的点
%~20倍[120 vs 2400] 转化为8bit图像 mask生成会有问题
%此 code 直接对16bit荧光图像进行处理 不转化， 阈值thresh需要自己设定
% Shuai Yang 2020.08.07
% beforeProcessingImages 16 bit fluo image
imageType = 'uint16';
% set threshold
if nargin<2 || isempty(thresh)
    thresh = 20;% 默认20
    thresh = thresh/(2^16-1);
else
    thresh = thresh/(2^16-1);
end
% here intiatlize the maskImages stacks
afterProcessingImages = zeros(size(beforeProcessingImages),imageType);
%建立一个边框 用于删除边缘识别错误的mask
edgeBorder = true(size(beforeProcessingImages(:,:,1)));
edgeBorder(1:7,:) = false;
edgeBorder(end-6:end,:) = false;
edgeBorder(:,1:7) = false;
edgeBorder(:,end-6:end) = false;
% Filter set
gaussianFilter = fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
edgeFilter = (ones(5,5)).*-1;edgeFilter(3,3) = 24;  %here create edageFilter
areaThreshold = 80;

parfor iframe = 1:size(beforeProcessingImages,3) 
    afterProcessingImages(:,:,iframe) = imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter);
    afterProcessingImages(:,:,iframe) = imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe) = imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe) = imbinarize(afterProcessingImages(:,:,iframe),thresh);
    afterProcessingImages(:,:,iframe) =  afterProcessingImages(:,:,iframe) & edgeBorder;
    afterProcessingImages(:,:,iframe) = bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    afterProcessingImages(:,:,iframe) = imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe) = imfill(afterProcessingImages(:,:,iframe),'holes');
    afterProcessingImages(:,:,iframe) = bwmorph(afterProcessingImages(:,:,iframe),'open');
end

afterProcessingImages = logical(afterProcessingImages);
end