function clusterMaskImages=clustingImageProcessing(beforeProcessingImages)
imageType='uint8';
gaussianFilter1=fspecial('gaussian', [7,7], 5);
gaussianFilter2=fspecial('gaussian', [3,3], 5);
edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;
clusterAreaThreshold1=1500;
clusterAreaThreshold2=800;
clusterMaskImages=zeros(size(beforeProcessingImages),imageType);
% clusterMaskImages=false(size(beforeProcessingImages));
parfor iframe=1:size(beforeProcessingImages,3)
    %
    clusterMaskImages(:,:,iframe)=imadjust(beforeProcessingImages(:,:,iframe));
    %
    %     clusterMaskImages(:,:,iframe)=imfilter(clusterMaskImages(:,:,iframe),gaussianFilter);
    %
    %     clusterMaskImages(:,:,iframe)=imfilter(clusterMaskImages(:,:,iframe),edgeFilter);
    %
    %     clusterMaskImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter);
    
    [clusterMaskImages(:,:,iframe),~]=clustingSegration(clusterMaskImages(:,:,iframe),0,1,1);
    
    %
    %     figure;imshow(clusterMaskImages(:,:,iframe));
    
    %
    clusterMaskImages(:,:,iframe)=im2bw_ent(clusterMaskImages(:,:,iframe));
    clusterMaskImages(:,:,iframe)=~clusterMaskImages(:,:,iframe);
    %
%     clusterMaskImages(:,:,iframe)=bwareaopen(clusterMaskImages(:,:,iframe),clusterAreaThreshold1); %remove the objecvies which area is smaller than threshold
    
%     clusterMaskImages(:,:,iframe)=bwmorph(clusterMaskImages(:,:,iframe),'dilate');
    
%     clusterMaskImages(:,:,iframe)=imfilter(clusterMaskImages(:,:,iframe),gaussianFilter1);
    %
    %     clusterMaskImages(:,:,iframe)=bwmorph(clusterMaskImages(:,:,iframe),'dilate');
    clusterMaskImages(:,:,iframe)=imfill(clusterMaskImages(:,:,iframe));
    clusterMaskImages(:,:,iframe)=bwmorph(clusterMaskImages(:,:,iframe),'erode');
     clusterMaskImages(:,:,iframe)=bwmorph(clusterMaskImages(:,:,iframe),'erode');
%     
    clusterMaskImages(:,:,iframe)=imfilter(clusterMaskImages(:,:,iframe),gaussianFilter1);
    
        clusterMaskImages(:,:,iframe)=bwareaopen(clusterMaskImages(:,:,iframe),clusterAreaThreshold2); %remove the objecvies which area is smaller than threshold
    
    
end
clusterMaskImages=logical(clusterMaskImages);
end


