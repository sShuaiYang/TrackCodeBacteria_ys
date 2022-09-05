function [peocessedImages,imageProcessingInfo]=myImageProcessing_cBFys(cBFImages,picType) %this function can segreate orignal images as you want and returen the mask images
% modified from myImageProcessing.m &myImageProcessing_cBF.m ,针对相关明场图像高对比对的特性，去掉了
%高斯滤波和边缘滤波，对图像直接进行二值化，然后根据参数判定进行bw图像的处理
%加入了边缘检测对图像生成mask 2020.01.01ys
imageType='uint8'; %here you can change your image type

cropInfo=true(size(cBFImages,1),size(cBFImages,2));
peocessedImages=zeros(size(cBFImages),imageType);% here intiatlize the maskImages stacks

if strcmp(picType,'w')
    grayThresh=100;
    areaThreshold=100;  % proper 60, overlap 60+40
    maxIntensity=45;
end
if strcmp(picType,'b')
    %     grayThresh=90;
    %     areaThreshold=10;
    %     maxIntensity=170;     %%% 100X black
    grayThresh=90;
    areaThreshold=100;
    maxIntensity=150;
end

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
cropArea=~cropInfo;

cBFImages = cleanupImageUsingReconstuction(cBFImages);

parfor iframe=1:size(cBFImages,3)    
    peocessedImages(:,:,iframe) = imbinarize(cBFImages(:,:,iframe),grayThresh/255);
    peocessedImages(:,:,iframe) = bwareaopen(peocessedImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open');
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'bridge'); % bridge process
    % edge detection
    tempBW = edge(cBFImages(:,:,iframe),'log');
    %     tempBW=bwmorph( tempBW,'fill');
    tempBW = imfill( tempBW,'holes');
    se = strel('disk',2);
    %     tempBW = imopen( tempBW, ones(2,2));
    tempBW = imopen( tempBW, se);
    tempBW = bwareaopen( tempBW, areaThreshold);
    tempBW = peocessedImages(:,:,iframe)| tempBW ;
    peocessedImages(:,:,iframe) = imopen(peocessedImages(:,:,iframe), se);
    peocessedImages(:,:,iframe) = tempBW;
    image=peocessedImages(:,:,iframe);
    cc=regionprops(logical(peocessedImages(:,:,iframe)),cBFImages(:,:,iframe),'MaxIntensity','PixelIdxList','MeanIntensity');
    for iCC = 1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    cc=regionprops(logical(1-peocessedImages(:,:,iframe)),cBFImages(:,:,iframe),'PixelIdxList','MeanIntensity','FilledArea');
    for iCC=1:size(cc,1)
        if cc(iCC).MeanIntensity==255 && cc(iCC).FilledArea<=areaThreshold
            image(cc(iCC).PixelIdxList)=1;
        end
    end
    peocessedImages(:,:,iframe)=image;
    peocessedImages(:,:,iframe)=imclearborder(peocessedImages(:,:,iframe) | cropArea);
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % fill holes process
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
end
peocessedImages=logical(peocessedImages);

parfor i=1:size(peocessedImages,3)
    peocessedImages(:,:,i)=imclearborder(peocessedImages(:,:,i));
    cc=regionprops(logical(peocessedImages(:,:,i)),cBFImages(:,:,i),'MaxIntensity','FilledArea','PixelIdxList');
    image=peocessedImages(:,:,i);
    for iCC=1:size(cc,1)
        %         if (cc(iCC).FilledArea<=areaThreshold+40 && cc(iCC).MaxIntensity<=maxIntensity) || cc(iCC).FilledArea<=areaThreshold || cc(iCC).MaxIntensity<=120
        if cc(iCC).FilledArea<=areaThreshold || cc(iCC).MaxIntensity<=maxIntensity
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    peocessedImages(:,:,i)=image;
    %     afterProcessingImages(:,:,i)=bwmorph(afterProcessingImages(:,:,i),'remove'); % fill holes process
    peocessedImages(:,:,i)=bwmorph(peocessedImages(:,:,i),'fill');%Fills isolated interior pixels (individual 0s that are surrounded by 1s).
    peocessedImages(:,:,i)=imfill(peocessedImages(:,:,i),'holes'); % fill holes
    %         afterProcessingImages(:,:,i)=imerode(afterProcessingImages(:,:,i),ones(2));%ys 明场相关图像改为ones(2),否则腐蚀太过
    %     afterProcessingImages(:,:,i)=bwmorph(afterProcessingImages(:,:,i),'open');
    
    peocessedImages(:,:,i)=imopen(peocessedImages(:,:,i),strel('disk',2));
    peocessedImages(:,:,i)=imopen(peocessedImages(:,:,i),ones(3));
    peocessedImages(:,:,i)=bwmorph(peocessedImages(:,:,i),'open');
end

%%
% 先腐蚀再骨架处理再重组细菌进行分割
parfor i=1:size(peocessedImages,3)
    image=peocessedImages(:,:,i);
    stats = regionprops(peocessedImages(:,:,i),'PixelIdxList');
    pixelNum=zeros(2,numel(stats));
    for j=1:numel(stats)
        pixelNum(1,j)=j;
        pixelNum(2,j)=numel(stats(j).PixelIdxList);
    end
    templogic=(pixelNum(2,:)>=600);
    if sum(templogic)>0
        stats=stats(templogic);
        for k=1:numel(stats)
            image(stats(k).PixelIdxList)=false;
            tempImage=false(size(image));
            tempImage(stats(k).PixelIdxList)=true;% 选取连在一起的细菌进行分割 通过面积筛选            
            tempImage=imerode(tempImage,ones(4));%腐蚀进行分割
            tempImage = bwareaopen(tempImage, 5);%删除小的零星的pixel
            skel = bwmorph(tempImage,'skel',Inf);%获取骨架
            E = bwmorph(skel, 'endpoints');
            skel=skel&(~E); %删除endpoint 要不膨胀会连在一起
            E = bwmorph(skel, 'endpoints');
            skel=skel&(~E); %删除两次
            %重组细菌
            se = strel('disk',4);
            bwSegNew=imdilate(skel,se);
            image=image|bwSegNew;
            peocessedImages(:,:,i)=image;
        end
    end
end

%% 试图通过‘open’操作分割连在一起的细菌
% se = strel('disk',2);
% rptTime=2;
% for n=1:rptTime
%     parfor i=1:size(afterProcessingImages,3)
%         image=afterProcessingImages(:,:,i);
%         stats = regionprops(afterProcessingImages(:,:,i),'PixelIdxList');
%         pixelNum=zeros(2,numel(stats));
%         for j=1:numel(stats)
%             pixelNum(1,j)=j;
%             pixelNum(2,j)=numel(stats(j).PixelIdxList);
%         end
%         templogic=(pixelNum(2,:)>=600);
%         if sum(templogic)>0
%             stats=stats(templogic);
%             for k=1:numel(stats)
%                 image(stats(k).PixelIdxList)=0;
%                 tempImage=false(size(image));
%                 tempImage(stats(k).PixelIdxList)=true;
%                 %             tempImage=logical(tempImage);
%                 %             tempImage=imopen(tempImage,se);
%                 tempImage=bwmorph(tempImage,'open');
%                 image=image|tempImage;
%                 afterProcessingImages(:,:,i)=image;
%             end
%         end
%     end
% end

end
%% 对连在一起的细菌通过imopen 进行分割
% function imageNew=splitBigRegions(image)
% % se = strel('disk',2);
% stats = regionprops(image,'PixelIdxList');
% pixelNum=zeros(2,numel(stats));
% for j=1:numel(stats)
%     pixelNum(1,j)=j;
%     pixelNum(2,j)=numel(stats(j).PixelIdxList);
% end
% templogic=(pixelNum(2,:)>=600);
% if sum(templogic)>0
%     stats=stats(templogic);
%     for i=1:numel(stats)
%         image(stats(i).PixelIdxList)=0;
%         localImage=false(size(image,1),size(image,2));
%         localImage(stats(i).PixelIdxList)=true;
%         %             tempImage=logical(tempImage);
% %         localImage=imopen(localImage,se);
%         localImage=bwmorph(localImage,'open');
%         localImage_new=splitBigRegions(localImage);
%         imageNew=image|localImage_new;
%     end
% else
%     imageNew=image;
% end
%
% end

