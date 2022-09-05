%% past code
function [afterProcessingImages,imageProcessingInfo]=myImageProcessing(beforeProcessingImages,picType) %this function can segreate orignal images as you want and returen the mask images
imageType='uint16'; %here you can chenge your image type
% % here are three types of crop method
%
% % % 1. x or y crop
% xCrop=270:2160;
% cropInfo=false(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
% cropInfo(xCrop,:)=1;
% beforeProcessingImages=beforeProcessingImages(xCrop,:,:);
%
% % 2.x and y crop
% xCrop=1000:1200;
% yCrop=700:900;
% cropInfo=false(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
% cropInfo(xCrop,yCrop)=1;
% beforeProcessingImages=beforeProcessingImages(xCrop,yCrop,:);
%
% 3. x and y remove
% xCrop=1:270;
% yCrop=1:2560;
% cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
% cropInfo(xCrop,yCrop)=0;
% beforeProcessingImages(xCrop,yCrop,:)=90;

% 4. no need to crop
cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));

gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
% edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
if strcmp(picType,'w')
    %     grayThresh=90;    
    %     areaThreshold=100;  % proper 60, overlap 60+40
    %     maxIntensity=230;      %%%%%%% long time 100X neo
    % orginal ys20191202
    %     grayThresh=65;
    %     areaThreshold=100;  % proper 60, overlap 60+40
    %   maxIntensity=180;      %%%%%%% agar 1
    %   correlated BF test ys 20191202
    grayThresh=40;
    areaThreshold=10;  % proper 60, overlap 60+40
    maxIntensity=10;
    %     grayThresh=90;
    %     areaThreshold=120;  % proper 60, overlap 60+40
    %     maxIntensity=240;      %%%%%%% agar 2
    %     grayThresh=70;
    %     areaThreshold=60;
    %     maxIntensity=180;     % general
end
if strcmp(picType,'b')
    %     grayThresh=90;
    %     areaThreshold=10;
    %     maxIntensity=170;     %%% 100X black
    grayThresh=50;
    areaThreshold=10;
    maxIntensity=100;
end

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
cropArea=~cropInfo;

parfor iframe=1:size(beforeProcessingImages,3)
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    %     thresholdImage=graythresh(afterProcessingImages(:,:,iframe)); %Otu thresholding images
    %     afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
    autoGrayThresh=graythresh(afterProcessingImages(:,:,iframe)); %Auto thresholding images
    afterProcessingImages(:,:,iframe)=imbinarize(afterProcessingImages(:,:,iframe),autoGrayThresh/1.2);
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open');
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'bridge'); % bridge process
    image=afterProcessingImages(:,:,iframe);
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList','MeanIntensity');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    cc=regionprops(logical(1-afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'PixelIdxList','MeanIntensity','FilledArea');
    for iCC=1:size(cc,1)
        if cc(iCC).MeanIntensity==255 && cc(iCC).FilledArea<=areaThreshold
            image(cc(iCC).PixelIdxList)=1;
        end
    end
    afterProcessingImages(:,:,iframe)=image;
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe) | cropArea);
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'dilate'); % dilate process
    %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % fill holes process
    %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
end
afterProcessingImages=logical(afterProcessingImages);

% % 这一步本来是想保留那些本不应分开的细菌之间的连接，但是也会带来很多错误，权衡后不要使用
% for i=2:size(afterProcessingImages,3)-1
%     currentImage=afterProcessingImages(:,:,i);
%     currentImage(afterProcessingImages(:,:,i-1)& afterProcessingImages(:,:,i+1))=1;
%     afterProcessingImages(:,:,i)=currentImage;
% end

parfor i=1:size(afterProcessingImages,3)
    afterProcessingImages(:,:,i)=imclearborder(afterProcessingImages(:,:,i));
    cc=regionprops(logical(afterProcessingImages(:,:,i)),beforeProcessingImages(:,:,i),'MaxIntensity','FilledArea','PixelIdxList');
    image=afterProcessingImages(:,:,i);
    for iCC=1:size(cc,1)
        %         if (cc(iCC).FilledArea<=areaThreshold+40 && cc(iCC).MaxIntensity<=maxIntensity) || cc(iCC).FilledArea<=areaThreshold || cc(iCC).MaxIntensity<=120
        if cc(iCC).FilledArea<=areaThreshold || cc(iCC).MaxIntensity<=120
            image(cc(iCC).PixelIdxList)=0;
        end
    end
    afterProcessingImages(:,:,i)=image;
    %     afterProcessingImages(:,:,i)=bwmorph(afterProcessingImages(:,:,i),'remove'); % fill holes process
end
end
% %% new code
% 新的程序虽然保留了很多不能分开的地方，但也带来了很多不必要的错误，权衡后还是用老的code
% function [afterProcessingImages,imageProcessingInfo]=myImageProcessing(beforeProcessingImages,picType) %this function can segreate orignal images as you want and returen the mask images
% imageType='uint8'; %here you can chenge your image type
% % here are three types of crop method
%
% % % 1. x or y crop
% % xCrop=270:2160;
% % cropInfo=false(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
% % cropInfo(xCrop,:)=1;
% % beforeProcessingImages=beforeProcessingImages(xCrop,:,:);
% %
% % % 2.x and y crop
% % xCrop=1000:1200;
% % yCrop=700:900;
% % cropInfo=false(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
% % cropInfo(xCrop,yCrop)=1;
% % beforeProcessingImages=beforeProcessingImages(xCrop,yCrop,:);
% %
% % 3. x and y remove
% xCrop=1900:2160;
% yCrop=2000:2560;
% cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
% cropInfo(xCrop,yCrop)=0;
% beforeProcessingImages(xCrop,yCrop,:)=90;
%
% % 4. no need to crop
% % cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
%
% gaussianFilter=fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
% % edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;  %here create edageFilter
% edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
% afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
% brightImages=false(size(beforeProcessingImages));
% cropArea=~cropInfo;
% if strcmp(picType,'w')
%     grayThresh=20;
%     areaThreshold=60;  % proper 60, overlap 60+40
%     maxIntensity=220;      %%%%%%% long time 100X neo
%     bacteriaLight=150;
%     spaceDist=6;
% %     grayThresh=65;
% %     areaThreshold=100;  % proper 60, overlap 60+40
% %   maxIntensity=180;      %%%%%%% agar 1
% %     grayThresh=90;
% %     areaThreshold=120;  % proper 60, overlap 60+40
% %     maxIntensity=240;      %%%%%%% agar 2
% %     grayThresh=80;
% %     areaThreshold=100;
% %     maxIntensity=220;     % general
% end
% if strcmp(picType,'b')
% %     grayThresh=90;
% %     areaThreshold=10;
% %     maxIntensity=170;     %%% 100X black
%     grayThresh=90;
%     areaThreshold=10;
%     maxIntensity=150;
% end
%
% imageProcessingInfo.cropInfo=cropInfo;
% imageProcessingInfo.grayThresh=grayThresh;
% imageProcessingInfo.areaThreshold=areaThreshold;
% imageProcessingInfo.maxIntensity=maxIntensity;
% imageProcessingInfo.bacteriaLight=bacteriaLight;
% imageProcessingInfo.spaceDist=spaceDist;
%
% parfor iframe=1:size(beforeProcessingImages,3)
%     afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
%     afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
%     afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
%     afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/255);%1.2 for the image with few particles
%     afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
%     %     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
%     brightImages(:,:,iframe)=afterProcessingImages(:,:,iframe)==1 & beforeProcessingImages(:,:,iframe)>=bacteriaLight;
%     brightImages(:,:,iframe)=bwareaopen(brightImages(:,:,iframe),10);
%     bwDistance=bwdist(brightImages(:,:,iframe));
%     afterProcessingImages(:,:,iframe)=logical(afterProcessingImages(:,:,iframe)) & bwDistance<=spaceDist;
%     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'open')
%     afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold);
%     image=afterProcessingImages(:,:,iframe);
%     cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList','MeanIntensity');
%     for iCC=1:size(cc,1)
%         if cc(iCC).MaxIntensity<=maxIntensity
%             image(cc(iCC).PixelIdxList)=0;
%         end
%     end
%     cc=regionprops(logical(1-afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'PixelIdxList','MeanIntensity','FilledArea');
%     for iCC=1:size(cc,1)
%         if cc(iCC).MeanIntensity==255 && cc(iCC).FilledArea<=areaThreshold
%             image(cc(iCC).PixelIdxList)=1;
%         end
%     end
%     afterProcessingImages(:,:,iframe)=image;
% %     afterProcessingImages(:,:,iframe)=bwmorph(afterProcessingImages(:,:,iframe),'remove'); % fill holes process
%     afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe) | cropArea);
% end
% afterProcessingImages=logical(afterProcessingImages);
% end
%
