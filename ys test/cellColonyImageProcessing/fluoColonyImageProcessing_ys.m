function [processedImages] = fluoColonyImageProcessing_ys(fluoImages,thresh)
% Shuai Yang 根据[processedImages] = fluoImageProcessing_ystest(fluoImages,thresh)
%修改来
% convert to 8bit images
convertedImages = uint8 (zeros (size(fluoImages)));
parfor iframe = 1:size(fluoImages,3)
    I0 = fluoImages(:,:, iframe);
    I0 = substractBackGround(I0);%扣除背景
    I0 = imgaussfilt(I0,2);
    I0 = rescale(double(I0))*255;
    I0 = uint8(I0);
    convertedImages (:,:,iframe) = I0;
end

fluoImages = convertedImages;
% imageType = 'uint8';

% thresh=20;% intensity threshold
if nargin<2 || isempty(thresh)
    I1 = fluoImages(:,:,1);
    thresh = imageThresholdGet(I1,'Otsu');
    thresh = thresh/1.3;% default equals to 2
else
    thresh = thresh/255;
end
se = strel('disk',20);
areaThreshold = 100;  % colony 面积pixel

processedImages = false(size(fluoImages));% here intiatlize the maskImages stacks
parfor iframe = 1:size(fluoImages,3)
    processedImages(:,:,iframe) = imbinarize(fluoImages(:,:,iframe),thresh);
    processedImages(:,:,iframe) = bwareaopen(processedImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    processedImages(:,:,iframe) = imfill(processedImages(:,:,iframe),'holes');
    processedImages(:,:,iframe) = bwmorph(processedImages(:,:,iframe),'open');
    processedImages(:,:,iframe) = bwareaopen(processedImages(:,:,iframe),areaThreshold);
    processedImages(:,:,iframe) = imerode(processedImages(:,:,iframe),se);
end

% result check
processedImages = logical(processedImages);
figure,imshow(imoverlay(uint8(rescale(double(fluoImages(:,:,1)))*255),processedImages(:,:,1),[1,0,0]));
bw_perm = bwperim(processedImages(:,:,1));
figure,imshowpair(bw_perm,uint8(rescale(double(fluoImages(:,:,1)))*255));

end
%%
function  T = imageThresholdGet(I,Method)
%I: 8bit image
%method 'Iterative','Otsu'
if nargin<2 || isempty(Method)
    Method = 'Otsu';
end

if strcmp(Method,'Otsu')
    T = graythresh(im2double(I));
end

if strcmp(Method,'Iterative')
    f = im2double(I);
    T = 0.5*(min(f(:))+max(f(:)));
    done=false;
    while ~done
        g = f>=T;
        Tn = 0.5*(mean(f(g))+mean(f(~g)));
        done = abs(T-Tn)<0.1;
        T = Tn;
    end
    T = T/255;
end
end
%%
function I = substractBackGround(I)
%用于减除图像的背景 Shuai Yang 2020.09.09
sz = size(I);
I0 = I(round(sz(1)/2-sz(1)/4):round(sz(1)/2+sz(1)/4),round(sz(2)/2-sz(2)/4):round(sz(2)/2+sz(2)/4));
I0 = double(I0);
I0 = sort(I0(:));
pixelSpace = 1:100:numel(I0);%每100个pixel计算平均值
M = zeros(numel(pixelSpace)-1,1);
for i = 1:numel(pixelSpace)-1
    M(i) = mean(I0(pixelSpace(i):pixelSpace(i+1)));
end

n = 5;
CV = std(M(1:n))/mean(M(1:n));
while CV > 0.03 % equal CV > 0.03 && n >= 1
    n = n-1;
    CV = std(M(1:n))/mean(M(1:n));
end
BG = mean(M(1:n));
I = I - BG;
end