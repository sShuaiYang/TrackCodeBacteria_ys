function [processedImages] = fluoImageProcessing_ystest(fluoImages,thresh)
%JZY�ɰ汾������ͼ��������޸ĵõ�ӫ��ͼ����code
% Shuai Yang 
% convert to 8bit images
convertedImages = uint8 (zeros (size(fluoImages)));
parfor iframe = 1:size(fluoImages,3)
    I0 = fluoImages(:,:, iframe);
    I0 = substractBackGround(I0);%�۳�����
    %     k = find(double(I0)>800);
    %     %ֻ�е�ͼ����ӫ��ǿ�ȴ���ĳһ�ض�ֵ ���й�һת��Ϊ8bitͼ��
    %     %  ����uint8ֱ��ת��16bitͼ��>255���ػḳֵ255��С��255�Ĳ���
    %     % Ŀ���Ǳ���û��ϸ������Ұ ӫ��ͼ����������Ա�ͼ ���߹���ϸ���ĳ��� ͼ��Ա�ͼ���ָı�
    %     if numel(k)> 1000 %���ص����Ŀ
    %         %I0 = uint8(double(I0-min(I0(:)))/double(max(I0(:))-min(I0(:)))*255);
    %         I0 = rescale(double(I0))*255;
    %         I0 = uint8(I0);
    %     else
    %         I0 = uint8(I0);
    %     end
    %maximum fliter
    %     I0 = colfilt(I0, [winLength winLength], 'sliding', @max);
    %minimum fliter
    %     I0 = colfilt(I0, [winLength winLength], 'sliding', @min);
    I0 = imgaussfilt(I0,2);
    I0 = rescale(double(I0))*255;
    I0 = uint8(I0);
    convertedImages (:,:,iframe) = I0;
end
fluoImages = convertedImages;

imageType = 'uint8';

% gaussianFilter = fspecial('gaussian',[5, 5],20); %here create Gaussian Blur Filter
gaussianFilter = fspecial('gaussian',[10, 10], 3);% suitable for fluo images
edgeFilter = (ones(5,5)).*-1; edgeFilter(3,3) = 24;  %here create edgeFilter

% thresh=20;% intensity threshold
if nargin < 2 || isempty(thresh)
    I1 = fluoImages(:,:,1);
    I1 = imfilter(I1,gaussianFilter);
    I1 = imfilter(I1,edgeFilter); %use edgeFilter process
    I1 = imfilter(I1,gaussianFilter);
    
    thresh = imageThresholdGet(I1,'Otsu');
    thresh = thresh/1.3;% default equals to 2
else
    thresh = thresh/255;
end
areaThreshold = 80;  % proper 60, overlap 60+40

processedImages = zeros(size(fluoImages),imageType);% here intiatlize the maskImages stacks
parfor iframe = 1:size(fluoImages,3)
    processedImages(:,:,iframe) = imfilter(fluoImages(:,:,iframe),gaussianFilter);
    processedImages(:,:,iframe) = imfilter(processedImages(:,:,iframe),edgeFilter); %use edgeFilter process
    processedImages(:,:,iframe) = imfilter(processedImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    processedImages(:,:,iframe) = imbinarize(processedImages(:,:,iframe),thresh);
    processedImages(:,:,iframe) = bwareaopen(processedImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    processedImages(:,:,iframe) = imclearborder(processedImages(:,:,iframe));
    processedImages(:,:,iframe) = imfill(processedImages(:,:,iframe),'holes');
    processedImages(:,:,iframe) = bwmorph(processedImages(:,:,iframe),'open');
end
processedImages = logical(processedImages);
%���ݳߴ����ɸѡ remove region �Ǳ���
% [afterProcessingImages] = imageMaskFilterBasedOnCellSize(afterProcessingImages);
imshowColorlabel_ys(processedImages(:,:,end));
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
%���ڼ���ͼ��ı��� Shuai Yang 2020.09.09
sz = size(I);
I0 = I(round(sz(1)/2-sz(1)/4):round(sz(1)/2+sz(1)/4),round(sz(2)/2-sz(2)/4):round(sz(2)/2+sz(2)/4));
I0 = double(I0);
I0 = sort(I0(:));
pixelSpace = 1:100:numel(I0);%ÿ100��pixel����ƽ��ֵ
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
%%
function [afterProcessingImages] = imageMaskFilterBasedOnCellSize(afterProcessingImages)
%afterProcessingImages binary image 
parfor i = 1:size(afterProcessingImages,3)
    image = afterProcessingImages(:,:,i);
    cc = regionprops(image,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList');
    for iCC = 1:size(cc,1)
        L = cc(iCC).MajorAxisLength;% cell length
        W = cc(iCC).MinorAxisLength; %cell width
        if cc(iCC).FilledArea <= 100 || L/W <= 1.3
            image(cc(iCC).PixelIdxList) = false;
        end
    end
    afterProcessingImages(:,:,i) = image;
end
end