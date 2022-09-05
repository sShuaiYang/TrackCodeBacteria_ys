function [maskImage,image1]=getMaskFromMultiImage_Ppsl_nopslA(image)
cutNum=32;
for i=1:cutNum
    image1(:,:,i)=image(:,:,i)-120;
    image1(:,:,i)=imadjust(image1(:,:,i),[double(min(min(image1(:,:,i))))/65535 double(max(max(image1(:,:,i))))/65535],[0 1]);
    [~,maskImage(:,:,i)]=getMaskFromMultiImage_HJDNewData(image1(:,:,i));
end
end
function maskImage=easyProcessing(oriImage)
backGround=150;
minArea=20;
threShold=80;
minIntensity=200;
image=oriImage;
gaussianFilter=fspecial('gaussian',[5,5],3); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=imfilter(image,edgeFilter); %use edgeFilter process
image=im2bw(image,threShold/65535);
image=bwmorph(image,'thin','inf');
image=imdilate(image,ones(3));
cc=regionprops(image,oriImage,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','MeanIntensity');
for i=1:numel(cc)
    if cc(i).FilledArea<=minArea || cc(i).MeanIntensity<=minIntensity
        image(cc(i).PixelIdxList)=0;
    end
end
maskImage=imclearborder(image);
end
function [maskInfo,maskImage]=getMaskFromMultiImage_HJDNewData(image)
% 针对EMCCD拍摄出的红绿channel荧光图像，由于光强明暗不均而导致的识别困难
% 尝试从多个channel同时得到mask信息
% 用稍大的高斯卷积，第一步圈出形态明显的细菌
% 再用稍小的高斯卷积，第二部圈出光强较弱的细菌
oriImage=image;
maskInfo.pixelIdxInfo=[];
maskInfo.centroidInfo=[];
maskImage=false(512,512);
n=0;
% for threShold=[2000,1500,1000,700,500,400,300,200,2000,1500,1000,700,500,400,300,200]
%     [maskImage1,oriImage]=imageProcessingForChannel1HighSigma(oriImage,threShold);
%     cc=regionprops(maskImage1,'PixelIdxList','Centroid');
%     for iPiece=1:numel(cc)
%         n=n+1;
%         maskInfo.pixelIdxInfo{n,1}=cc(iPiece).PixelIdxList;
%         maskInfo.centroidInfo(n,:)=cc(iPiece).Centroid;
%     end
%     maskImage1=bwmorph(maskImage1,'thin','inf');
%     maskImage=maskImage | maskImage1;
% end
% for threShold=[200]
%     [maskImage1,oriImage]=imageProcessingForChannel1LowSigma(oriImage,threShold);
%     cc=regionprops(maskImage1,'PixelIdxList','Centroid');
%     for iPiece=1:numel(cc)
%         n=n+1;
%         maskInfo.pixelIdxInfo{n,1}=cc(iPiece).PixelIdxList;
%         maskInfo.centroidInfo(n,:)=cc(iPiece).Centroid;
%     end
%     maskImage1=bwmorph(maskImage1,'thin','inf');
%     maskImage=maskImage | maskImage1;
% end
for threShold=[10000]
    [maskImage1,oriImage]=imageProcessingForChannel1LowSigma1(oriImage,threShold);
    cc=regionprops(maskImage1,'PixelIdxList','Centroid');
    for iPiece=1:numel(cc)
        n=n+1;
        maskInfo.pixelIdxInfo{n,1}=cc(iPiece).PixelIdxList;
        maskInfo.centroidInfo(n,:)=cc(iPiece).Centroid;
    end
%     maskImage1=bwmorph(maskImage1,'thin','inf');
    maskImage=maskImage | maskImage1;
end
for threShold=[10000]
    [maskImage1,oriImage]=imageProcessingForChannel1LowSigma2(oriImage,threShold);
    cc=regionprops(maskImage1,'PixelIdxList','Centroid');
    for iPiece=1:numel(cc)
        n=n+1;
        maskInfo.pixelIdxInfo{n,1}=cc(iPiece).PixelIdxList;
        maskInfo.centroidInfo(n,:)=cc(iPiece).Centroid;
    end
%     maskImage1=bwmorph(maskImage1,'thin','inf');
    maskImage=maskImage | maskImage1;
end
end
function [image,oriImage]=imageProcessingForChannel1HighSigma(oriImage,threShold)
backGround=150;
minArea=30;
image=oriImage;
gaussianFilter=fspecial('gaussian',[5,5],20); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=imfilter(image,edgeFilter); %use edgeFilter process
% image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=im2bw(image,threShold/65535);
image=bwmorph(image,'open');
cc=regionprops(image,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList');
for i=1:numel(cc)
    if cc(i).FilledArea<=minArea || abs(cc(i).MajorAxisLength*cc(i).MinorAxisLength/4*pi-cc(i).FilledArea)/cc(i).FilledArea>=0.3;
        image(cc(i).PixelIdxList)=0;
    end
end
image1=image;
image1=bwmorph(image1,'dilate',1);
oriImage(image1)=0;
end
function [image,oriImage]=imageProcessingForChannel1LowSigma(oriImage,threShold)
backGround=150;
minIntensity=300;
minArea=30;
image=oriImage;
gaussianFilter=fspecial('gaussian',[3,3],0.5); %here create Gaussian Blur Filter
edgeFilter=(ones(3,3)).*-1;edgeFilter(2,2)=8;
image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=imfilter(image,gaussianFilter);
image=imfilter(image,gaussianFilter);
image=imfilter(image,gaussianFilter);
image=imfilter(image,gaussianFilter);
image=imfilter(image,edgeFilter); %use edgeFilter process
% image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=im2bw(image,threShold/65535);
image=bwareaopen(image,minArea);
image=bwmorph(image,'thin','inf');
image=imdilate(image,ones(3));
cc=regionprops(image,oriImage,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','MeanIntensity');
for i=1:numel(cc)
    if cc(i).MeanIntensity<=minIntensity
        image(cc(i).PixelIdxList)=0;
    end
end
image1=image;
image1=imdilate(image1,ones(5));
oriImage(image1)=0;
image=imclearborder(image);
% oriImage=imadjust(oriImage,[double(min(min(oriImage)))/65535 double(max(max(oriImage)))/65535],[0,1]);
end
function [image,oriImage]=imageProcessingForChannel1LowSigma1(oriImage,threShold)
backGround=150;
minIntensity=6000;
minArea=20;
image=oriImage;
gaussianFilter=fspecial('gaussian',[5,5],1); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=imfilter(image,edgeFilter); %use edgeFilter process
% image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=im2bw(image,threShold/65535);
image=bwmorph(image,'thin','inf');
image=imdilate(image,ones(3));
cc=regionprops(image,oriImage,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','MeanIntensity');
for i=1:numel(cc)
    if cc(i).FilledArea<=minArea || cc(i).MeanIntensity<=minIntensity
        image(cc(i).PixelIdxList)=0;
    end
end
image1=image;
image1=bwmorph(image1,'dilate',2);
oriImage(image1)=0;
image=imclearborder(image);
end
function [image,oriImage]=imageProcessingForChannel1LowSigma2(oriImage,threShold)
backGround=150;
minIntensity=3000;
minArea=20;
image=oriImage;
gaussianFilter=fspecial('gaussian',[5,5],5); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=imfilter(image,edgeFilter); %use edgeFilter process
% image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=im2bw(image,threShold/65535);
cc=regionprops(image,oriImage,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','MeanIntensity');
for i=1:numel(cc)
    if cc(i).FilledArea<=minArea || cc(i).MeanIntensity<=minIntensity
        image(cc(i).PixelIdxList)=0;
    end
end
image=bwmorph(image,'thin','inf');
image=imdilate(image,ones(3));
image1=image;
image1=bwmorph(image1,'dilate',2);
oriImage(image1)=0;
image=imclearborder(image);
end