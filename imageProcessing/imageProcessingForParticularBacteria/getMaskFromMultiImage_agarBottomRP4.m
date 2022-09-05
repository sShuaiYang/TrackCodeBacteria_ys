function [maskInfo,maskImage]=getMaskFromMultiImage_agarBottomRP4(image)
% ���EMCCD������ĺ���channelӫ��ͼ�����ڹ�ǿ�������������µ�ʶ������
% ���ԴӶ��channelͬʱ�õ�mask��Ϣ
% ���Դ�ĸ�˹�������һ��Ȧ����̬���Ե�ϸ��
% ������С�ĸ�˹������ڶ���Ȧ����ǿ������ϸ��
oriImage=image-2000;
maskInfo.pixelIdxInfo=[];
maskInfo.centroidInfo=[];
maskImage=false(size(image));
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
for threShold=[80]
    [maskImage1,oriImage]=imageProcessingForChannel1LowSigma(oriImage,threShold);
    cc=regionprops(maskImage1,'PixelIdxList','Centroid');
    for iPiece=1:numel(cc)
        n=n+1;
        maskInfo.pixelIdxInfo{n,1}=cc(iPiece).PixelIdxList;
        maskInfo.centroidInfo(n,:)=cc(iPiece).Centroid;
    end
    maskImage1=bwmorph(maskImage1,'thin','inf');
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
minArea=50;
image=oriImage;
gaussianFilter=fspecial('gaussian',[5,5],1); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=imfilter(image,edgeFilter); %use edgeFilter process
% image=imfilter(image,gaussianFilter); % use Guassian blur filter process
image=im2bw(image,threShold/65535);
image=imclearborder(image);
image=imfill(image,'holes');
% image=bwmorph(image,'close');
image=bwmorph(image,'thin','inf');
image=imdilate(image,ones(3));
cc=regionprops(image,oriImage,'MajorAxisLength','MinorAxisLength','FilledArea','PixelIdxList','MeanIntensity');
for i=1:numel(cc)
    if cc(i).FilledArea<=minArea || cc(i).MeanIntensity<=minIntensity
        image(cc(i).PixelIdxList)=0;
    end
end
image1=image;
image1=imdilate(image1,ones(4));
oriImage(image1)=0;
end