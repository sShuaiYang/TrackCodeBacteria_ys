function  [gfpImage]=backGroundSubstract(gfpImage)
%����Ҫ�Ժ�ɫӫ�����ɫӫ��ͼ��ȥ������,���ܶ�ʱ���ܻ����Щ������
for i=1:size(gfpImage,3)
    gImage=gfpImage(:,:,i);
    gIndex=gImage(:);
    gIndex=sort(gIndex);
    gBack=mean(gIndex(1:end*1/3));
    gfpImage(:,:,i)=gfpImage(:,:,i)-gBack;
end
end

% function  [gfpImage,maskImage]=backGroundSubstract(gfpImage)
% %��������������ٳ�ȥ����
% [maskImage,gfpImage]=fluoImageProcessing(gfpImage);
% end
% function [maskImage,imageStack]=fluoImageProcessing(imageStack)
% for iStack=1:size(imageStack,3)
%     image1=imageStack(:,:,iStack);
%     pixelInfo=image1(:);
%     pixelInfo=sort(pixelInfo);
%     backGround=mean(pixelInfo(1:end/3));
%     image=image1-backGround;
%     edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
%     gaussianFilter1=fspecial('gaussian',[5, 5],20);
%     gaussianFilter2=fspecial('gaussian',[3, 3],2);
%     image=uint8(double(image)/double(max(max(image)))*255);
%     image=imfilter(image,gaussianFilter1); % use Guassian blur filter process
%     image=imfilter(image,edgeFilter); %use edgeFilter process
%     image=imfilter(image,gaussianFilter2);
%     image=im2bw(image,10/255);
%     image=bwareaopen(image,15);
%     image=logical(image);
%     imageBack=logical(true-image);
%     backGround=mean(image1(imageBack));
%     imageStack(:,:,iStack)=imageStack(:,:,iStack)-backGround;
%     image=imclearborder(image);
%     image=imfill(image,'holes');
%     maskImage(:,:,iStack)=image;
% end
% end