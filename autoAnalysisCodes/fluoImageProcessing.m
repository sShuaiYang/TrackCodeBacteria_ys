function [maskImage,imageStack]=fluoImageProcessing(imageStack)
% minIntensity=300;
% imageStack=imageStack-200;
for iStack=1:size(imageStack,3)
    image=imageStack(:,:,iStack);
    pixelInfo=image(:);
    pixelInfo=sort(pixelInfo);
%     backGround=mean(pixelInfo(1:end/3));
    maxPixel=max(pixelInfo);
    image=image-100;
    edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
    gaussianFilter1=fspecial('gaussian',[10, 10],3);
    gaussianFilter2=fspecial('gaussian',[10, 10],2);
    image=uint8(double(image)/double(maxPixel)*255);
    image=imfilter(image,gaussianFilter1); % use Guassian blur filter process
    image=imfilter(image,edgeFilter); %use edgeFilter process
%     image=imfilter(image,gaussianFilter2);
    thre=8;
    image=im2bw(image,thre/255);
    image=imclearborder(image);
    image=bwareaopen(image,100,4);
    image=logical(image);
    image=imfill(image,'holes');
    maskImage(:,:,iStack)=image;
    %     maskImage(:,:,iStack)=bwmorph(image,'close');
end
end