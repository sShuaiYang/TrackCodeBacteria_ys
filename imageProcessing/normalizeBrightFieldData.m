function [imageStack1] = normalizeBrightFieldData(imageStack)
% 用于把不同的imageStack光强归一后用于图像处理的程序
for i=1:size(imageStack,3)
    image=imageStack(:,:,i);
    imageStack1(:,:,i)=imadjust(image,[0,max(max(double(image)))/65535],[0,1000/65535]);
end
end

