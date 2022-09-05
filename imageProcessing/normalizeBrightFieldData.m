function [imageStack1] = normalizeBrightFieldData(imageStack)
% ���ڰѲ�ͬ��imageStack��ǿ��һ������ͼ����ĳ���
for i=1:size(imageStack,3)
    image=imageStack(:,:,i);
    imageStack1(:,:,i)=imadjust(image,[0,max(max(double(image)))/65535],[0,1000/65535]);
end
end

