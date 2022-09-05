function [image1,image2]=imCameraTransform(image1,image2,transCamera,bestPosition)
imageTran=imtransform(image2,transCamera,'XData',[1 size(image2,1)],'YData',[1 size(image2,2)]);
imageStack=cat(3,image1,imageTran);
imageStack=imageCorrectionWithBestPosition(imageStack,bestPosition);
image2=imageStack(:,:,2);
end