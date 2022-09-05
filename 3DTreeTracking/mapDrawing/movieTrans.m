function movieTrans(fileName,scale)
% 将需要压缩的tiff文件拖到新建的文件夹中，改名为ori.tiff
image=import_tiff_stack([fileName,'\ori.tif']);
for i=1:size(image,4)
    image1(:,:,:,i)=imresize(image(:,:,:,i),scale);  
end
mkdir([fileName,'\result'])
for i=1:size(image,4)
    imwrite(image1(:,:,:,i),[fileName,'\result\',num2str(i),'.tif'],'WriteMode','append');
end
end
function imageStack= import_tiff_stack( fname )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
warning off all
infoImage=imfinfo(fname);
frameNum=size(infoImage,1);
imageWith=infoImage(1).Width;
imageHeight=infoImage(1).Height;
imageBit=infoImage(1).BitsPerSample;
imageBit=strcat('uint',num2str(imageBit));
imageStack=zeros(imageHeight,imageWith,3,'uint8');
imageCurrent=Tiff(fname,'r');
for iframe=1:frameNum
    imageCurrent.setDirectory(iframe);
    imageStack(:,:,:,iframe)= imageCurrent.read();
end
imageCurrent.close();
end