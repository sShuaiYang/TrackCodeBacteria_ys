function rudeMaskForMultiChannelImage()
dirFile=uigetdir();
dirImage1=[dirFile,'\image1.tif'];
dirImage2=[dirFile,'\image2.tif'];
image1=import_tiff_stack(dirImage1)-110;
image2=import_tiff_stack(dirImage2)-110;
for i=1:size(image1,3)
    image1(:,:,i)=imadjust(image1(:,:,i),[0,max(max(double(image1(:,:,i))))/65535],[0,0.6]);
    image2(:,:,i)=imadjust(image2(:,:,i),[0,max(max(double(image2(:,:,i))))/65535],[0,0.6]);
end
imageAll=image1+image2;
imageAll=im2uint8(imageAll);
% imageAll=255-imageAll;
testFile=[dirFile,'\properResult'];
mkdir(testFile)
cd(testFile)
fluoImageProcessing(imageAll,testFile);
% for i=1:size(image,3)
%     imageAll(:,:,i)=imadjust(imageAll(:,:,i));
% end
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
imageStack=zeros(imageHeight,imageWith,imageBit);
imageCurrent=Tiff(fname,'r');
for iframe=1:frameNum
    imageCurrent.setDirectory(iframe);
    imageStack(:,:,iframe)= imageCurrent.read();
end
imageCurrent.close();
end
function fluoImageProcessing(imageStack,testFile)
for iStack=22:size(imageStack,3)
    image=imageStack(:,:,iStack);
    edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;
    gaussianFilter1=fspecial('gaussian',[5, 5],5);
    gaussianFilter2=fspecial('gaussian',[3, 3],2);
    %     image=uint8(double(image)/double(max(max(image)))*65535);
    image=imfilter(image,gaussianFilter1); % use Guassian blur filter process
    image=imfilter(image,edgeFilter); %use edgeFilter process
    image=imfilter(image,gaussianFilter1);
    profileImage=im2bw(image,12/255);
    profileImage=imclearborder(profileImage);
    profileImage=imfill(profileImage,'holes');
    profileImage=bwareaopen(profileImage,40);
    image(~profileImage)=255;
    finalImage=bwmorph(profileImage,'remove');
    
    n=0;
    for thre=[240,230,220]
        n=n+1;
        image1=255-image;
        image1=im2bw(image1,thre/255);
        finalImage(image1)=1;
        imageA=imadjust(imageStack(:,:,iStack),[0,double(max(max(imageStack(:,:,iStack))))/255/2],[0,1]);
        imageB=imageA;
        imageC=imageA;
        imageA(finalImage)=255;
        imageB(finalImage)=0;
        imageC(finalImage)=0;
%         imageC(finalImage)=255;
        imageF=cat(3,imageA,imageB,imageC);
        imageF=uint8(imageF);
        imwrite(imageF,[num2str(iStack),'-',num2str(n),'.tif']);
    end
end
end