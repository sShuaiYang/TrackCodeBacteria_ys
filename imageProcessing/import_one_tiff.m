function imageStack= import_one_tiff( fname , frameNum)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
warning off all
if nargin==1
    frameNum=1;
end
infoImage=imfinfo(fname);
imageWith=infoImage(1).Width;
imageHeight=infoImage(1).Height;
imageBit=infoImage(1).BitsPerSample;
imageBit=strcat('uint',num2str(imageBit));
imageStack=zeros(imageHeight,imageWith,imageBit);
imageCurrent=Tiff(fname,'r');
imageCurrent.setDirectory(frameNum);
imageStack= imageCurrent.read();
imageCurrent.close();
end