function batchImageProcessJob()
dirImage='K:\FanJIN\raw data\1028202010_jds20_filM\tiff2matlab';
addpath(dirImage);
clc;
nameList=dir(dirImage);
for i=1:(length(nameList)-2)    
    saveName='maskImagesTiffStack.tif';
    dirResultSave='K:\FanJIN\raw data\1028202010_jds20_filM\maskImage'; %put the folder of the save results
    dirImage='K:\FanJIN\raw data\1028202010_jds20_filM\tiff2matlab';%put the folder of the matlab image
    addpath(dirImage);
    nameList=dir(dirImage);
    jobNumber=strcat('start Job', num2str(i));
    fileName=nameList(i+2).name;
    disp(jobNumber);
    disp(strcat('load file',fileName));
    tic;disp(strcat('1.load:',fileName,'...'));tempImage=load(fileName);toc;
    imageStack=tempImage.imageStack;
    clear tempImage;
    disp('image Processing')
    maskImages=myImageProcessing(imageStack);
    disp('save MaskImage');
    saveFile1=strcat(dirResultSave,'\',saveName);
    tic;
    for iFrame=1:size(maskImages,3)              
       imwrite(maskImages(:,:,iFrame),saveFile1,'tif','Compression','none','WriteMode','append');
    end
    toc;
    clear all;
end
end