function overlapBFandRGB(originalPicURL,trajectoryImage,startFrame,stepFrame,endFrame)
dirImages=originalPicURL;
savePath=uigetdir();
cd(dirImages);
nameList=dir(dirImages);
vars = whos('-file', nameList(3).name);
stackSize=vars.size(3);
stackNum=fix(startFrame/(stackSize));
iImage=1;
if stackNum~=startFrame/(stackSize)
    fileName=nameList(stackNum+3).name;
else
    fileName=nameList(stackNum+2).name;
end
imageStack=loadImageStack(fileName);
for iframe=startFrame:stepFrame:endFrame;
    disp(iframe);
    stackNumNext=fix(iframe/(stackSize));
    if stackNumNext==iframe/(stackSize)
        stackNumNext=stackNumNext-1;
    end
    if stackNumNext~=stackNum
        clear imageStack;
        fileName=nameList(stackNumNext+3).name;
        imageStack=loadImageStack(fileName);
        stackNum=stackNumNext;
    end
    overlapAndSave(imageStack(:,:,iframe-stackNum*stackSize),trajectoryImage(:,:,:,iImage),iframe,savePath);
    iImage=iImage+1;
end
rmpath(dirImages);
end
function overlapAndSave(originPic,trajectoryPic,i,savePath)
finalImage1=originPic;
finalImage2=originPic;
finalImage3=originPic;
red=trajectoryPic(:,:,1);
green=trajectoryPic(:,:,2);
blue=trajectoryPic(:,:,3);
finalImage1(blue~=0)=0;
finalImage2(blue~=0)=0;
finalImage3(blue~=0)=255;
finalImage1(red~=0)=255;
finalImage2(red~=0)=0;
finalImage3(red~=0)=0;
finalImage1(green~=0&red==0)=0;
finalImage2(green~=0)=255;
finalImage3(green~=0)=0;
finalImage=cat(3,finalImage1,finalImage2,finalImage3);
imwrite(finalImage,strcat(savePath,'\demo',num2str(i)),'tif');
end
function imageStack=loadImageStack(fileName)
tempImage=load(fileName);
imageStack=tempImage.imageStack;
end

