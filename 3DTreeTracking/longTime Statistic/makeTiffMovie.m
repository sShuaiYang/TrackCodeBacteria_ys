function makeTiffMovie(dirFile,beginFrame,endFrame,interval)
% 2014-12-1
% 从背景矫正后的图片tiff2matlab中，每隔100张（5分钟保存一张图片）,选择出来保存
if nargin==0
    dirFile=uigetdir();
end
dirTiff=strcat(dirFile,'\original data');
dirSave=strcat(dirFile,'\movie');
dirBackGround=strcat(dirFile,'\backGround.tif');
backGround=import_tiff_stack(dirBackGround);
mkdir(dirSave)
nameList=dir(dirTiff);
cd(dirTiff);
imageStack=import_tiff_stack(nameList(end-1).name);
stackSize=size(imageStack,3);  % 理论上最大视野一个stack是189张，防止有修改视野大小导致图像大小改变的情况
imageStack=import_tiff_stack(nameList(end).name);
allNum=stackSize*(numel(nameList)-3)+size(imageStack,3);
if nargin~=4
    beginFrame=1;
    endFrame=allNum;
    interval=100;   % 若不设置，默认为间隔5分钟
end
timeNum=beginFrame:interval:endFrame;
stackNum=fix((timeNum-1)/stackSize)+1;
smallNum=timeNum-(stackNum-1)*stackSize;
for i=1:max(stackNum)
   if all(stackNum~=i)
       continue
   end
   order=1:size(imageStack,3);
   smallNumChoose=smallNum(stackNum==i);
   for iFrame=1:numel(smallNumChoose)
       image=import_one_tiff(nameList(i+2).name,smallNumChoose(iFrame));
       [image,~]=backGroundCorrection(image,backGround,'16bit');
       clc
       name=strcat(dirSave,'\00000.tif');
       bacNum=(i-1)*stackSize+smallNumChoose(iFrame);
       numDigit=fix(log10(bacNum))+1;
       name(end-3-numDigit:end-4)=num2str(bacNum);
       image=imresize(image,0.25);
       imwrite(image,name);
   end
end
end
