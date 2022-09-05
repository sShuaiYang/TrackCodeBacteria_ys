function montageImageGet(dirFile)
dirOpen=strcat(dirFile,'\','original data');
dirSave=strcat(dirFile,'\','allResult');
mkdir(dirSave);
filein=strcat(dirFile,'\','information\montageInfo.txt');
nameList=dir(dirOpen);

stackNum=700;
backGroundGFP=import_tiff_stack(strcat(dirFile,'\backGroundGFP.tif'));
backGroundRFP=import_tiff_stack(strcat(dirFile,'\backGroundRFP.tif'));

montageInfoRead=scaleInforeadMontage(filein);
wavelength=montageInfoRead.WavelengthNum;
montageNum=montageInfoRead.MontageNum;
timeNum=montageInfoRead.TimeNum;

if (length(nameList)-2)/(wavelength*montageNum)~=timeNum
    timeNum=(length(nameList)-2)/(wavelength*montageNum);
end

for iMontage=1:montageNum
    disp(iMontage);
    tic;
    motangeSave=strcat(dirSave,'\',nameList((iMontage-1)*wavelength+3).name(end-14:end-10));
    mkdir(motangeSave);
    for iStack=1:ceil(timeNum/stackNum)
        
        dirSaveGFP=strcat(motangeSave,'\imageGFP');
        dirSaveRFP=strcat(motangeSave,'\imageRFP');
        mkdir(dirSaveGFP);
        mkdir(dirSaveRFP);
        if iStack==ceil((timeNum/stackNum))
            for iTime=1:timeNum-(iStack-1)*stackNum
                timeStackNum=(iStack-1)*stackNum+iTime;
                imageGFP(:,:,timeStackNum)=import_tiff_stack(strcat(dirOpen,'\',nameList((timeStackNum-1)*montageNum*wavelength+(iMontage-1)*wavelength+3).name));
                imageRFP(:,:,timeStackNum)=import_tiff_stack(strcat(dirOpen,'\',nameList((timeStackNum-1)*montageNum*wavelength+(iMontage-1)*wavelength+4).name));
            end
            imageGFP=FloBackGroundCorretion(imageGFP(1:490,:,:),backGroundGFP(1:490,:,:));
            imageRFP=FloBackGroundCorretion(imageRFP(1:490,:,:),backGroundRFP(1:490,:,:));
            save(strcat(dirSaveGFP,'\imageGFP',num2str(iStack,'%02d'),'.mat'),'imageGFP','-v7.3');
            save(strcat(dirSaveRFP,'\imageRFP',num2str(iStack,'%02d'),'.mat'),'imageRFP','-v7.3');
            clear imageGFP;
            clear imageRFP;
        else
            for iTime=1:stackNum
                timeStackNum=(iStack-1)*stackNum+iTime;
                imageGFP(:,:,timeStackNum)=import_tiff_stack(strcat(dirOpen,'\',nameList((timeStackNum-1)*montageNum*wavelength+(iMontage-1)*wavelength+3).name));
                imageRFP(:,:,timeStackNum)=import_tiff_stack(strcat(dirOpen,'\',nameList((timeStackNum-1)*montageNum*wavelength+(iMontage-1)*wavelength+4).name));
            end
            imageGFP=FloBackGroundCorretion(imageGFP(1:490,:,:),backGroundGFP(1:490,:,:));
            imageRFP=FloBackGroundCorretion(imageRFP(1:490,:,:),backGroundRFP(1:490,:,:));
            save(strcat(dirSaveGFP,'\imageGFP',num2str(iStack,'%02d'),'.mat'),'imageGFP','-v7.3');
            save(strcat(dirSaveRFP,'\imageRFP',num2str(iStack,'%02d'),'.mat'),'imageRFP','-v7.3');
            clear imageGFP;
            clear imageRFP;
        end
    end
    toc;
end

end
function dataOut=scaleInforeadMontage(filein)
% input filein is the URL of the information.txt
fidin=fopen(filein,'r');
nline=0;
while ~feof(fidin) % 判断是否为文件末尾
    tline=fgetl(fidin); % 从文件读行
    nline=nline+1;
    switch nline
        case 5
            % find the line *** x : 1721 * 0.061728 : um ***
            info=textscan(tline,'%s%s%f%s%f');
            scaleInfo=info{5};
            xSize=info{3};
        case 6
            info=textscan(tline,'%s%s%f%s%f');
            ySize=info{3};
        case 7
            info=textscan(tline,'%s%s%f%s%f');
            WavelengthNum = info{3};
        case 8
            info=textscan(tline,'%s%s%f%s%f');
            MontageNum = info{3};
        case 9
            info=textscan(tline,'%s%s%f%s%f');
            TimeNum = info{3};
    end
    %find the line *** Repeat T - 20000 times (5 sec) ***
    isRepeatLine=strfind(tline,'Repeat');
    if ~isempty(isRepeatLine)
        info=textscan(tline,'%s %s %s %f %s %s %s');
        if ~isempty(info{7})
        switch info{7}{1}(1:end-1)
            case 'sec'
                timeInterval=str2num(info{6}{1}(2:end));                
            case 'ms'
                timeInterval=str2num(info{6}{1}(2:end))/1000;
            case 'min'
                timeInterval=str2num(info{6}{1}(2:end))*60;
            case 'h'
                timeInterval=str2num(info{6}{1}(2:end))*3600;
        end
        end
    end
end
timeInterval=1;
fclose(fidin);
dataOut.rcSize=[ySize,xSize];
dataOut.scaleInfo=scaleInfo;
dataOut.timeInterval=timeInterval;
dataOut.WavelengthNum=WavelengthNum;
dataOut.MontageNum=MontageNum;
dataOut.TimeNum=TimeNum;
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
function imageGFP=FloBackGroundCorretion(imageGFP,backGround)
backGround=double(backGround);
fluoBack=backGroundGet(imageGFP);
averageBack=mean(mean(backGround));
for iframe=1:size(imageGFP,3)
%     dispFrame(iframe);
    imageGFP=imageGFP(1:size(backGround,1),1:size(backGround,2),:);
    imageTemp=double(imageGFP(:,:,iframe));
%     averageGFP=mean(mean(imageTemp));
    imageGFP(:,:,iframe)=(imageTemp-fluoBack)./backGround*averageBack;    % change by jzy 10.14
end
% fprintf('\n');
imageGFP=uint16(imageGFP);
end
function dispFrame(iConnect)
lengthB=length(num2str(iConnect-1));
back='';
for i=1:lengthB
    back=strcat(back,'\b');
end
fprintf(strcat(back,num2str(iConnect)));
end
function  [gBack]=backGroundGet(gfpImage)
%首先要对红色荧光和绿色荧光图像去除背景
    gImage=gfpImage;
    gIndex=gImage(:);
    gIndex=sort(gIndex);
    gBack=mean(gIndex(1:round(end/6)));
end