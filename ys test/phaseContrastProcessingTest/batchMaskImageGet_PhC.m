function batchMaskImageGet_PhC()
%此函数暂时适用于micromanger 拍摄的相差图像stack处理
dirFile='F:\2020-04-02 pvdA_IP31_100x Temp 30 fastshot PHC ys';
fieldList=dir(dirFile);
for iField=1:length(fieldList)-2
    if ~strcmp(fieldList(iField+2).name(1:5),'field')
        continue
    end
    t0=clock;
    disp('start PhC image processing ')
    dirField=[dirFile,'\',fieldList(iField+2).name];
%     dirTrack=[dirField,'\Tracking'];
%     mkdir(dirTrack);
%     dirMaskCheck= [dirField,'\maskCheck'];
%     mkdir(dirMaskCheck);
    dirMaskSave=strcat(dirField,'\tiff2matlab');
    mkdir(dirMaskSave);
    dirImage=[dirField,'\','BF'];
    imageList=dir(dirImage);
    %去除Thumbs.db文件影响文件读取
    if isfile(strcat(dirImage,'\Thumbs.db'))
        imageList = thumbsFileRemove(imageList);
    end
    
    disp(fieldList(iField+2).name);
    
    j = 1;
    for iImage = 1:length(imageList)-2
        if ~strcmp(imageList(iImage+2).name(end-2:end),'tif')
            continue
        end
        imageStack = import_tiff_stack( strcat(dirImage,'\',imageList(iImage+2).name) );
        [maskImages] = phaseContrastImageProcessing_ys(imageStack);
        save([dirMaskSave,'\',num2str(j,'%05.f')],'maskImages','-v7.3')
        j = j+1;
        clear imageStack maskImages
    end
    disp(['耗时:',num2str(etime(clock,t0)),'秒']);
end
end
%%
function  fileList = thumbsFileRemove(fileList)
templogic = true (1, numel (fileList));
for iFile = 1: numel (fileList)
    
    if strcmp (fileList(iFile).name, 'Thumbs.db')
        templogic(iFile) = false;
        fileList = fileList(templogic);
        return
    end
    
end
end