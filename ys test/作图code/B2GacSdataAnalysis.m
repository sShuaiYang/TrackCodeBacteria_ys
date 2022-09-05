%%B2GacSdataAnalysis_ImTransform
fiexedChannel = 'sfGFP';
dirAll = 'F:\GacS\cxy20181205';


channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道

folderList = dir(dirAll);
for ifolder = 8:numel(folderList)-2
    dirFile =  [dirAll,filesep,folderList(ifolder+2).name];
    disp(folderList(ifolder+2).name);
    if ~isfile(strcat(dirFile,'\transformInfo.mat'))
        save(strcat(dirFile,'\transformInfo.mat'),'transformInfo');
    end
    load([dirFile,'\transformInfo.mat'])
    existedCameraTformCheck(dirFile,channels_cam1,channels_cam2,transformInfo);
    saveas(gcf,[dirFile,filesep,'transformCheck.fig']);
    close all
    imageFluoReImageTransform_ys(transformInfo,dirFile);
    
end

%%
function existedCameraTformCheck(dirFile,channels_cam1,channels_cam2,transformInfo)
%check 现有的tform是否适用
% shuai yang 2020.05.28
dirField_check = strcat(dirFile,'\feild0001');
% if ismember(fixedChannel,channels_cam1)
%     chName_cam1= fixedChannel;
%     image_Name = strcat(dirField_check,'\',chName_cam1,'\image',chName_cam1,'00001.tif');
%     image_cam1 = import_tiff_stack(image_Name);
% else
% end
for iCh_cam1 = 2:numel(channels_cam1)
    chName_cam1 = channels_cam1{iCh_cam1};
    image_Name = strcat(dirField_check,'\',chName_cam1,'\image',chName_cam1,'00001.tif');
    try
        image_cam1 = import_tiff_stack(image_Name);
        break
    catch ME
        continue
    end
end
% 先试图import荧光图像，如果不存在荧光通道图像，再试图import明场图像
if ~exist('image_cam1','var')
    chName_cam1 = channels_cam1{1};
    image_Name = strcat(dirField_check,'\',chName_cam1,'\image',chName_cam1,'00001.tif');
    image_cam1 = import_tiff_stack(image_Name);
end

for iCh_cam2 = 1:numel(channels_cam2)
    chName_cam2 = channels_cam2{iCh_cam2};
    image_Name = strcat(dirField_check,'\',chName_cam2,'\image',chName_cam2,'00001.tif');
    try
        image_cam2 = import_tiff_stack(image_Name);
        break
    catch ME
        continue
    end
end

if ~isa(transformInfo.(chName_cam2),'double')
    tform = transformInfo.(chName_cam2);
%     image_cam2 = flip(image_cam2,2);
    image_cam2 = imwarp(image_cam2,tform,'OutputView',imref2d(size(image_cam2)));
else
    tform = transformInfo.(chName_cam1);
%     image_cam1 = flip(image_cam1,2);
    image_cam1 = imwarp(image_cam1,tform,'OutputView',imref2d(size(image_cam1)));
end
if strcmp(chName_cam1,'BF1')
    I0 = image_cam1;
    I0 = uint8(double(I0-min(I0(:)))/double(max(I0(:))-min(I0(:)))*255);
    I = imcomplement(I0);% equal to I = 255 - I0;
    I = imsubtract(I,I0);
    image_cam1 = I;
else
    I = image_cam1-100;
    image_cam1 = uint8(double(I-min(I(:)))/double(max(I(:))-min(I(:)))*255);
end

I = image_cam2-100;
image_cam2 = uint8(double(I-min(I(:)))/double(max(I(:))-min(I(:)))*255);

% C = imfuse(image_cam1,image_cam2,'falsecolor','Scaling','joint','ColorChannels',[2 1 0]);
C = cat(3,image_cam1,image_cam2,image_cam1*0);
figure, imshow(C,[])
end
%%
function imageFluoReImageTransform_ys(transformInfo,dirFile)
%此函数可以用imtransformFluoImage2FL_basedOnCam或者imtransformFluoImage2FL获得de
%transformInfo信息 对图像进行校正
% 根据旧函数imageFluoReImageTransform修改得到
%获得transformInfo后直接对用imwarp图像进行处理
%Apply geometric transformation to imagecollapse

disp('Apply geometric transformation of fluoImages.');
allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
% channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
% channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道

fieldList=dir(dirFile);
for iField = 1:length(fieldList)-2
    if ~strcmp(fieldList(iField+2).name(1:5),'feild')
        continue
    end
    dirField=strcat(dirFile,'\',fieldList(iField+2).name);
    disp(fieldList(iField+2).name)
    
    for iChannel = 1:numel(allChannels)

        dirImage = [dirField,'\',allChannels{iChannel}];
        % 利用try 判断是否存在此通道的image
        try
            load([dirImage,'\frameInfo.mat'])
        catch ME
            continue
        end
        tform = transformInfo.(allChannels{iChannel});
        %判断当前图像与fixedChannel是一个camera拍的 tform不计算
        if isa(tform,'double')%~isa(tform,'affine2d')
            continue
        end
        disp(strcat('Now transformation of',32, allChannels{iChannel},' images'));
        
        %                 dirNewImage=strcat(dirField,'\',channelList(iChannel).name,'_transform');
        %                 mkdir(dirNewImage);
        
        imageList = dir(dirImage);
        for iImage = 1:length(imageList)-2
            if ~strcmp(imageList(iImage+2).name(1:5),'image')
                continue
            end
            myImage = import_tiff_stack(strcat(dirImage,'\',imageList(iImage+2).name));
            
%             myImage = flip(myImage,2);%两个相机拍摄图像y轴镜像对称
            myImage = imwarp(myImage,tform,'OutputView',imref2d(size(myImage)));
            imwrite(myImage,strcat(dirImage,'\',imageList(iImage+2).name));
            %                         imwrite(myImage,strcat(dirNewImage,'\',imageList(iImage+2).name));
        end

    end
      
%  save(strcat(dirFile,'\transformInfo.mat'),'transformInfo');
end
end
