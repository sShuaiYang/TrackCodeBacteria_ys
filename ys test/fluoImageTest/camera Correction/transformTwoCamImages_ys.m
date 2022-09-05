function transformTwoCamImages_ys (transformInfo,dirFile)
% 获得 transformInfo信息 对图像进行校正
% 根据旧函数imageFluoReImageTransform修改得到
%获得transformInfo后直接对用imwarp图像进行处理
%Apply geometric transformation to image 

disp('Apply geometric transformation of two cam Images.');
% allChannels = {'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1拍摄的通道 短波长
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2拍摄的通道

fieldList = dir([dirFile,filesep,'field*']);
for iField = 1:length(fieldList)
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    disp(fieldList(iField).name)

    for iCh_cam2 = 1:numel(channels_cam2)

        dirImage = [dirField,'\',channels_cam2{iCh_cam2}];
        % 利用try 判断是否存在此通道的image
        try
            load([dirImage,'\frameInfo.mat'])
        catch ME
            continue
        end

        tform = transformInfo.(channels_cam2{iCh_cam2});
        %判断当前图像与fixedChannel是一个camera拍的 tform不计算
        if isa(tform,'double') %~isa(tform,'affine2d')
            continue
        end
        disp(strcat('Now transformation of',32, channels_cam2{iCh_cam2},' images'));
        %  dirNewImage=strcat(dirField,'\',channels_cam2{iCh_cam2}.name,'_transform');
        %  mkdir(dirNewImage);

        imageList = dir([dirImage,filesep,'image*.tif']);
        for iImage = 1:length(imageList)
            myImage = imread(strcat(dirImage,'\',imageList(iImage).name));
            myImage = flip(myImage,1);%两个相机拍摄图像y轴镜像对称
            myImage = imwarp(myImage,tform,'OutputView',imref2d(size(myImage)));
            if ~exist('imPixelsTF','var')
                % 获取赋值为0得图像得pixel的位置 用于将cam1相应的区域赋值为0
                imPixelsTF = double(myImage)== 0;
            end
            imwrite(myImage,strcat(dirImage,'\',imageList(iImage).name));
            %  imwrite(myImage,strcat(dirNewImage,'\',imageList(iImage).name));
        end
    end

    % 将cam2 赋值为0图像区域，也给cam1的图像赋值为0;PhC的BF1通道赋值中值
    disp('Now specfied values for outside pixels in Cam1 images')
    disp('Fill MEDIAN values for PhC images')
    disp('Fill ZERO values for fluo images')

    for iCh_cam1 = 1:numel(channels_cam1)

        dirImage = [dirField,'\',channels_cam1{iCh_cam1}];
        % 利用try 判断是否存在此通道的image
        try
            load([dirImage,'\frameInfo.mat'])
        catch ME
            continue
        end

        tform = transformInfo.(channels_cam1{iCh_cam1});
        % cam1 的tform为0
        if ~isa(tform,'double')
            continue
        end
        imageList = dir([dirImage,filesep,'image*.tif']);
        disp(['Cam1 Channel: ',channels_cam1{iCh_cam1}]);

        if strcmp(channels_cam1{iCh_cam1},'BF1')            
            for iImage = 1:length(imageList)
                myImage = imread(strcat(dirImage,'\',imageList(iImage).name));
                myImage(imPixelsTF) = median(myImage(:));
                imwrite(myImage,strcat(dirImage,'\',imageList(iImage).name));
            end
        else            
            for iImage = 1:length(imageList)
                myImage = imread(strcat(dirImage,'\',imageList(iImage).name));
                myImage(imPixelsTF) = 0;
                imwrite(myImage,strcat(dirImage,'\',imageList(iImage).name));
            end
        end
    end

end