function  oneShotBatchImageMaskGet_ys(dirFile,fixedChannel)
%多通道荧光图像处理code Shuai Yang 2020.04.02
% fixedChannel = 'sfGFP';
% 适用于多视野oneShotmontage图像批处理的mask获得
% 针对oneshot实验
% 将每个视野第一张要做mask的图像拿出来 批处理生成mask 再存储
% Shuai Yang 2021.08.23

disp('oneShot Multifield image mask generation')
fieldList = dir([dirFile,filesep,'field*']);
% all field fixChannel Image Get
maskNameList  = cell([length(fieldList),1]) ;% 创建一个空的cell体 存储字符串
fixChTF = true(1,length(fieldList));% 用于判断视野中是否存在fixChannel
fieldNum = length(fieldList);

for iField = 1:fieldNum
    
    dirField = strcat(dirFile,'\',fieldList(iField).name);
    disp (fieldList(iField).name);
    
    dirImage = [dirField,filesep,fixedChannel];
    %
    if ~isfolder(dirImage)
        fixChTF(iField) = false;
        warning(['No images were captured in ',fixedChannel]);
        disp(['No images were captured in ',fixedChannel]);
        continue
    end
    
    dirTrack = [dirField,'\Tracking'];
    %如果已有tracking文件夹 删除里面的内容
    if isfolder(dirTrack)
        %                     rmdir (dirTrack, 's') %尝试删除 folderName 中的所有子文件夹和文件
        delete ([dirTrack,'\','*.*']);
        rmdir(dirTrack);
    end
    mkdir(dirTrack);
    
    load([dirImage,'\frameInfo.mat'])
    save([dirTrack,'\','frameInfo.mat'],'frameInfo');
    
    imageList = dir([dirImage,filesep,'image*.tif']);
    im_Seq = 1;% image frame sequence 可以取中间某个frame计算 如果不止一张图像

    imageStack(:,:,iField) = imread( strcat(dirImage,filesep, ...
        imageList(im_Seq).name) );

    %     fileName = [dirTrack,'\','imageTracking',fixedChannel,...
    %         imageList(im_Seq).name(end-8:end-4),'.mat'];


    %  oneShot 改为存储为tif 2022/5/25 Shuai Yang
    fileName = [dirTrack,'\','imageTracking',fixedChannel,...
        imageList(im_Seq).name(end-8:end-4),'.tif'];


    maskNameList{iField} = fileName;

end

% 以防field太多的montage 无法一次处理
if strcmp(fixedChannel,'BF1')
    stackGap = 64;
else
    stackGap = 200;
end

n = ceil(fieldNum/stackGap);%subImageStack的数目
for iSubStk = 1:n
    if iSubStk == n
        subImStack = imageStack(:,:,(iSubStk-1)*stackGap+1:end);
        submaskNameList = maskNameList((iSubStk-1)*stackGap+1:end);
    else
        subImStack = imageStack(:,:,(iSubStk-1)*stackGap+1:iSubStk*stackGap);
        submaskNameList = maskNameList((iSubStk-1)*stackGap+1:iSubStk*stackGap);
    end

    % maskImageGet
    switch fixedChannel
        case 'BF1'          
            [imageTrackings,~] = phCImProcessing_supperSegger(subImStack);% PhC图像处理 supperSegger 处理
        otherwise
            [imageTrackings] = fluoImCellSeg_ys(subImStack);% fluo 图像
    end

    % save maskImages
    for iIm = 1:size(subImStack,3)
        if ~fixChTF(iIm)
            continue
        end
        imageTracking = imageTrackings(:,:,iIm);
        %     save(maskNameList{iIm},'imageTracking');
        %     imwrite(uint8(imageTracking*255),[maskNameList{iIm}(1:end-4),'.tif']);
        imwrite(uint8(imageTracking*255),submaskNameList{iIm});
    end
    clear subImStack

end

end