function  batchFluoImageMaskGet_ys(dirFile,fixedChannel,thresh)
%多通道荧光图像处理code Shuai Yang 2020.04.02
% fixedChannel = 'sfGFP';
% dirFile = 'F:\2020-04-22 promoterlibrary B217 test ys';
disp('Multifield fluo image mask generation')
fieldList = dir([dirFile,filesep,'field*']);
for iField = 1:length(fieldList)

    dirField = strcat(dirFile,'\',fieldList(iField).name);
    channelList = dir(dirField);
    disp (fieldList(iField).name);

    for iChannel = 1:length(channelList)
        if(isequal(channelList(iChannel).name,'.')||... % 去除系统自带的两个隐文件夹
                isequal(channelList(iChannel).name,'..')||...
                ~channelList(iChannel).isdir)                % 去除遍历中不是文件夹的
            continue;
        end

        if strcmp(channelList(iChannel).name,fixedChannel)
            dirImage = strcat(dirField,'\',channelList(iChannel).name);
            imageList = dir([dirImage,filesep,'image*.tif']);

            dirTrack = [dirField,'\Tracking'];
            %如果已有tracking文件夹 删除里面的内容
            if isfolder(dirTrack)
                %                     rmdir (dirTrack, 's') %尝试删除 folderName 中的所有子文件夹和文件
                delete ([dirTrack,'\','*.*']);
                rmdir(dirTrack);
            end
            mkdir(dirTrack);

            %%如果已有tiff2matlab文件夹 删除里面的内容
            dirMaskSave = [dirField,'\tiff2matlab'];
            %如果已有tracking文件夹 删除里面的内容
            if isfolder(dirMaskSave)
                %                     rmdir (dirMaskSave, 's') %尝试删除 folderName 中的所有子文件夹和文件
                delete ([dirMaskSave,'\','*.*']);
                rmdir(dirMaskSave);
            end

            load([dirImage,'\frameInfo.mat'])
            save([dirTrack,'\','frameInfo.mat'],'frameInfo');
            imageNum = length(imageList);
            if imageNum ~= 1 %不等于1 说明timelaspe 可能需要bioTree 生成
                dirMaskSave = strcat(dirField,'\tiff2matlab');
                mkdir(dirMaskSave);
            end

            stackGap = 500;% 500frame进行stack处理
            tic;
            n = 0;
            j = 0;
            for iImage = 1:imageNum
                n = n + 1;
                tempImages(:,:,n) = import_tiff_stack( strcat(dirImage,'\',imageList(iImage).name) );
                if n == stackGap||(j == floor(imageNum/stackGap) && n == mod(imageNum,stackGap))
                    disp (['imageStack',num2str((j+1),'%05.f')]);
                    switch fixedChannel
                        case 'BF1'
                            %                                 [imageTrackings] = phaseContrastImageProcessing_ys(tempImages);% PhC图像处理
                            [imageTrackings,~] = phCImProcessing_supperSegger(tempImages);% PhC图像处理 supperSegger 处理
                        otherwise
                            %                                 [imageTrackings] = fluoImageProcessing_ystest(tempImages,thresh);% fluo 图像
                            [imageTrackings] = fluoImCellSeg_ys(tempImages);% fluo 图像
                    end
                    % save single maskImages
                    for iframe = 1:n
                        imageTracking = imageTrackings(:,:,iframe);
                        %     filename = [dirTrack,'\','imageTracking',fixedChannel,imageList(j*stackGap+iframe).name(end-8:end-4),'.mat'];
                        %     save(filename,'imageTracking');
                        % 存储为8bit tif 图像 2022/6/8 shuai Yang
                        filename = [dirTrack,'\','imageTracking',fixedChannel,imageList(j*stackGap+iframe).name(end-8:end-4),'.tif'];
                        imwrite(uint8(imageTracking*255),filename);
                    end
                    j = j + 1;
                    n = 0;
                    if imageNum ~= 1
                        save([dirMaskSave,'\',num2str(j,'%05.f')],'imageTrackings','-v7.3');
                    end
                    clear imageTrackings imageTracking tempImages
                end

            end
        end
    end
end
end

