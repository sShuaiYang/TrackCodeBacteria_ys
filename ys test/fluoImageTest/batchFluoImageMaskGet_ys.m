function  batchFluoImageMaskGet_ys(dirFile,fixedChannel,thresh)
%��ͨ��ӫ��ͼ����code Shuai Yang 2020.04.02
% fixedChannel = 'sfGFP';
% dirFile = 'F:\2020-04-22 promoterlibrary B217 test ys';
disp('Multifield fluo image mask generation')
fieldList = dir([dirFile,filesep,'field*']);
for iField = 1:length(fieldList)

    dirField = strcat(dirFile,'\',fieldList(iField).name);
    channelList = dir(dirField);
    disp (fieldList(iField).name);

    for iChannel = 1:length(channelList)
        if(isequal(channelList(iChannel).name,'.')||... % ȥ��ϵͳ�Դ����������ļ���
                isequal(channelList(iChannel).name,'..')||...
                ~channelList(iChannel).isdir)                % ȥ�������в����ļ��е�
            continue;
        end

        if strcmp(channelList(iChannel).name,fixedChannel)
            dirImage = strcat(dirField,'\',channelList(iChannel).name);
            imageList = dir([dirImage,filesep,'image*.tif']);

            dirTrack = [dirField,'\Tracking'];
            %�������tracking�ļ��� ɾ�����������
            if isfolder(dirTrack)
                %                     rmdir (dirTrack, 's') %����ɾ�� folderName �е��������ļ��к��ļ�
                delete ([dirTrack,'\','*.*']);
                rmdir(dirTrack);
            end
            mkdir(dirTrack);

            %%�������tiff2matlab�ļ��� ɾ�����������
            dirMaskSave = [dirField,'\tiff2matlab'];
            %�������tracking�ļ��� ɾ�����������
            if isfolder(dirMaskSave)
                %                     rmdir (dirMaskSave, 's') %����ɾ�� folderName �е��������ļ��к��ļ�
                delete ([dirMaskSave,'\','*.*']);
                rmdir(dirMaskSave);
            end

            load([dirImage,'\frameInfo.mat'])
            save([dirTrack,'\','frameInfo.mat'],'frameInfo');
            imageNum = length(imageList);
            if imageNum ~= 1 %������1 ˵��timelaspe ������ҪbioTree ����
                dirMaskSave = strcat(dirField,'\tiff2matlab');
                mkdir(dirMaskSave);
            end

            stackGap = 500;% 500frame����stack����
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
                            %                                 [imageTrackings] = phaseContrastImageProcessing_ys(tempImages);% PhCͼ����
                            [imageTrackings,~] = phCImProcessing_supperSegger(tempImages);% PhCͼ���� supperSegger ����
                        otherwise
                            %                                 [imageTrackings] = fluoImageProcessing_ystest(tempImages,thresh);% fluo ͼ��
                            [imageTrackings] = fluoImCellSeg_ys(tempImages);% fluo ͼ��
                    end
                    % save single maskImages
                    for iframe = 1:n
                        imageTracking = imageTrackings(:,:,iframe);
                        %     filename = [dirTrack,'\','imageTracking',fixedChannel,imageList(j*stackGap+iframe).name(end-8:end-4),'.mat'];
                        %     save(filename,'imageTracking');
                        % �洢Ϊ8bit tif ͼ�� 2022/6/8 shuai Yang
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

