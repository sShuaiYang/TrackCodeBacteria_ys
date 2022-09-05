function imageFluoReImageTransform_ysOld(transformInfo,dirFile,fixedChannel)
%�˺���������imtransformFluoImage2FL_basedOnCam����imtransformFluoImage2FL���de
%transformInfo��Ϣ ��ͼ�����У��
% ���ݾɺ���imageFluoReImageTransform�޸ĵõ�
%���transformInfo��ֱ�Ӷ���imwarpͼ����д���
%Apply geometric transformation to imagecollapse

disp('Apply geometric transformation of fluoImages.');
allChannels={'BF1','PVD','Venus','sfGFP','mScarletI','CyOFP','TDsmURFP'};
channels_cam1 = {'BF1','PVD','sfGFP'};% camera 1�����ͨ�� �̲���
channels_cam2 = {'mScarletI','CyOFP','Venus','TDsmURFP'};% camera 2�����ͨ��

if ismember(fixedChannel,channels_cam1)
    fixedChanLib = channels_cam1;%�̶���channels library
else
    fixedChanLib = channels_cam2;
end

fieldList=dir(dirFile);
for iField = 1:length(fieldList)-2
    if strcmp(fieldList(iField+2).name(1:5),'field')
        dirField=strcat(dirFile,'\',fieldList(iField+2).name);
        disp(fieldList(iField+2).name)
        channelList=dir(dirField);
        for iChannel=1:length(channelList) 
            if(isequal(channelList(iChannel).name,'.')||... % ȥ��ϵͳ�Դ����������ļ���
                    isequal(channelList(iChannel).name,'..')||...
                    ~channelList(iChannel).isdir)                % ȥ�������в����ļ��е�
                continue;
            end
            if ismember(channelList(iChannel).name,allChannels ) && ~strcmp(channelList(iChannel).name,fixedChannel)
                fluoChannel = channelList(iChannel).name;
                disp(strcat('Now transformation of',32,fluoChannel,' images'));
                dirImage = strcat(dirField,'\',channelList(iChannel).name);
%                 dirNewImage=strcat(dirField,'\',channelList(iChannel).name,'_transform');
%                 mkdir(dirNewImage);
                imageList = dir(dirImage);
                for iImage = 1:length(imageList)-2
                    if strcmp(imageList(iImage+2).name(1:5),'image')
                        imageFluo=import_tiff_stack(strcat(dirImage,'\',imageList(iImage+2).name));
                        switch fluoChannel
                            case 'BF1'
                                tform = transformInfo.BF1;
                                %�жϵ�ǰfixedChannel��һ��camera�ĵ� tform������
                                if isa(tform,'double')%~isa(tform,'affine2d')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                            case 'mScarletI'
                                tform = transformInfo.mScarletI;
                                if isa(tform,'double')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                            case 'sfGFP'
                                tform = transformInfo.sfGFP;
                                if isa(tform,'double')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                            case 'CyOFP'
                                tform = transformInfo.CyOFP;
                                if isa(tform,'double')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                            case 'Venus'
                                tform = transformInfo.Venus;
                                if isa(tform,'double')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                            case 'PVD'
                                tform = transformInfo.PVD;
                                if isa(tform,'double')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                            case'TDsmURFP'
                                tform = transformInfo.TDsmURFP;
                                if isa(tform,'double')
                                    break
                                end
                                if ~ismember(fluoChannel,fixedChanLib)
                                    imageFluo = flip(imageFluo,2);%�����������ͼ��y�᾵��Գ�
                                end
                                imageFluo = imwarp(imageFluo,tform,'OutputView',imref2d(size(imageFluo)));
                        end
                        imwrite(imageFluo,strcat(dirImage,'\',imageList(iImage+2).name));
                        %                         imwrite(imageFluo,strcat(dirNewImage,'\',imageList(iImage+2).name));
                    end
                end
            end
        end
    end
    
 save(strcat(dirFile,'\transformInfo.mat'),'transformInfo');   
end