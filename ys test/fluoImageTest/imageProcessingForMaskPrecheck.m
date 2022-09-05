function  [imageTrackings,imIdx] = imageProcessingForMaskPrecheck(dirFile,fixedChannel,thresh)
%thresh =[];or other numerical values 
%example: [imageTrackings,imIdx] = imageProcessingForMaskPrecheck(dirFile,fixedChannel,[])
%用于图像处理时Mask生成的预先检查code Shuai Yang 2020.08.07
%fixedChannel 不局限于荧光 如果是BF 需要是PHC 图像
% fixedChannel = 'sfGFP';
disp('image mask generation for precheck')
fieldList = dir(dirFile);

iField = 1;%数字代表取第几个视野进行check

if strcmp(fieldList(iField+2).name(1:5),'field')
    dirField = strcat(dirFile,'\',fieldList(iField+2).name);
    channelList = dir(dirField);
    disp (fieldList(iField+2).name);
    
    for iChannel = 1:length(channelList)
        if(isequal(channelList(iChannel).name,'.')||... % 去除系统自带的两个隐文件夹
                isequal(channelList(iChannel).name,'..')||...
                ~channelList(iChannel).isdir)                % 去除遍历中不是文件夹的
            continue;
        end
        
        if strcmp(channelList(iChannel).name,fixedChannel)
            dirImage = strcat(dirField,'\',channelList(iChannel).name);
            imageList = dir(dirImage);
            imageNum = length(imageList)-3;
            
            if imageNum <= 400 
                imIdx = 1:imageNum;
            else
                imIdx = 1:floor(imageNum/200):imageNum;
            end
            % image 类型的获得 例如'uint16'
            tempImage = import_tiff_stack( strcat(dirImage,'\',imageList(4).name) );
            imageStack = zeros(size(tempImage,1),size(tempImage,2),numel(imIdx),class(tempImage));
            
            parfor iImage = 1:numel(imIdx)
                imageStack(:,:,iImage) = import_tiff_stack( strcat(dirImage,'\',imageList(imIdx(iImage)+3).name) );
            end
            
            switch fixedChannel
                case 'BF1'
                    [imageTrackings] = phaseContrastImageProcessing_ys(imageStack);% PhC图像处理
                otherwise
                    [imageTrackings] = fluoImageProcessing_ystest(imageStack,thresh);% fluo 图像
            end
            
        end
    end
end
end

