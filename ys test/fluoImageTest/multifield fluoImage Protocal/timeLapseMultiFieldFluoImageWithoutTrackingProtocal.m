function [allDataCollect]= timeLapseMultiFieldFluoImageWithoutTrackingProtocal(dirFile,fixedChannel)
% shuai Yang 2020.05.15
% fixedChannel,example: 'sfGFP';
% 针对只拍摄一次的montage图像的处理，统计视野每个通道细菌荧光强度，菌长等基本信息
% %1.荧光图像视场不均匀校正
% fluoImageViewNonuniformCorrection(dirFile)
% 2.相机的校正信息的获得以及新图像的生成
cameraCorrectionForMultiFluoChannels(dirFile,fixedChannel);
% 3.根据fixedChannel生成Tracking信息 
% 如果需要自定义阈值 请在[]中填入数值即可 如 [5]
batchImageMaskGet_ys(dirFile,fixedChannel,[]);
% 4.获得bioInfo信息
[~] = multiFieldFluoImageBioInfoGetWithoutTracking(dirFile);
% 5.细菌绝对位置的获得 unit um 以及荧光光强校正
timeLapseMultiFieldAbsCentroidGetAndIntCorrection(dirFile);
% 6.bioInfo reshape 按时间点存储
timeLapseMultiFieldBioInfoReshape(dirFile);
% 7.plot function; data visualization
% 加入了计算光强随时间的变化
[~] = timeLapseMultiFieldFluoImageWithoutTrackingPlot(dirFile,fixedChannel);
% 每个视野都计算了荧光强度的变化
% timeLapseMultiFieldFluoImagePlot_singleField(dirFile,fixedChannel);
% 用校正后的荧光强度作图
try 
    [allDataCollect] = timeLapseMultiFieldFluoImageWithoutTrackingPlot_2(dirFile,fixedChannel);
%     timeLapseMultiFieldFluoImagePlot_singleField_2(dirFile,fixedChannel);
catch ME
    rethrow(ME)
end
end
