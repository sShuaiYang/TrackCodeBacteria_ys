function [dataCollect]= timeLapseMultiFieldFluoImageProtocal(dirFile,fixedChannel)
% shuai Yang 2020.05.15
% fixedChannel,example: 'sfGFP';
% %1.荧光图像视场不均匀校正
% fluoImageViewNonuniformCorrection(dirFile)
% 2.相机的校正信息的获得以及新图像的生成
% cameraCorrectionForMultiFluoChannels(dirFile,fixedChannel);
% 3.time lapse xy shift cooretion
% timeLapsexyShiftCorrectionForMultiFluoChannels(dirFile,fixedChannel);
% [~] = batchXYShiftCorrectionUsingFTPhaseCorrelation(dirFile,fixedChannel);
% 4.根据fixedChannel生成Tracking信息
% 如果需要自定义阈值 请在[]中填入数值即可 如 [5]
batchImageMaskGet_ys(dirFile,fixedChannel,[]);
% 5.生成bioTree
multiFieldFluoImageBioTreeGeneration(dirFile)
% 6.获得bioInfo信息和添加荧光背景信号
multiFieldFluoImageBioInfoGetByBioTree(dirFile);
addingBKGSignal2BioInfo(dirFile);
% 7.细菌绝对位置的获得 unit um 以及荧光光强校正
timeLapseMultiFieldAbsCentroidGetAndIntCorrection(dirFile);
% 8.bioInfo reshape 按时间点存储
timeLapseMultiFieldBioInfoReshape(dirFile)
% 9.plot function; data visualization
[~] = timeLapseMultiFieldFluoImagePlot(dirFile,fixedChannel);
% 用校正后的荧光强度作图
try 
    [dataCollect] = timeLapseMultiFieldFluoImagePlot_2(dirFile,fixedChannel);
catch ME
    rethrow(ME)
end
end
