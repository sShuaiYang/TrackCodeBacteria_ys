function [dataCollect] = oneShotMontageFluoImageProtocal(dirFile,fixedChannel)
% shuai Yang 2020.05.15
%fixedChannel,example: 'sfGFP';
%针对只拍摄一次的montage图像的处理，统计视野每个通道细菌荧光强度，菌长等基本信息
% %1.荧光图像视场不均匀校正
% fluoImageViewNonuniformCorrection(dirFile);
%2.相机的校正信息的获得以及新图像的生成
% cameraCorrectionForMultiFluoChannels(dirFile,fixedChannel);
%3.如果需要自定义阈值 请在[]中填入数值即可 如 [5]
% oneShotBatchImageMaskGet_ys(dirFile,fixedChannel);
%4.获得bioInfo信息
[bioInfo] = oneShotMontageFluoImageBioInfoGet(dirFile);
% %5.细菌绝对位置的获得 unit um
% [bioInfo] = absCentroidGetInMultiFields(bioInfo,dirFile);
% %6.荧光光强校正
[bioInfo] = cellFluoIntensityCorrection(bioInfo,dirFile);
%7.plot function; data visualization
[~] = oneShotMontageFluoImagePlot(dirFile,bioInfo,fixedChannel);
[dataCollect] = oneShotMontageFluoImagePlot_2(dirFile,bioInfo,fixedChannel);
end

