function fluoColonyImageResultGetProtocal(dirFile,fixedChannel)
%体式显微镜数据处理
%Shuai Yang 2021.03.29
%1.mask 获得
batchFluoColonyImageMaskGet_ys(dirFile,fixedChannel,[]);
%2.获得bioInfo信息
[~] = fluoColonyImageBioInfoGet(dirFile);
end