function makeMovie(dirFile,lowIn,highIn,fluoChannel)
fluoFile=[dirFile,'\',fluoChannel,'new'];
movieFile=[dirFile,'\',fluoChannel,'movie'];
mkdir(movieFile)
nameList=dir(fluoFile);
fileNum=numel(nameList)-3;
for i=1:fileNum
    load([fluoFile,'\',nameList(i+3).name])
    imageI=uint8(255*(imageI-lowIn)/(highIn-lowIn));
    imwrite(imageI,[movieFile,'\',nameList(i+3).name(1:end-4),'.tif'])
end
end