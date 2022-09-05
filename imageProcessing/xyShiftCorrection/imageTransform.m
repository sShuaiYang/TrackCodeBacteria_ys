function image=imageTransform()
dirFile=uigetdir();
nameList=dir(dirFile);
se=translate(strel(1),[18,58]);
parfor i=1:numel(nameList)-2
    image(:,:,i)=import_tiff_stack([dirFile,'\',nameList(i+2).name]);
    image(:,:,i)=imdilate(image(:,:,i),se);
end
end