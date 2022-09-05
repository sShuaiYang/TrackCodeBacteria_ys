dirFile='E:\2015-01-10 jzy & xag xyScan';
dirAll='E:\2015-01-10 jzy & xag xyScan\allResult';
dirSmallFile='E:\2015-01-10 jzy & xag xyScan\allResult\m0000';
dirMaskImagePre='E:\2015-01-10 jzy & xag xyScan\maskImage';
dirBioTreeResult='E:\2015-01-10 jzy & xag xyScan\bioTreeResult';
for i=1:399
    dirSmallFileNew=dirSmallFile;
    dirSmallFileNew(end-numel(num2str(i-1))+1:end)=num2str(i-1);
    mkdir(dirSmallFileNew)
    cd(dirSmallFileNew)
    mkdir('bioTreeResult')
    mkdir('maskImage')
    dirResult1=strcat(dirSmallFileNew,'\bioTreeResult\1');
    dirMaskImage1=strcat(dirSmallFileNew,'\maskImage\1');
    load(strcat(dirBioTreeResult,'\',num2str(i)))
    save(dirResult1,'bioTree');
    load(strcat(dirMaskImagePre,'\',num2str(i)))
    save(dirMaskImage1,'imageGFP','imageRFP','maskImage')
end