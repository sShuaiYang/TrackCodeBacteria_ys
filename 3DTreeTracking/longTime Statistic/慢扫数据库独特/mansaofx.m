function mansaofx(Xrange)
for i=2:62
    group1=[1,i];
    figure(i)
    subplot(2,3,1);name=mansao(group1,'beforeDivision',Xrange);
    subplot(2,3,2);name=mansao(group1,'afterDivision',Xrange);
    subplot(2,3,3);name=mansao(group1,'lengthTime',Xrange);
    subplot(2,3,4);name=mansao(group1,'divisionTime',Xrange);
    subplot(2,3,5);name=mansao(group1,'detachingResult',Xrange);
    subplot(2,3,6);name=mansao(group1,'velocityResult',Xrange);
    dirFile='C:\Users\user\Desktop\2015.8.11ÂýÉ¨\picture';
    name1=strcat(dirFile,'\',name{2},'.tiff');
    saveas(i,name1,'tiffn');
    clear name
end
end
 
 
 
 