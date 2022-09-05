function sumStatisticMap(Xrange)
for i=2:62
    group1=[1,i];
    figure(i)
    subplot(2,3,1);name=mansao(group1,'beforeDivision',Xrange);
    subplot(2,3,2);name=mansao(group1,'afterDivision',Xrange);
    subplot(2,3,3);name=mansao(group1,'lengthTime',Xrange);
    subplot(2,3,4);name=mansao(group1,'divisionTime',Xrange);
    subplot(2,3,5);name=mansao(group1,'detachingResult',Xrange);
    subplot(2,3,6);name=mansao(group1,'velocityResult',Xrange);
    dirFile='C:\Users\user\Desktop\2015.8.11慢扫\picture';
    name1=strcat(dirFile,'\',name{2},'.tiff');
    saveas(i,name1,'tiffn');
    clear name
end
end
function name1=mansao(group1,choice,Xrange)
% 对比慢扫数据，group1是选中要比较的mat文件，Xrange是在afterDivision或者beforeDivision设置X轴区间，choice是选择要比较的种类
[xdata,ydata,unit1,unit2,name1]=getXY(group1,choice,Xrange);
if strcmp('afterDivision',choice)==1
    [xdata,ydata,fengX]=Dhist(ydata,Xrange);
    disp(fengX)
end
if strcmp('beforeDivision',choice)==1
    [xdata,ydata,fengX]=Dhist(ydata,Xrange);
    disp(fengX)
end

hold on
colormap=jet(size(group1,2));
for i=1:size(group1,2)
    if size(xdata{i},1)==1&&size(ydata{i},1)==1
        v=[xdata{i};ydata{i}];
    else
        v1=xdata{i}';v2=ydata{i}';
        v=[v1;v2];
    end
    v=spcrv(v,3);
    %       plot(xdata{i},ydata{i},'Color',[0.3+1/(2*i),0.1*i,1-0.1*i],'LineWidth',0.5+0.2*i);
    plot(v(1,:),v(2,:),'Color',colormap(i,:),'LineWidth',1);
    str{i}=name1{i};
    clear v
end
xlabel(unit1);ylabel(unit2);title(choice);
legend(str)
hold off
end

function [X,Y,fengzhiX]=Dhist(dlength,Xrange)
%对before或者afterDivision数据进行统计
linSpaceX=linspace(Xrange(1),Xrange(2),101);
number=size(dlength,2);
for inum=1:number   
    for i=2:2:size(linSpaceX,2)
        Y{inum}(i/2)=sum(dlength{inum}>=linSpaceX(i-1)&dlength{inum}<linSpaceX(i+1));
    end
    X{inum}=linSpaceX(2:2:end);
    Y{inum}=Y{inum}/sum(Y{inum});
    ymax=find(Y{inum}==max(Y{inum}));
    fengzhiX(inum)=X{inum}(ymax(1));
end
end
