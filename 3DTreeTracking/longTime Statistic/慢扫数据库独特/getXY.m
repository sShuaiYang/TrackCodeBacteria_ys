function [xdata,ydata,unit1,unit2,name1]=getXY(group1,choice,Xrange)
% 获得mat文件里的数据
dirallResult='E:\test';
for i=1:size(group1,2)
   dirFile=strcat(dirallResult,'\newallResult\allResult',num2str(group1(i)),'.mat'); 
   load(dirFile);
   if isfield(allResult,choice)==0
       eval(['allResult.',choice,'=zeros(100,2)']);
   end
   switch choice
    case 'afterDivision'
        result=[allResult.afterDivision];unit1='length（um）';unit2=' ';
       case 'beforeDivision'
           result=[allResult.beforeDivision];unit1='length(um)';unit2=' ';
       case 'divisionTime'
           result=[allResult.divisionTime];unit1='time(min)';unit2='time(min)';
       case 'velocityResult'
           result=[allResult.velocityResult];unit1='time(min)';unit2='velocityrate(um/min)';
       case 'lengthTime'
           result=[allResult.lengthTime];unit1='time(min)';unit2=' length(um)';
       case 'detachingResult'
           result=[allResult.detachingResult];unit1='bacterianum';unit2='detachaingrate ';
   end
   name1{i}=allResult.name;
   if size(result,2)==1
       linSpaceX=linspace(Xrange(1),Xrange(2),101);
       xdata{i}=linSpaceX(2:2:end);
       ydata{i}=result;
   else
       xdata{i}=result(:,1);
       ydata{i}=result(:,2);
   end
   
end
end

