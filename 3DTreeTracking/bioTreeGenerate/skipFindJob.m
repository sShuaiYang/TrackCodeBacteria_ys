function [Result1,Result2]=skipFindJob(velocityData)
                   
Result1(:,1)=find(velocityData(:,4)>1);
Result1(:,2)=velocityData(Result1(:,1),1);
Result1(:,3)=velocityData(Result1(:,1),4);
Result1(:,4)=velocityData(Result1(:,1),9);

Result2(:,1)=find(velocityData(:,9)>1);
Result2(:,2)=velocityData(Result2(:,1),1);
Result2(:,3)=velocityData(Result2(:,1),4);
Result2(:,4)=velocityData(Result2(:,1),9);


end