for i=1:size(result.head,2)
    meanVelocity.low(i,1)=sum(result.head{i}.low(:,2))/sum(result.head{i}.low(:,1));
    meanVelocity.high(i,1)=sum(result.head{i}.high(:,2))/sum(result.head{i}.high(:,1));
    meanVelocity.all(i,1)=(sum(result.head{i}.low(:,2))+sum(result.head{i}.high(:,2)))/(sum(result.head{i}.low(:,1))+sum(result.head{i}.high(:,1)));
end
meanVelocity.lowvar=var(meanVelocity.low(i,1));
meanVelocity.highvar=var(meanVelocity.high(i,1));
meanVelocity.allvar=var(meanVelocity.all(i,1));