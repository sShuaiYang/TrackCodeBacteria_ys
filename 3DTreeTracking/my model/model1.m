function bacteria=model1(num)
% this function is a simulation for bacteria survival strategy
% some basic hypothesis
% survival for getting contact with more bacteria
% the ability of motility is inversely proportional to the degree
% a. bacteria with low degree should have a high probability for detaching
%    and they would be washed away by flow
% b. the probability to know a new one is direct proportional to the
%    degree, and the link is correlated with both node
% c. know each other & face to detaching
tao1=num/100;
tao2=10;
possibility=0.0002;
basicDetachingRate=0.1;
for i=1:num
    bacteria(i).linkInfo=[];
    bacteria(i).linkNum=0;
    bacteria(i).linkRate=1;
    bacteria(i).DetachingRate=[];
end
bacteria=connectEach(bacteria,tao1);
bacteria=bacteriaCopy(bacteria,possibility,tao1);
bacteria=detachingEach(bacteria,basicDetachingRate,tao2);
end
function bacteria=connectEach(bacteria,tao1)
num=numel(bacteria);
for i=1:num/10
    numI=fix(rand*num);
    numJ=fix(rand*num);
    if numI==0 || numJ==0 || numI==numJ
        i=i-1;
        continue
    end
    linkRate=bacteria(numI).linkRate*bacteria(numJ).linkRate;
    if rand<linkRate
        bacteria(numI).linkInfo=[bacteria(numI).linkInfo;numJ];
        bacteria(numI).linkNum=bacteria(numI).linkNum+1;
        bacteria(numI).linkRate=0.2+0.8*exp(-bacteria(numI).linkNum/tao1);
        bacteria(numJ).linkInfo=[bacteria(numJ).linkInfo;numI];
        bacteria(numJ).linkNum=bacteria(numJ).linkNum+1;
        bacteria(numJ).linkRate=0.2+0.8*exp(-bacteria(numJ).linkNum/tao1);
    end
end
end
function bacteria=detachingEach(bacteria,basicDetachingRate,tao2)
num=numel(bacteria);
for i=1:num
    bacteria(i).DetachingRate=basicDetachingRate*exp(-bacteria(i).linkNum/tao2)+0.01;
    if rand<bacteria(i).DetachingRate
        bacteria(i).linkInfo=[];
        bacteria(i).linkNum=0;
        bacteria(i).linkRate=1;
        bacteria(i).DetachingRate=basicDetachingRate+0.01;
    end
end
end
function regionHist=getDegreeDistri(bacteria)
for i=1:numel(bacteria)
    degree(i)=bacteria(i).linkNum;
end
maxDegree=max(degree);
regionHist(1,:)=1:maxDegree;
for i=1:maxDegree
    regionHist(2,i)=numel(degree(degree>=i));
end
end
function bacteria=bacteriaCopy(bacteria,possibility,tao1)
for i=1:numel(bacteria)
    if bacteria(i).linkNum<=tao1
        bacteria(i).copyRate=possibility*(1/2+exp(-(tao1-bacteria(i).linkNum)/tao1));
    else
        bacteria(i).copyRate=1.5*possibility;
    end
    if rand<bacteria(i).copyRate
        bacteria(end+1)=bacteria(i);
        bacteria(i).linkInfo=[bacteria(i).linkInfo;numel(bacteria)];
        bacteria(i).linkNum=bacteria(i).linkNum+1;
        bacteria(i).linkRate=0.2+0.8*exp(-bacteria(i).linkNum/tao1);
        bacteria(end).linkInfo=[bacteria(end).linkInfo;i];
        bacteria(end).linkNum=bacteria(end).linkNum+1;
        bacteria(end).linkRate=0.2+0.8*exp(-bacteria(end).linkNum/tao1);
    end
end
end