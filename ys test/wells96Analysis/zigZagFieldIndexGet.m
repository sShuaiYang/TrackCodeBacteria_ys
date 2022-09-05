function [moveIndex]=zigZagFieldIndexGet(m,n)
% fieldNum=100;% n*n field n=2,4,6,8,10...
% n=sqrt(fieldNum);
% m 12 longth
% n 8 width
% m,n must be 
index=cell(1,n+1);
index{1}=1:m;
index{n+1}=[];
for i=n:-1:2
    if mod(i,2)==0
        index{n+1}(end+1)=i*m;
    else
        index{n+1}(end+1)=(i-1)*m+1;
    end    
end
for i=2:n
    if mod(i,2)==0
        index{i}=(i-1)*m+1:1:i*m-1;
    else
        index{i}=(i-1)*m+2:1:i*m;
    end
end
moveIndex=[];
for i=1:n+1
    moveIndex=[moveIndex,index{i}];
end
end