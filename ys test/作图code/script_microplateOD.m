

% PAO1(:,1)=amean_1211(:,2)-amean_1211(:,1);
% PAO1(:,2)=amean_1211(:,8)-amean_1211(:,7);
% PAO1(1:110,3)=amean_1212(1:110,2)-amean_1212(1:110,1);
% PAO1(1:110,4)=amean_1212(1:110,8)-amean_1212(1:110,7);
% 
% YZ(:,1)=amean_1211(:,3)-amean_1211(:,1);
% YZ(:,2)=amean_1211(:,9)-amean_1211(:,7);
% YZ(1:110,3)=amean_1212(1:110,3)-amean_1212(1:110,1);
% YZ(1:110,4)=amean_1212(1:110,9)-amean_1212(1:110,7);
% 
% RsmA(:,1)=amean_1211(:,4)-amean_1211(:,1);
% RsmA(:,2)=amean_1211(:,10)-amean_1211(:,7);
% RsmA(1:110,3)=amean_1212(1:110,4)-amean_1212(1:110,1);
% RsmA(1:110,4)=amean_1212(1:110,10)-amean_1212(1:110,7);
% 
% GacA(:,1)=amean_1211(:,5)-amean_1211(:,1);
% GacA(:,2)=amean_1211(:,11)-amean_1211(:,7);
% GacA(1:110,3)=amean_1212(1:110,5)-amean_1212(1:110,1);
% GacA(1:110,4)=amean_1212(1:110,11)-amean_1212(1:110,7);
% 
% RetS(:,1)=amean_1211(:,6)-amean_1211(:,1);
% RetS(:,2)=amean_1211(:,12)-amean_1211(:,7);
% RetS(1:110,3)=amean_1212(1:110,6)-amean_1212(1:110,1);
% RetS(1:110,4)=amean_1212(1:110,12)-amean_1212(1:110,7);

%取1-4h进行拟合

figure,line(F,log2(PAO1),'lineWidth',2,'Color', c(2,:))
for i=1:4
    f=fit(F(7:27)',log2(PAO1(7:27,i)),'poly1');
    DT(1,i)=1/f.p1;
    hold on
    plot(f,F(7:27)',log2(PAO1(7:27,i)))
    
end

figure,line(F,log2(YZ),'lineWidth',2,'Color', c(3,:))
for i=1:4
    f=fit(F(7:27)',log2(YZ(7:27,i)),'poly1');
    DT(2,i)=1/f.p1;
    hold on
    plot(f,F(7:27)',log2(YZ(7:27,i)))
    
end

figure,line(F,log2(RsmA),'lineWidth',2,'Color', c(4,:))
for i=1:4
    f=fit(F(7:27)',log2(RsmA(7:27,i)),'poly1');
    DT(3,i)=1/f.p1;
    hold on
    plot(f,F(7:27)',log2(RsmA(7:27,i)))
    
end

figure,line(F,log2(GacA),'lineWidth',2,'Color', c(5,:))
for i=1:4
    f=fit(F(7:27)',log2(GacA(7:27,i)),'poly1');
    DT(4,i)=1/f.p1;
    hold on
    plot(f,F(7:27)',log2(GacA(7:27,i)))
    
end

figure,line(F,log2(RetS),'lineWidth',2,'Color', c(6,:))
for i=1:4
    f=fit(F(7:27)',log2(RetS(7:27,i)),'poly1');
    DT(5,i)=1/f.p1;
    hold on
    plot(f,F(7:27)',log2(RetS(7:27,i)))
    
end

figure,
for i=1:5
    plot(F,log2(amean1(:,i)),'lineWidth',2,'Color', c(i+1,:));
    hold on
    f=fit(F(13:43)',log2(amean1(13:43,i)),'poly1');
    DT(1,i)=1/f.p1;
    hold on
    plot(f,F(13:43)',log2(amean1(13:43,i)))
    
end

