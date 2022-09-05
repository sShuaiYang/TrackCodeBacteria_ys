function RgVsTime=getAveRg(imageStack,scale)
for iStack=1:size(imageStack,3);
    [xCo,yCo]=find(imageStack(:,:,iStack)==1);
%     meanCo=[1136.8 1326.3];
    if iStack==1
        meanCo=[mean(xCo),mean(yCo)];
    end
    allCoo=cat(2,xCo,yCo);
    RgVsTime(iStack,2)=mean(pdist2(allCoo,meanCo).^2).^0.5*scale;
    RgVsTime(iStack,3)=numel(xCo);
end
RgVsTime(:,1)=0:10:10*53;
RgVsTime(:,4)=(RgVsTime(:,3).*scale.^2)/pi./(RgVsTime(:,2).^2);
end
%4.24 test
% 40x �ﾵ��0:10:220��0.160256
% 20x �ﾵ 240:10:(240+10*70) 0.325380

% 40x ��ʼ�� [1093.6 1268.2]  ���һ��ƽ�� [1342.1,1218.7]

% 20x ��һ��ƽ�� [1089.5,1293.6] ������ʼ����[965.25,1318.3]

%5.18 F6
% 40x  0:10:33*10, 0.160256   ��ʼ�� [1082.1 1306.3]  ���һ��ƽ��[1021.7 1054.5]
% 20x  340:10:(340+10*79) 1.32538  ��һ��ƽ�� [1106.6 1200.4] ������ʼ����
% [1106.6+��1082.1-1021.7��/2,1200.4+(1306.3-1054.5)/2] =[1136.8 1326.3]