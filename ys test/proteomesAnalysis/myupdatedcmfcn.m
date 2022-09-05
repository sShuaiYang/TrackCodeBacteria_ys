function [ txt ] = myupdatedcmfcn(empt,event_obj,label)
%自定义datacursormode 使其显示自定义的标签
%ys 2020.08.28
%refer to:https://ww2.mathworks.cn/help/matlab/ref/datacursormode.html#bsawkea-7
%https://zhidao.baidu.com/question/809091143326523852.html

% pos = get(dcm_obj,'Position');
% txt = {['Time: ',num2str(pos(1))],...
% 	      ['Amplitude: ',num2str(pos(2))]};
dcm_obj = datacursormode(gcf);
info = getCursorInfo(dcm_obj);
ind = info.DataIndex;
txt = label{ind};
end

%% demo
% x=rand(1,10)*10;
% y=rand(1,10)*10;
% L={'A','A1','A2','A3','B','B1','B2','B3','C','C1'};
%  
% plot(x,y,'.');figure(gcf);
% dcm_obj = datacursormode(gcf);
% set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on');
% set(dcm_obj,'UpdateFcn',{@dcmcallback,L});
% function [ txt ] = dcmcallback(empt,event_obj,label)
% dcm_obj = datacursormode(gcf);
% info=getCursorInfo(dcm_obj);
% ind = info.DataIndex;
% txt=label{ind};
% end

