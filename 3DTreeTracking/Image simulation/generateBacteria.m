function model=generateBacteria(length,oritation)
% if length<=1
%     p=1
% end
length=ceil(length);
oritation=ceil(oritation);
model=false(7,length);
model(3:5,1)=true;
model(3:5,length)=true;
model(2:6,2)=true;
model(2:6,length-1)=true;
model(:,3:length-2)=true;
model=imrotate(model,oritation);
% model=bwmorph(model,'remove');
end

