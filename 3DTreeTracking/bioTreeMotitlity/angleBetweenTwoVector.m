function angle=angleBetweenTwoVector(vector1,vector2) %return the angle ([-180,180]) betwwen the vector1 and vector2

angle1=atan(vector1(:,2)./abs(vector1(:,1))).*(180/pi);
angle1(vector1(:,1)<0&vector1(:,2)>=0)= 180-angle1(vector1(:,1)<0&vector1(:,2)>=0);
angle1(vector1(:,1)<0&vector1(:,2)<0)= -180-angle1(vector1(:,1)<0&vector1(:,2)<0);

angle2=atan(vector2(:,2)./abs(vector2(:,1))).*(180/pi);
angle2(vector2(:,1)<0&vector2(:,2)>=0)= 180-angle2(vector2(:,1)<0&vector2(:,2)>=0);
angle2(vector2(:,1)<0&vector2(:,2)<0)= -180-angle2(vector2(:,1)<0&vector2(:,2)<0);

angle=angle2- angle1;
angle(angle>180)=angle(angle>180)-360;
angle(angle<-180)=angle(angle<-180)+360;
end