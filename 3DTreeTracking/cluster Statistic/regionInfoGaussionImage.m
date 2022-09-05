function imageEnd=regionInfoGaussionImage(degreeRadAll)
imageEnd=zeros(size(degreeRadAll.data{1}.Image));
for i=1:size(degreeRadAll.data,2)
    imageEnd=imageEnd+degreeRadAll.data{i}.Image;
end
threShould=max(max(imageEnd))/20;
imageEnd(imageEnd>=threShould)=1;
imageEnd(imageEnd~=1)=0;
imageEnd=logical(imageEnd);
imageEnd=bwmorph(imageEnd,'remove');
end