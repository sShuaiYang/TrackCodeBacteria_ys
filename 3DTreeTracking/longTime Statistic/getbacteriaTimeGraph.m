function getbacteriaTimeGraph(bacteriaTime,dirFile)
dirSave=strcat(dirFile,'\longTimeResult\length');
mkdir(dirSave)
dirSave1=strcat(dirSave,'\jpgPicture');
dirSave2=strcat(dirSave,'\figPicutre');
dirSaveW=strcat(dirFile,'\longTimeResult\width');
mkdir(dirSaveW)
dirSave3=strcat(dirSaveW,'\jpgPicture');
dirSave4=strcat(dirSaveW,'\figPicutre');
mkdir(dirSave1)
mkdir(dirSave2)
mkdir(dirSave3)
mkdir(dirSave4)
for t=1:size(bacteriaTime,2)
    for i=1:numel(bacteriaTime)
        if bacteriaTime(1,i).branchIndex==t
            % figure;
            a=bacteriaTime(1,i).LengthSeries;
            
            length=a(:,2);
            [~,sorh,keepapp] = ddencmp('den','wv',length);
            thr=150;
            NewData=wdencmp('gbl',length,'db10',5,thr,sorh,keepapp);
            
            [regionMax,regionMin]=findRegionPeak(NewData);
            figure
            plot(a(:,1),a(:,2),'g');
            frameInfo=a(:,1);
            hold on
%             plot(frameInfo,NewData,'r');
            plot(frameInfo(regionMax),length(regionMax),'bo')
            plot(frameInfo(regionMin),length(regionMin),'bx')
            plot(frameInfo(logical(a(:,3))),NewData(logical(a(:,3))),'rd')
            title(strcat('branch list is',num2str(bacteriaTime(1,i).branchIndex)));
            h=gcf;
            saveas(h,strcat(dirSave1,'\',num2str(i)),'jpg');
            saveas(h,strcat(dirSave2,'\',num2str(i)),'fig');
            
            
            figure
            width=a(:,4);
            plot(a(:,1),width,'g');
            hold on
            plot(frameInfo(regionMax),width(regionMax),'bo')
            plot(frameInfo(regionMin),width(regionMin),'bx')
            plot(frameInfo(logical(a(:,3))),width(logical(a(:,3))),'rd')
            h=gcf;
            saveas(h,strcat(dirSave3,'\',num2str(i)),'jpg');
            saveas(h,strcat(dirSave4,'\',num2str(i)),'fig');
            title(strcat('branch list is',num2str(bacteriaTime(1,i).branchIndex)));
            close all
            break
        end
    end
end
end
function [regionMaxOri,regionMinOri]=findRegionPeak(NewData)
regionMaxOri=imregionalmax(NewData);
regionMinOri=imregionalmin(NewData);
regionMinOri(1)=0;
regionMinOri(end)=0;
for i=1:numel(regionMaxOri)
    if regionMaxOri(i)==1
        whetherhavemin=0;
        for j=1:numel(regionMinOri)-i
            if regionMinOri(i+j)==1
                whetherhavemin=1;
                if NewData(i)-NewData(i+j)<=4
                    regionMaxOri(i)=0;
                    regionMinOri(i+j)=0;
                else
                    for u=1:min(100,numel(regionMaxOri)-i-j)
                        if regionMinOri(i+j+u)==1
                            if NewData(i+j+u)-NewData(i+j)<-4
                                regionMinOri(i+j)=0;
                                for p=1:u
                                    if regionMaxOri(i+j+p)==1
                                        if NewData(i+j+p)<NewData(i)
                                            regionMaxOri(i+j+p)=0;
                                            if NewData(i)-NewData(i+j+u)<=20
                                                regionMaxOri(i)=0;
                                                regionMinOri(i+j+u)=0;
                                            end
                                        else
                                            regionMaxOri(i)=0;
                                        end
                                        break
                                    end
                                end
                            else
                                if NewData(i)-NewData(i+j)<=20
                                    regionMaxOri(i)=0;
                                    regionMinOri(i+j)=0;
                                end
                            end
                            break
                        end
                        if u==min(100,numel(regionMaxOri)-i-j)
                            if NewData(i)-NewData(i+j)<=20
                                regionMaxOri(i)=0;
                                regionMinOri(i+j)=0;
                            end
                        end
                    end
                end
                break
            end
        end
        if whetherhavemin==0
            regionMaxOri(i)=0;
        end
    end
end
regionMaxOri(end)=0;
[regionMaxOri,regionMinOri]=findDivitionPart(regionMaxOri,regionMinOri,NewData);
end
function [regionMaxOri,regionMinOri]=findDivitionPart(regionMaxOri,regionMinOri,NewData)
maxRegion=find(regionMaxOri==1);
minRegion=find(regionMinOri==1);
for i=1:numel(maxRegion)
    if NewData(minRegion(i))<=20
        regionMaxOri(maxRegion(i))=0;
        regionMinOri(minRegion(i))=0;
    end
end
end