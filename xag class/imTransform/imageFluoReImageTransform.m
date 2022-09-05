function imageFluoReImageTransform(B3TransformInfo,dirFile)

disp('Used to transform mScarletI Venus TDsmURFP sfGFP image to CyOFP.');
nameList=dir(dirFile);
for iField=1:length(nameList)-2
    if strcmp(nameList(iField+2).name(1:5),'feild')
        dirImageFile=strcat(dirFile,'\',nameList(iField+2).name);
        dirImageList=dir(dirImageFile);
        for jChannel=1:length(dirImageList)-2
            if strcmp(dirImageList(jChannel+2).name,'mScarletI')||strcmp(dirImageList(jChannel+2).name,'sfGFP')||strcmp(dirImageList(jChannel+2).name,'Venus')||strcmp(dirImageList(jChannel+2).name,'TDsmURFP')
                flouProtein=dirImageList(jChannel+2).name;
                dirImage=strcat(dirImageFile,'\',dirImageList(jChannel+2).name);
%                 dirNewImage=strcat(dirImageFile,'\',dirImageList(jChannel+2).name,'transform');
%                 mkdir(dirNewImage);
                imageList=dir(dirImage);
                for iImage=1:length(imageList)-2
                    if strcmp(imageList(iImage+2).name(1:5),'image')
                        imageFluo=import_tiff_stack(strcat(dirImage,'\',imageList(iImage+2).name));
                        switch flouProtein
                            case 'mScarletI'
                                [~,imageFluo]=imCameraTransform(imageFluo,imageFluo,B3TransformInfo.transformInfoZylaCyOFP2mScarletI,B3TransformInfo.bestPositionCyOFP2mScarletI);
                            case 'sfGFP'
                                [~,imageFluo]=imCameraTransform(imageFluo,imageFluo,B3TransformInfo.transformInfoZylaCyOFP2sfGFP,B3TransformInfo.bestPositionCyOFP2sfGFP);
                            case 'Venus'
                                [~,imageFluo]=imCameraTransform(imageFluo,imageFluo,B3TransformInfo.transformInfoZylaCyOFP2Venus,B3TransformInfo.bestPositionCyOFP2Venus);
                            case'TDsmURFP'
                                [~,imageFluo]=imCameraTransform(imageFluo,imageFluo,B3TransformInfo.transformInfoZylaCyOFP2TDsmURFP,B3TransformInfo.bestPositionCyOFP2TDsmURFP);
                        end
                        imwrite(imageFluo,strcat(dirImage,'\',imageList(iImage+2).name),'tiff');
                    end
                end
            end
        end
    end
    
    
end