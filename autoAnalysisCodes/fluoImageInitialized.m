function fluoImageInitialized()
dirFile=uigetdir();
load([dirFile,'\mip.mat'])
fluo_calibration=mip.Calib.fluo_protein_calibration1_6;
fluoChannel=mip.Calib.fluo_protein_info;
fieldNameList=dir(dirFile);
sita=144.9;
for iField=1:numel(fieldNameList)-2
    fieldNum=0;
    dirFieldFile=[dirFile,'\',fieldNameList(iField+2).name];
    if numel(fieldNameList(iField+2).name)>=9 && strcmp(fieldNameList(iField+2).name(end-8:end-4),'feild')
        fieldNum=fieldNum+1;
        if fieldNum==1
            fluoHere=[];
            n=0;
            for iFluo=1:numel(fluoChannel)
                try load([dirFieldFile,'\',fluoChannel{iFluo},'\frameInfo.mat']);
                    n=n+1;
                    fluoChannelChoose(n)=iFluo;
                    fluoHere=[fluoHere,fluoChannel(iFluo)];
                    mkdir([dirFieldFile,'\',fluoChannel{iFluo},'new']);
                    switch iFluo
                        case 1
                            backGround(:,:,n)=mip.Calib.LumencorBlueField;
                        case 2
                            backGround(:,:,n)=mip.Calib.LumencorTealField;
                        case 3
                            backGround(:,:,n)=mip.Calib.LumencorVioletField;
                        case 4
                            backGround(:,:,n)=mip.Calib.LumencorCyanField;
                        case 5
                            backGround(:,:,n)=mip.Calib.LumencorCyanField;
                        case 6
                            backGround(:,:,n)=mip.Calib.LumencorGreenField;
                        case 7
                            backGround(:,:,n)=mip.Calib.LumencorRedField;
                    end
                end
            end
        end
        for i=1:6
            for iFluo=1:numel(fluoHere)
                load([dirFieldFile,'\',fluoHere{iFluo},'\frameInfo.mat']);
                image=import_tiff_stack([dirFieldFile,'\',fluoHere{iFluo},'\image',fluoHere{iFluo},num2str(i,'%05.f'),'.tif']);
                image=double(image)-100;
                image=image./backGround(:,:,iFluo);
                image=image/frameInfo(i,14)*frameInfo(i,13)/frameInfo(i,9);
                image=image(:);
                newMatrix(iFluo,:)=double(image);
            end
            image=[];caliMatrix=fluo_calibration(fluoChannelChoose,fluoChannelChoose);
            newMatrix=inv(caliMatrix)'*newMatrix;
            for iFluo=1:numel(fluoHere)
                image(:,:,iFluo)=reshape(newMatrix(iFluo,:),2048,2048)*sita/1000;
            end
        end
    end
end
end