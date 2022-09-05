while(1)
    HostI=SearchAndProcess(HostI);
    if HostI.ProcessedImageNumber>=HostI.TotalImageNumber
        break;
    end
    pause(3);
end