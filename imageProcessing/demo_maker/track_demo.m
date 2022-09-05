function track_demo(tracks,startframe,endframe)
nbeads=length(tracks);
allframes=startframe:endframe;
for iframe=allframes, % loop over all frame
    for i=1:nbeads; % loop over all beads
        %         if(tracks(i).startFrame<=(iframe+1));
        if(tracks(i).startFrame<=(iframe+1) && tracks(i).endFrame>=(iframe-1));
            trace_i=tracks(i).frameDetail;
            d_frame=iframe-tracks(i).startFrame+1;
            if (d_frame==0)
                trace_l=1;
            else
                trace_l=min([d_frame,length(trace_i)]);
            end
            C_xy=[trace_i(1:trace_l).Centroid];
            x_trace=[C_xy([1:2:length(C_xy)]),0];y_trace=[C_xy([2:2:length(C_xy)]),0];
            if ~isempty(tracks(i).DivideFrame)
                divdeindex=find(tracks(i).DivideFrame==iframe, 1);
                if ~isempty(divdeindex)
                    dvideFrame=[tracks(i).DivideFrame(divdeindex),0];
                else
                    dvideFrame=[0,0];
                end
            else
                dvideFrame=[0,0];
            end
            if(tracks(i).endFrame>=iframe)
                num=[iframe,i,trace_l,0,tracks(i).startFrame,tracks(i).endFrame,tracks(i).possibleParent,tracks(i).possibleChildren];
            else
                num=[iframe,i,trace_l,1,tracks(i).startFrame,tracks(i).endFrame,tracks(i).possibleParent,tracks(i).possibleChildren];
            end
            MIJ.setColumn('num',num);
            MIJ.setColumn('dvideFrame',dvideFrame);
            MIJ.setColumn('xtrace',x_trace); MIJ.setColumn('ytrace',y_trace);
            MIJ.run('Trace Maker');
        end
    end;
    %fprintf(1,'.');
end
end
