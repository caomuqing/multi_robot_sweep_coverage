colours1=["#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
fig1=figure(q)
if ishandle(fig1)
clf(fig1)
end
fig1
axis([0  la    0  lb+0.5])
fullregion=rectangle('EdgeColor','k');
fullregion.Position=[0 0    la  lb];
drawnow
%h = rectangle('EdgeColor','k','FaceColor','b');
kk=1;
%     for i=1:n
%         final_workload(i)=rectangle('EdgeColor','k');
%         final_workload(i).Position=[mark(i,j,k_max1(j)), (j-1)*d, mark(i+1,j,k_max1(j))-mark(i,j,k_max1(j)),d];
%     end
    drawnow
    for i=1:n
        for j=1:q    
            partition_progress(i,j)=rectangle('EdgeColor','b','FaceColor',colours1(i));
        end
    end
    for i=1:n
        sweeping_progress(i)=rectangle('EdgeColor','b','FaceColor','g'); 

    end
    for k=1:k_max

        for li=1:n
           %sweeping_progress(li)=rectangle('EdgeColor','b','FaceColor','g'); 
           sweeping_progress(li).Position=[bars(li,pos_level(li,k),k), (pos_level(li,k)-1)*d, pos(li,k)-bars(li,pos_level(li,k),k),d];
           for j=1:q
               %partition_progress(li)=rectangle('EdgeColor','b','FaceColor',colours1(li)); 
               %partition_progress(li,j)=rectangle('EdgeColor','b','FaceColor',colours1(li));
               partition_progress(li,j).Position=[bars(li,j,k), (j-1)*d, bars(li+1,j,k)-bars(li,j,k),d];
           end
        end
        Frame1(kk) = getframe(gcf);
        kk=kk+1;
        drawnow
        pause(0.01)
    %drawnow
    end


%   % create the video writer with 1 fps
%   writerObj = VideoWriter('algo1.avi');
%   writerObj.FrameRate = 10;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(Frame1)
%     % convert the image to a frame
%     frame = Frame1(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);