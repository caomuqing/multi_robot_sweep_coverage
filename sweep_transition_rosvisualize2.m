close all hidden
clear Frame3 Frame4
str3=["$\psi^{j+1}_1(k)$";"$\psi^{j+1}_2(k)$";"$\psi^{j+1}_3(k)$";"$\psi_4(k)$";"$\psi^q_5(k)$";"$\psi^q_6(k)$";"$\psi^q_7(k)$";"$\psi^q_8(k)$"];

uav1plot=plot(1,1,'Color','#B12C21', 'linewidth',3);
hold on
uav2plot=plot(1,1,'Color','#275791', 'linewidth',3);
uav3plot=plot(1,1,'Color','#E49600', 'linewidth',3);
xlabel('Time(s) Since Start of Operation');
ylabel('Operation cycle')
change_time(1)=1;
j=1;
residueee=4;
for i=1:size(uav1Msg,1)-1
    uav1psi(i)=uav1Msg{i}.Point.X;
    uav2psi(i)=uav2Msg{i}.Point.X;
    uav3psi(i)=uav3Msg{i}.Point.X;
    time(i)=uav1Msg{i}.Header.Stamp.Sec+ uav1Msg{i}.Header.Stamp.Nsec/1e9 - timeStart;
    
    if uav1Msg{i+1}.Point.Z~=uav1Msg{i}.Point.Z
        change_time(j+1)=i+1+residueee;
        j=j+1;
    end
end
change_time(j+1)=size(uav1Msg,1);
sannotate=annotation('textbox', [0.3, 0.8, 0.1, 0.1],'EdgeColor', 'white','FontSize',15);
for i=1:size(uav1Msg,1)
    for ii=1:j
        if i>=change_time(ii)&&i<change_time(ii+1)-residueee
            legend({str3(1),str3(2),str3(3)},'Interpreter','latex','FontSize',14, 'Orientation','vertical',...
                'Location','southeast')
            sannotate.String= ['partitioning stripe  ',int2str(ii+1)];
            uav1plot.XData=time(change_time(ii):i);
            uav1plot.YData=uav1psi(change_time(ii):i);
            uav2plot.XData=time(change_time(ii):i);
            uav2plot.YData=uav2psi(change_time(ii):i);
            uav3plot.XData=time(change_time(ii):i);
            uav3plot.YData=uav3psi(change_time(ii):i);
            drawnow
            break;
        end
    end
    Frame3(i) = getframe(gcf);
    pause(0.1);
end

for wer=1:residueee+5
    Frame_tmp(wer)=Frame3(change_time(2)-residueee);
end
Frame4=[Frame3(1:change_time(2)-residueee) Frame_tmp];
for wer=1:residueee+10
    Frame_tmp(wer)=Frame3(change_time(3)-residueee);
end
Frame4=[Frame4 Frame3(change_time(2):change_time(3)-residueee) Frame_tmp];
clear Frame_tmp;
for wer=1:residueee+6
    Frame_tmp(wer)=Frame3(change_time(4)-residueee);
end
Frame4=[Frame4 Frame3(change_time(3):change_time(4)-residueee) Frame_tmp];
clear Frame_tmp;
for wer=1:residueee+8
    Frame_tmp(wer)=Frame3(change_time(5)-residueee);
end
Frame4=[Frame4 Frame3(change_time(4):change_time(5)-residueee) Frame_tmp];
Frame4=[Frame4 Frame3(change_time(5):change_time(6)-residueee)];

  % create the video writer with 1 fps
  writerObj = VideoWriter('rosvisualize2.avi');
  writerObj.FrameRate = 4;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(Frame4)
    % convert the image to a frame
    frame = Frame4(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
