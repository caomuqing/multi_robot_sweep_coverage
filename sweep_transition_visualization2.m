close all hidden
clear Frame2 Frame3 Frame4 frame2_tmp
colours1=["#E8E8E8";"#B8B8B8";"#8C8C8C";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
statescolor=["#ED7D31";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];
statescolor2=["#7030A0";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];
str3=["$\psi^{j+1}_1(k)$";"$\psi^{j+1}_2(k)$";"$\psi^{j+1}_3(k)$";"$\psi^{j+1}_4(k)$";"$\psi^{j+1}_5(k)$";"$\psi^{j+1}_6(k)$";"$\psi^q_7(k)$";"$\psi^q_8(k)$"];


figure

hold on

sannotate=annotation('textbox', [0.3, 0.8, 0.1, 0.1],'EdgeColor', 'white','FontSize',15);

uav1plot=plot(1,1,'linewidth',3);
hold on
uav2plot=plot(1,1,'linewidth',3);
uav3plot=plot(1,1, 'linewidth',3);
uav4plot=plot(1,1,'linewidth',3);
uav5plot=plot(1,1, 'linewidth',3)
xlabel('Sampling Time k','fontweight','bold','fontsize',12);
ylabel('Operation cycle','fontweight','bold','fontsize',12)
j=1;
j_prev=2;
Frame2=getframe(gcf);
for k=1:size(pos_x(1,:),2)

    
    for num=2:p.num_vert+1
        if k<=sum(kmax(2:num))+(num-1)*(k_vert-1)
            j=num;
            break;
%         else
%             j=0;
        end
    end

    if j==2
        kk=k;
    else 
        kk=k-sum(kmax(2:j-1))-(j-2)*(k_vert-1);
    end
    if j~=p.num_vert+1 && kk<=kmax(j)
        legend({str3(1),str3(2),str3(3),str3(4),str3(5)},'Interpreter','latex','FontSize',14,...
            'Orientation','vertical','Location','southeast')
        sannotate.String= ['partitioning stripe  ',int2str(j)];
        uav1plot.XData=1:kk;
        uav1plot.YData=m(1,j,1:kk);
        uav2plot.XData=1:kk;
        uav2plot.YData=m(2,j,1:kk);
        uav3plot.XData=1:kk;
        uav3plot.YData=m(3,j,1:kk);
        uav4plot.XData=1:kk;
        uav4plot.YData=m(4,j,1:kk);
        uav5plot.XData=1:kk;
        uav5plot.YData=m(5,j,1:kk);
    end

    drawnow
%     if j_prev~=j
%         for wh=1:k_vert-10
%             frame2_tmp(wh)=Frame2(size(Frame2,2));
%         end
%         Frame2=[Frame2 frame2_tmp getframe(gcf)];
%     else
        Frame2=[Frame2 getframe(gcf)];
%     end
    j_prev=j;

    %pause(0.01)  %control speed, if desired    
end


  % create the video writer with 1 fps
  writerObj = VideoWriter('test456.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(Frame2)
    % convert the image to a frame
    frame = Frame2(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
