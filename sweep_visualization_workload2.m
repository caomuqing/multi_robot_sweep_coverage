close all hidden
clear Frame2 Frame3 Frame4 frame2_tmp
colours1=["#E8E8E8";"#B8B8B8";"#8C8C8C";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
statescolor=["#ED7D31";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];
statescolor2=["#7030A0";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];
str3=["$\mu^{j+1}_1(k)$";"$\mu^{j+1}_2(k)$";"$\mu^{j+1}_3(k)$";"$\mu^{j+1}_4(k)$";"$\mu^{j+1}_5(k)$";"$\mu^{j+1}_6(k)$";"$\psi^q_7(k)$";"$\psi^q_8(k)$"];
colours2=["#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];


figure

hold on

sannotate=annotation('textbox', [0.3, 0.8, 0.1, 0.1],'EdgeColor', 'white','FontSize',15);

uav1plot=plot(1,1,'linewidth',3,'Color',colours2(1));
hold on
uav2plot=plot(1,1,'linewidth',3,'Color',colours2(2));
uav3plot=plot(1,1, 'linewidth',3,'Color',colours2(3));
uav4plot=plot(1,1,'linewidth',3,'Color',colours2(4));
uav5plot=plot(1,1, 'linewidth',3,'Color',colours2(5))
xlabel('Sampling Time k','fontweight','bold','fontsize',12);
ylabel('Workload Amount','fontweight','bold','fontsize',12)
j=1;
j_prev=2;
Frame3=getframe(gcf);
for j=2:q
    y_lim=max(m2(:,j,k_max2(j)))+0.1;
    x_lim=k_max2(j)+20;
    xlim([0 x_lim])
    ylim([0 y_lim])
    
    for kk=1:k_max2(j)
        legend({str3(1),str3(2),str3(3),str3(4),str3(5)},'Interpreter','latex','FontSize',14,...
            'Orientation','vertical','Location','southeast')
        sannotate.String= ['partitioning stripe  ',int2str(j)];
        uav1plot.XData=1:kk;
        uav1plot.YData=m2(1,j,1:kk);
        uav2plot.XData=1:kk;
        uav2plot.YData=m2(2,j,1:kk);
        uav3plot.XData=1:kk;
        uav3plot.YData=m2(3,j,1:kk);
        uav4plot.XData=1:kk;
        uav4plot.YData=m2(4,j,1:kk);
        uav5plot.XData=1:kk;
        uav5plot.YData=m2(5,j,1:kk);
        drawnow
        Frame3=[Frame3 getframe(gcf)];
    end

    
%     if j_prev~=j
%         for wh=1:k_vert-10
%             frame2_tmp(wh)=Frame2(size(Frame2,2));
%         end
%         Frame2=[Frame2 frame2_tmp getframe(gcf)];
%     else
%     end
    j_prev=j;

    %pause(0.01)  %control speed, if desired    
end


  % create the video writer with 1 fps
  writerObj = VideoWriter('algo2_workload.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(Frame3)
    % convert the image to a frame
    frame = Frame3(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
