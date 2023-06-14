fig=figure(q+1)
if ishandle(fig)
clf(fig)
end
fig

colormap(flipud(gray(256)))
x = 0:0.1:la;
y = 0:0.1:lb;
[X,Y] = meshgrid(x,y);
Z = (rou_up-rou_low)/2*sin(X+Y)+(rou_up+rou_low)/2;
contourf(X,Y,Z,100,'LineColor','none')
c = colorbar;
c.Label.String = 'workload density $\rho(x,y)$';
c.Label.Interpreter='latex';
c.Label.FontSize=13;

posnn = get(gca, 'Position');
posnn(1) = 0.045;
posnn(3) = 0.50;
set(gca, 'Position', posnn)
set(gcf, 'Position',  [10, 100, 800, 450])


colours=['r', 'c', 'm', 'b'];
colours1=["#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
axis([0  la    0  lb+0.5])
fullregion2=rectangle('EdgeColor','k');
fullregion2.Position=[0 0    la  lb];
drawnow
%h = rectangle('EdgeColor','k','FaceColor','b');
kk=1;


partition_progress2=gobjects(n,q);
extra_workload=gobjects(n,q);
sweeping_progress2=gobjects(n,q);

partition_bar2=gobjects(n,q);
robot_pos=gobjects(n,q);

for i=1:n

    for j=1:q

        hold on
            partition_progress2(i,j)=rectangle('EdgeColor',colours1(i),'FaceColor',colours1(i));        
        partition_progress2(i,j).Position=[0,...
        0,1/15000000,1/100000];
    if rem(j,2)==1
        extra_workload(i,j)=rectangle('EdgeColor',colours1(i+1),'FaceColor',colours1(i+1));
    else
        if i==1
             extra_workload(i,j)=rectangle('EdgeColor',colours1(i),'FaceColor',colours1(i));
        else
        extra_workload(i,j)=rectangle('EdgeColor',colours1(i-1),'FaceColor',colours1(i-1));
        end
             
    end
    
                extra_workload(i,j).Position=[0,...
        0,1/15000000,1/100000];
    
                 sweeping_progress2(i,j)=rectangle('EdgeColor','b','FaceColor','g'); 
        sweeping_progress2(i,j).Position=[0,...
        0,1/15000000,1/100000];  
    end
end
  for i=1:n
    for j=1:q       
        partition_bar2(i,j)=rectangle('EdgeColor','#009004','FaceColor','#009004');
        partition_bar2(i,j).Position=[0,...
        0,1/15000000,1/100000];
        hold on
%         robot_pos(i,j)=rectangle('EdgeColor','b','FaceColor','g');
        robot_pos(i,j)=rectangle('EdgeColor','blue','FaceColor','blue');
        robot_pos(i,j).Position=[0,...
        0,1/15000000,1/100000];
    end
end  

for j=1:q
    for i=1:n
        final_workload2(i)=rectangle('EdgeColor','k');
        final_workload2(i).Position=[mark2(i,j,k_max2(j)), (j-1)*d, mark2(i+1,j,k_max2(j))-mark2(i,j,k_max2(j)),d];
    end
    if j~=q
       fullstripe=rectangle('EdgeColor','k');
       fullstripe.Position=[0 j*d la d];
    end
    drawnow

    for k=1:k_max2(j+1)
        if rem(j,2)==1
            for i=1:n
%                 sweeping_progress2(i)=rectangle('EdgeColor','b','FaceColor','g');
               if j~=q
%                    partition_progress2(i,j)=rectangle('EdgeColor','b','FaceColor',colours1(i)); 
                   if mark2(i,j+1,k+1)>mark2(i,j,k_max2(j))
                        partition_progress2(i,j).Position=[mark2(i,j+1,k+1), (j)*d, ...
                                                    mark2(i+1,j+1,k+1)-mark2(i,j+1,k+1),d];
                   else
                        partition_progress2(i,j).Position=[mark2(i,j,k_max2(j)), (j)*d, ...
                                                    mark2(i+1,j+1,k+1)-mark2(i,j,k_max2(j)),d];
                   end
                   if pos2(i,j,k)> mark2(i+1,j+1,k+1)
%                        extra_workload(i,j)=rectangle('EdgeColor','b','FaceColor',colours1(i+1));    
                       extra_workload(i,j).Position=[mark2(i+1,j+1,k+1), (j)*d, ...
                                    pos2(i,j,k)-mark2(i+1,j+1,k+1),d];
                   end
                    partition_bar2(i,j).Position=[mark2(i+1,j+1,k)-1/60,(j)*d,1/30,d];    
                    robot_pos(i,j).Position=[pos2(i,j,k)-1/60,(j)*d,1/30,d];  
               end
                sweeping_progress2(i,j).Position=[mark2(i,j,k_max2(j)), (j-1)*d, pos2(i,j,k)-mark2(i,j,k_max2(j)),d];
               
            end
        else
            for i=1:n
%                 sweeping_progress2(i)=rectangle('EdgeColor','b','FaceColor','g');
               if j~=q
%                    partition_progress2(i)=rectangle('EdgeColor','b','FaceColor',colours1(i)); 
                   if mark2(i,j+1,k+1)<mark2(i,j,k_max2(j))
                        partition_progress2(i,j).Position=[mark2(i,j+1,k+1), (j)*d, ...
                                                    mark2(i+1,j+1,k+1)-mark2(i,j+1,k+1),d];
                   else
                        partition_progress2(i,j).Position=[mark2(i,j+1,k+1), (j)*d, ...
                                                    -mark2(i,j+1,k+1)+mark2(i+1,j,k_max2(j)),d];
                   end
                   if pos2(i,j,k)< mark2(i,j+1,k+1)&&i~=1
%                        extra_workload(i,j)=rectangle('EdgeColor','b','FaceColor',colours1(i-1));    
                       extra_workload(i,j).Position=[pos2(i,j,k), (j)*d, ...
                                    -pos2(i,j,k)+mark2(i,j+1,k+1),d];
                   end
                    partition_bar2(i,j).Position=[mark2(i,j+1,k)-1/60,(j)*d,1/30,d];    
                    robot_pos(i,j).Position=[pos2(i,j,k)-1/60,(j)*d,1/30,d];  
                   
               end
                sweeping_progress2(i,j).Position=[pos2(i,j,k), (j-1)*d, -pos2(i,j,k)+mark2(i+1,j,k_max2(j)),d];
               
            end            
        end
        Frame2(kk) = getframe(gcf);
        kk=kk+1;
        drawnow

        %pause(0.001)
    %drawnow
    end
    if j~=q
        for i=1:n
           extra_workload(i,j).Position=[0,0,1/15000000,1/100000];
            partition_progress2(i,j).Position=[mark2(i,j+1,k_max2(j+1)), (j)*d, ...
                                        -mark2(i,j+1,k_max2(j+1))+mark2(i+1,j+1,k_max2(j+1)),d];           
        end    
    end
    
    for i=1:n
        robot_pos(i,j).Position=[0,0,1/15000000,1/100000];
    end
end

  % create the video writer with 1 fps
  writerObj = VideoWriter('algo2.avi');
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