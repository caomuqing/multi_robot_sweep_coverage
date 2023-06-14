close all
colours1=["#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
fig1=figure(q)
if ishandle(fig1)
clf(fig1)
end
fig1

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

axis([0  la    0  lb+0.5])
fullregion=rectangle('EdgeColor','k');
fullregion.Position=[0 0    la  lb];
drawnow
%h = rectangle('EdgeColor','k','FaceColor','b');

partition_bar=gobjects(n-1,q);
% robot_pos=gobjects(n,p.num_vert);
% finished_boundary=gobjects(2,p.num_vert);
sweeping_progress=gobjects(n,q);
partition_progress=gobjects(n,q);

for i=1:n

    for j=1:q
            sweeping_progress(i,j)=rectangle('EdgeColor','b','FaceColor','g'); 
        sweeping_progress(i,j).Position=[0,...
        0,1/15000000,1/100000];
        hold on
            partition_progress(i,j)=rectangle('EdgeColor','b','FaceColor',colours1(i));        
        partition_progress(i,j).Position=[0,...
        0,1/15000000,1/100000];
        hold on
    end
end
  for i=1:n-1
    for j=1:q       
        partition_bar(i,j)=rectangle('EdgeColor','#009004','FaceColor','#009004');
        partition_bar(i,j).Position=[0,...
        0,1/15000000,1/100000];
        hold on
%         robot_pos(i,j)=rectangle('EdgeColor','b','FaceColor','g');
    end
end      
kk=1;
for j=1:q
    for i=1:n
        final_workload(i)=rectangle('EdgeColor','k');
        final_workload(i).Position=[mark(i,j,k_max1(j)), (j-1)*d, mark(i+1,j,k_max1(j))-mark(i,j,k_max1(j)),d];
    end
    drawnow
    for k=1:k_max1(j+1)

        for li=1:n
           %sweeping_progress(li)=rectangle('EdgeColor','b','FaceColor','g'); 
           if (rem(j,2)==1)
               sweeping_progress(li,j).Position=[mark(li,j,k_max1(j)), (j-1)*d, pos1(li,j,k)-mark(li,j,k_max1(j)),d];
           else
               sweeping_progress(li,j).Position=[pos1(li,j,k), (j-1)*d, -pos1(li,j,k)+mark(li+1,j,k_max1(j)),d];               
           end
           
           if j~=q
               %partition_progress(li)=rectangle('EdgeColor','b','FaceColor',colours1(li)); 
               partition_progress(li,j).Position=[mark(li,j+1,k), (j)*d, mark(li+1,j+1,k)-mark(li,j+1,k),d];
           end
           if j~=q
           if li~=n
            partition_bar(li,j).Position=[mark(li+1,j+1,k)-1/30,(j)*d,1/15,d];    
           end
           end
           
        end
        Frame1(kk) = getframe(gcf);
        kk=kk+1;
        drawnow
        %pause(0.001)
    %drawnow
    end
end

  % create the video writer with 1 fps
  writerObj = VideoWriter('algo1.avi');
  writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(Frame1)
    % convert the image to a frame
    frame = Frame1(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);