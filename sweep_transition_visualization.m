close all hidden
clear Frame2
colours1=["#E8E8E8";"#B8B8B8";"#8C8C8C";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
statescolor=["#ED7D31";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];
statescolor2=["#7030A0";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];

states_name=["transiting";"sweeping";"waiting";"moving up";"finished"];
robots=["robot 1:";"robot 2:";"robot 3:";"robot 4:";"robot 5:";"robot 6:";"robot 7:";"robot 8:";"robot 9:";"robot 10:"];

figure
skip_a_robot=false;

regions=gobjects(p.num_vert,p.num_horiz);
polygons=gobjects(p.num_vert,2);
boundary_polygons=gobjects(p.num_vert,2);

for i=1:p.num_vert
    for j=1:p.num_horiz
        for k=1:size(finite_omega_values,2)
            if p.building(i,j)==finite_omega_values(k)
                if j==1
                    regions(i,j)=rectangle('EdgeColor','black','FaceColor',colours1(k));
                    regions(i,j).Position=[x_0(p.num_vert-i+1,1), (p.num_vert-i)*p.building_block_height,...
                        p.building_block_length-x_0(p.num_vert-i+1,1),p.building_block_height]; 
                    myvertices=[0 (p.num_vert-i+1)*p.building_block_height; 0 (p.num_vert-i)*p.building_block_height];
                    for ii=(p.num_vert-i)*p.building_block_height/0.1+1:(p.num_vert-i+1)*p.building_block_height/0.1+1
                        myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                    end
                    f=1:size(myvertices,1);
                    polygons(i,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','white');
                    break;                    
                elseif j==p.num_horiz
                    regions(i,j)=rectangle('EdgeColor','black','FaceColor',colours1(k));
                    regions(i,j).Position=[(j-1)*p.building_block_length, (p.num_vert-i)*p.building_block_height,...
                        x_n(p.num_vert-i+1,2)-(j-1)*p.building_block_length,p.building_block_height];
                    myvertices=[p.num_horiz*p.building_block_length (p.num_vert-i+1)*p.building_block_height;...
                        p.num_horiz*p.building_block_length (p.num_vert-i)*p.building_block_height];
                    for ii=(p.num_vert-i)*p.building_block_height/0.1+1:(p.num_vert-i+1)*p.building_block_height/0.1+1
                        myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                    end
                    f=1:size(myvertices,1);
                    polygons(i,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','white');
                    break;                    
                else
                    regions(i,j)=rectangle('EdgeColor','black','FaceColor',colours1(k));
                    regions(i,j).Position=[(j-1)*p.building_block_length, (p.num_vert-i)*p.building_block_height,...
                        p.building_block_length,p.building_block_height];
                    break;
                end
            end
        end
    end
    boundary_polygons(i,1)=patch('Faces',1:3,'Vertices',[0 0;0 0;0 0],'FaceColor','green');
    boundary_polygons(i,2)=patch('Faces',1:3,'Vertices',[0 0;0 0;0 0],'FaceColor','green');
end


drawnow
hold on
plot(boundary(:,2),boundary(:,1),'Color','black','LineWidth',2.0);
plot(boundary(:,3),boundary(:,1),'Color','black','LineWidth',2.0);

posnn = get(gca, 'Position');
posnn(1) = 0.045;
posnn(3) = 0.65;
set(gca, 'Position', posnn)
set(gcf, 'Position',  [10, 100, 800, 450])

partition_bar=gobjects(n,p.num_vert);
finished_part=gobjects(n,p.num_vert);
finished_boundary=gobjects(2,p.num_vert);
for i=1:n
    for j=1:p.num_vert
        partition_bar(i,j)=rectangle('EdgeColor','blue','FaceColor','blue');
        partition_bar(i,j).Position=[p.building_block_length*p.num_horiz/n*i-p.building_block_length/3000000,...
        0,p.building_block_length/15000000,p.building_block_height/100000];
        hold on
        finished_part(i,j)=rectangle('EdgeColor','b','FaceColor','g');
    end
end
for i=2:n
    partition_bar(i,1)=rectangle('EdgeColor','blue','FaceColor','blue');
    partition_bar(i,1).Position=[partition(i,1,1)-p.building_block_length/30,...
        0,p.building_block_length/15,p.building_block_height];
end
drawnow
hold on
for i=1:n
    annotation('textbox', [0.75, n*0.06-(i-1)*0.06, 0.1, 0.1], 'String', robots(i), 'EdgeColor', '#F0F0F0','FontSize',12);
    hold on
end
for i=1:n
    statesannotate(i)=annotation('textbox', [0.85, n*0.06-(i-1)*0.06, 0.1, 0.1],'EdgeColor', '#F0F0F0','FontSize',12);
end

ph =gobjects(n,1);
for i=1:n
    if (skip_a_robot&&i==3)
        continue
    end
    ph(i) = plot(pos_x(i,1),pos_y(i,1),'or','MarkerSize',15,'MarkerFaceColor','b','MarkerEdgeColor','b');  %plot initial point
    hold on
end
set(gca,'Xlim',[0 p.building_block_length*p.num_horiz], 'Ylim',[0 p.building_block_height*p.num_vert]);
%axis([0 p.building_block_length*p.num_horiz 0 p.building_block_height*p.num_vert])
j=1;
for k=1:size(pos_x(1,:),2)

    
    for num=2:p.num_vert+1
        if k<=sum(kmax(2:num))+(num-1)*(k_vert-1)
            j=num;
            break;
%         else
%             j=0;
        end
    end
    for i=1:n
        if rem(j,2)==1
            statesannotate(i).String=states_name(states_all(i,k));
            statesannotate(i).Color=statescolor(states_all(i,k));   
        else
            statesannotate(i).String=states_name(states_all(i,k));
            statesannotate(i).Color=statescolor2(states_all(i,k)); 
        end
    end
    if j==2
        kk=k;
    else 
        kk=k-sum(kmax(2:j-1))-(j-2)*(k_vert-1);
    end
    if j~=p.num_vert+1 && kk<=kmax(j)
        for i=2:n
            %partition_bar(i,j)=rectangle('EdgeColor','g','FaceColor','g');
            partition_bar(i,j).Position=[partition(i,j,kk)-p.building_block_length/30,...
                (j-1)*p.building_block_height,p.building_block_length/15,p.building_block_height];
        end        
    end
    for i=1:n
        if (skip_a_robot&&i==3)
            continue
        end
        if rem(j-1,2)==1&& kk<=kmax(j)
            if kk>finish_transition_time(i,j-1)
                if i==1
                    if pos(i,j-1,kk)<x_0(j-1,2)
                        bound_low=boundary((j-2)*p.building_block_height/0.1+1,2);
                        bound_up=boundary((j-1)*p.building_block_height/0.1+1,2);
                        if pos(i,j-1,kk)>=bound_low&&pos(i,j-1,kk)<bound_up
                            myvertices=[pos(i,j-1,kk) (j-2)*p.building_block_height];
                            for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                                if boundary(ii,2)<=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;                            
                            boundary_polygons(j-1,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');          
                        elseif pos(i,j-1,kk)<bound_low&&pos(i,j-1,kk)>=bound_up
                            myvertices=[pos(i,j-1,kk) (j-1)*p.building_block_height];
                            for ii=(j-1)*p.building_block_height/0.1+1:(j-2)*p.building_block_height/0.1+1
                                if boundary(ii,2)<=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(j-1,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');          
                        else
                            myvertices=[];
                            for ii=(j-1)*p.building_block_height/0.1+1:(j-2)*p.building_block_height/0.1+1
                                myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(j-1,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');                                 
                        end
                    else
                        if abs(pos(i,j-1,kk)-x_0(j-1,2))>=p.building_block_length/30
                            myvertices=[pos(i,j-1,kk)-p.building_block_length/30 (j-1)*p.building_block_height;...
                                pos(i,j-1,kk)-p.building_block_length/30 (j-2)*p.building_block_height];
                        else
                            myvertices=[pos(i,j-1,kk) (j-1)*p.building_block_height;...
                                pos(i,j-1,kk) (j-2)*p.building_block_height];
                        end
                        for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                            myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                        end
                        f=1:size(myvertices,1);
%                         boundary_polygons(i,1).Faces=f;                            
%                         boundary_polygons(i,1).Vertices=myvertices;      
                        boundary_polygons(j-1,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');                          
                    end
                elseif i==n
                    if pos(i,j-1,kk)>x_n(j-1,1)
                        bound_low=boundary((j-2)*p.building_block_height/0.1+1,3);
                        bound_up=boundary((j-1)*p.building_block_height/0.1+1,3);
                        if pos(i,j-1,kk)>=bound_low&&pos(i,j-1,kk)<bound_up
                            myvertices=[...
                                pos(i,j-1,kk) (j-1)*p.building_block_height;...
                                partition(i,j-1,kmax(j-1)) (j-1)*p.building_block_height;...
                                partition(i,j-1,kmax(j-1)) (j-2)*p.building_block_height];
                            for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                                if boundary(ii,3)<=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;                            
                            boundary_polygons(j-1,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');          
                        elseif pos(i,j-1,kk)<bound_low&&pos(i,j-1,kk)>=bound_up
                            myvertices=[...
                                pos(i,j-1,kk) (j-2)*p.building_block_height;...
                                partition(i,j-1,kmax(j-1)) (j-2)*p.building_block_height;...
                                partition(i,j-1,kmax(j-1)) (j-1)*p.building_block_height];
                            for ii=(j-1)*p.building_block_height/0.1+1:-1:(j-2)*p.building_block_height/0.1+1
                                if boundary(ii,3)<=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(j-1,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');          
                        else
                            myvertices=[partition(i,j-1,kmax(j-1)) (j-1)*p.building_block_height;...
                                partition(i,j-1,kmax(j-1)) (j-2)*p.building_block_height];
                            for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                                myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(i,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');                                 
                        end
                    else
                        myvertices=[partition(i,j-1,kmax(j-1)) (j-1)*p.building_block_height;...
                                partition(i,j-1,kmax(j-1)) (j-2)*p.building_block_height;...
                                pos(i,j-1,kk) (j-2)*p.building_block_height;...
                                pos(i,j-1,kk) (j-1)*p.building_block_height];
                        f=1:size(myvertices,1);
%                         boundary_polygons(i,1).Faces=f;                            
%                         boundary_polygons(i,1).Vertices=myvertices;      
                        boundary_polygons(i,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');                          
                    end
                else
                    finished_part(i,j-1).Position=[partition(i,j-1,kmax(j-1)),(j-2)*p.building_block_height...
                        abs(pos(i,j-1,kk)-partition(i,j-1,kmax(j-1))),p.building_block_height];
                end
            end
        elseif kk<=kmax(j)
            if kk>finish_transition_time(i,j-1)
                if i==1
                    if pos(i,j-1,kk)<x_0(j-1,2)
                        bound_low=boundary((j-2)*p.building_block_height/0.1+1,2);
                        bound_up=boundary((j-1)*p.building_block_height/0.1+1,2);
                        if pos(i,j-1,kk)>=bound_low&&pos(i,j-1,kk)<bound_up
                            myvertices=[...
                                pos(i,j-1,kk) (j-2)*p.building_block_height;...
                                partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-2)*p.building_block_height;...
                                partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-1)*p.building_block_height];
                            for ii=-(j-1)*p.building_block_height/0.1-1:-(j-2)*p.building_block_height/0.1-1
                                if boundary(-ii,2)>=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(-ii,2) boundary(-ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;                            
                            boundary_polygons(j-1,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green','EdgeColor','none');          
                        elseif pos(i,j-1,kk)<bound_low&&pos(i,j-1,kk)>=bound_up
                            myvertices=[...
                                pos(i,j-1,kk) (j-1)*p.building_block_height;...
                                partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-1)*p.building_block_height;...
                                partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-2)*p.building_block_height];
                            for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                                if boundary(ii,2)>=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(j-1,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green','EdgeColor','none');          
                        else
                            myvertices=[partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-1)*p.building_block_height;...
                                partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-2)*p.building_block_height];
                            for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                                myvertices=[myvertices;boundary(ii,2) boundary(ii,1)];
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(i,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green','EdgeColor','none');                                 
                        end
                    else
                        if abs(partition(i+1,j-1,kmax(j-1))-pos(i,j-1,kk))<=p.building_block_length/30
                            myvertices=[partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-1)*p.building_block_height;...
                                    partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-2)*p.building_block_height;...
                                    pos(i,j-1,kk)-p.building_block_length/30 (j-2)*p.building_block_height;...
                                    pos(i,j-1,kk)-p.building_block_length/30 (j-1)*p.building_block_height];
                        else
                            myvertices=[partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-1)*p.building_block_height;...
                                    partition(i+1,j-1,kmax(j-1))-p.building_block_length/30 (j-2)*p.building_block_height;...
                                    pos(i,j-1,kk) (j-2)*p.building_block_height;...
                                    pos(i,j-1,kk) (j-1)*p.building_block_height];
                        end
                        f=1:size(myvertices,1);
%                         boundary_polygons(i,1).Faces=f;                            
%                         boundary_polygons(i,1).Vertices=myvertices;      
                        boundary_polygons(i,1)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green','EdgeColor','none');                          
                    end
                elseif i==n
                    if pos(i,j-1,kk)>x_n(j-1,1)
                        bound_low=boundary((j-2)*p.building_block_height/0.1+1,3);
                        bound_up=boundary((j-1)*p.building_block_height/0.1+1,3);
                        if pos(i,j-1,kk)>=bound_low&&pos(i,j-1,kk)<bound_up
                            myvertices=[pos(i,j-1,kk) (j-1)*p.building_block_height];
                            for ii=(j-1)*p.building_block_height/0.1+1:-1:(j-2)*p.building_block_height/0.1+1
                                if boundary(ii,3)>=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;                            
                            boundary_polygons(j-1,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');          
                        elseif pos(i,j-1,kk)<bound_low&&pos(i,j-1,kk)>=bound_up
                            myvertices=[pos(i,j-1,kk) (j-2)*p.building_block_height];
                            for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                                if boundary(ii,3)>=pos(i,j-1,kk)
                                    myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                                else
                                    break;
                                end
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(j-1,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');          
                        else
                            myvertices=[];
                            for ii=(j-1)*p.building_block_height/0.1+1:(j-2)*p.building_block_height/0.1+1
                                myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                            end
                            f=1:size(myvertices,1);
%                             boundary_polygons(i,1).Faces=f;                            
%                             boundary_polygons(i,1).Vertices=myvertices;      
                            boundary_polygons(j-1,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');                                 
                        end
                    else
                        myvertices=[pos(i,j-1,kk) (j-1)*p.building_block_height;...
                            pos(i,j-1,kk) (j-2)*p.building_block_height];
                        for ii=(j-2)*p.building_block_height/0.1+1:(j-1)*p.building_block_height/0.1+1
                            myvertices=[myvertices;boundary(ii,3) boundary(ii,1)];
                        end
                        f=1:size(myvertices,1);
%                         boundary_polygons(i,1).Faces=f;                            
%                         boundary_polygons(i,1).Vertices=myvertices;      
                        boundary_polygons(j-1,2)=patch('Faces',f,'Vertices',myvertices,'FaceColor','green');                          
                    end
                else
                    finished_part(i,j-1).Position=[pos(i,j-1,kk),(j-2)*p.building_block_height...
                         abs(pos(i,j-1,kk)-partition(i+1,j-1,kmax(j-1))),p.building_block_height];
                end
            end
        end
    end
    for i=1:n
        if (skip_a_robot&&i==3)
            continue
        end        
        ph(i).XData = pos_x(i,k);         %change x coordinate of the point
        ph(i).YData = pos_y(i,k);         %change y coordinate of the point
%                     uistack(ph(i),'top') 

    end    
    drawnow
%     uistack(ph(1),'top') 
%     uistack(ph(5),'top') 
    Frame2(k) = getframe(gcf);
%     pause(0.01)  %control speed, if desired    
end


  % create the video writer with 1 fps
  writerObj = VideoWriter('test345.avi');
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
