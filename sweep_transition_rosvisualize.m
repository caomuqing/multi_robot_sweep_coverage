close all hidden
clear variables



finite_omega_values=[1, 2];
n= 3;
p.sampling_period=0.5; %has to be one
p.transition_speed= 1.4;
p.sweeping_speed= 1;
p.building_block_length = 0.3;
p.building_block_height= 0.2;
p.num_horiz= 12;
p.num_vert= 6;
kappa= 0.1;
p.building=zeros(p.num_vert, p.num_horiz);
p.building=[1,1,2,1,1,2,1,1,2,1,1,1;
            2,2,1,1,1,1,1,1,1,1,2,2;
            1,1,1,1,1,1,1,2,2,2,2,2;
            1,1,1,1,2,2,2,2,1,1,1,1;
            2,2,2,2,2,1,1,1,1,1,1,1;
            1,1,1,1,1,1,1,1,1,1,1,1];
    
colours1=["#E8E8E8";"#8C8C8C";"#B8B8B8";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
statescolor=["#ED7D31";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];
statescolor2=["#7030A0";"#24AC51";"#77AC30";"#000000";"#A2142F";"#7030A0";"#D95319"];

states_name=["transiting";"sweeping";"waiting";"moving up";"finished"];
robots=["robot 1:";"robot 2:";"robot 3:";"robot 4:";"robot 5:";"robot 6:";"robot 7:";"robot 8:";"robot 9:";"robot 10:"];

sampling_time = 0.5;
bag = rosbag('sweep4.bag');
bagselect1 = select(bag, 'Topic', '/uav1/workload');
bagselect2 = select(bag, 'Topic', '/uav2/workload');
bagselect3 = select(bag, 'Topic', '/uav3/workload');
uav1Msg = readMessages(bagselect1);
uav2Msg = readMessages(bagselect2);
uav3Msg = readMessages(bagselect3);
timeStart = uav1Msg{1}.Header.Stamp.Sec+ uav1Msg{1}.Header.Stamp.Nsec/1e9;

data = readtable('log1734.txt');
div=5;
q=fix(size(data(:,2))/div);
videoxmin=table2array(data(1,2));
videox_1_2=table2array(data(1,5));
videox_2_4=table2array(data(1,8));
videox_3_6=max(table2array(data(1:1000,8)));
videoxmax=max(table2array(data(:,8)));
videoymin=table2array(data(1,3));
videoymax=162;


for kk=1:q
    uav1X(kk)=positionxfrompixel(data{(kk-1)*div+1,2},videoxmin,videox_1_2,videox_2_4,videox_3_6);
    uav1Y(kk)=(data{(kk-1)*div+1,3}-videoymin)/(videoymax-videoymin)*p.num_vert*p.building_block_height;
    uav2X(kk)=positionxfrompixel(data{(kk-1)*div+1,5},videoxmin,videox_1_2,videox_2_4,videox_3_6);
    uav2Y(kk)=(data{(kk-1)*div+1,6}-videoymin)/(videoymax-videoymin)*p.num_vert*p.building_block_height;
    uav3X(kk)=positionxfrompixel(data{(kk-1)*div+1,8},videoxmin,videox_1_2,videox_2_4,videox_3_6);
    uav3Y(kk)=(data{(kk-1)*div+1,9}-videoymin)/(videoymax-videoymin)*p.num_vert*p.building_block_height;
end
figure
    xlim([-0.3 p.num_horiz*p.building_block_length+0.3]);
    ylim([-0.2 p.num_vert*p.building_block_height]+0.1);
    posnn = get(gca, 'Position');
posnn(1) = 0.10;
posnn(3) = 0.85;
set(gca, 'Position', posnn)
    set(gcf, 'Position',  [100, 100, 1200, 400])
%     title('UAV Trajectory from Recorded Data')
    
regions=gobjects(p.num_vert,p.num_horiz);
for i=1:p.num_vert
    for j=1:p.num_horiz
        for k=1:size(finite_omega_values,2)
            if p.building(i,j)==finite_omega_values(k)
                regions(i,j)=rectangle('EdgeColor','black','FaceColor',colours1(k));
                regions(i,j).Position=[(j-1)*p.building_block_length, (p.num_vert-i)*p.building_block_height,...
                    p.building_block_length,p.building_block_height];
                break;
            end
        end
    end
end

drawnow
hold on

partition_bar=gobjects(n,p.num_vert);
for i=1:n
    for j=1:p.num_vert
        partition_bar(i,j)=rectangle('EdgeColor','blue','FaceColor','blue');
        partition_bar(i,j).Position=[p.building_block_length*p.num_horiz/n*i-p.building_block_length/3000000,...
        0,p.building_block_length/15000000,p.building_block_height/100000];
        hold on
    end
end
for i=1:n
    partition_bar(i,1)=rectangle('EdgeColor','blue','FaceColor','blue');
    partition_bar(i,1).Position=[p.building_block_length*p.num_horiz/n*i-p.building_block_length/10,...
        0,p.building_block_length/5,p.building_block_height];
end
drawnow
% tnow1=tic;
% tnow2=tic;
kk=1;
ii=1;
% while true
%     if toc(tnow1)>=p.sampling_period/2 && kk<=size(uav1Msg,1)
%         tnow1=tic;
%         for j=1:p.num_vert-1
%             if uav1Msg{kk}.Point.Z==j
%                 partition_bar(1,j+1).Position=[uav1Msg{kk}.Point.Y-p.building_block_length/20,...
%                     (j)*p.building_block_height,p.building_block_length/10,p.building_block_height];       
%             end
%             if uav2Msg{kk}.Point.Z==j
%                 partition_bar(2,j+1).Position=[uav2Msg{kk}.Point.Y-p.building_block_length/20,...
%                     (j)*p.building_block_height,p.building_block_length/10,p.building_block_height];       
%             end
%             if uav3Msg{kk}.Point.Z==j
%                 partition_bar(3,j+1).Position=[uav3Msg{kk}.Point.Y-p.building_block_length/20,...
%                     (j)*p.building_block_height,p.building_block_length/10,p.building_block_height];       
%             end
%         end
%         kk=kk+1;
% 
%     end
%     if toc(tnow2)>=1/30*div
%         tnow2=tic;
%         plot(uav1X(1:ii),uav1Y(1:ii),'b-', 'linewidth',3);
%         plot(uav2X(1:ii),uav2Y(1:ii),'r-', 'linewidth',3);
%         plot(uav3X(1:ii),uav3Y(1:ii),'y-', 'linewidth',3);
%         ii=ii+1;
%         drawnow
%         if ii>size(uav1X)
%             break;
%         end
%     end
%     pause(0.0001)
% end

ratio=30/div/(2/p.sampling_period);
kk=1;
for j=1:p.num_vert-1
    if uav1Msg{kk}.Point.Z==j
        partition_bar(1,j+1).Position=[uav1Msg{kk}.Point.Y-p.building_block_length/10,...
            (j)*p.building_block_height,p.building_block_length/5,p.building_block_height];       
    end
    if uav2Msg{kk}.Point.Z==j
        partition_bar(2,j+1).Position=[uav2Msg{kk}.Point.Y-p.building_block_length/10,...
            (j)*p.building_block_height,p.building_block_length/5,p.building_block_height];       
    end
    if uav3Msg{kk}.Point.Z==j
        partition_bar(3,j+1).Position=[uav3Msg{kk}.Point.Y-p.building_block_length/10,...
            (j)*p.building_block_height,p.building_block_length/5,p.building_block_height];       
    end
end
base_value=4.2*30/div;        
for ii=1:size(uav1X,2)-10*30/div
    plot(uav1X(1:ii),uav1Y(1:ii),'Color','#B12C21', 'linewidth',3);
    plot(uav2X(1:ii),uav2Y(1:ii),'Color','#00863D', 'linewidth',3);
    plot(uav3X(1:ii),uav3Y(1:ii),'Color','#E49600', 'linewidth',3);
    if ii>kk*ratio+base_value && kk+2<=size(uav1Msg,1)        
        kk=kk+1;
        %uav1Msg{kk}.Header.Stamp.Sec+ uav1Msg{kk}.Header.Stamp.Nsec/1e9 - timeStart
        %timeStart=uav1Msg{kk}.Header.Stamp.Sec+ uav1Msg{kk}.Header.Stamp.Nsec/1e9;
        for j=1:p.num_vert-1
            if uav1Msg{kk}.Point.Z==j
                partition_bar(1,j+1).Position=[uav1Msg{kk}.Point.Y-p.building_block_length/10,...
                    (j)*p.building_block_height,p.building_block_length/5,p.building_block_height];       
            end
            if uav2Msg{kk}.Point.Z==j
                partition_bar(2,j+1).Position=[uav2Msg{kk}.Point.Y-p.building_block_length/10,...
                    (j)*p.building_block_height,p.building_block_length/5,p.building_block_height];       
            end
            if uav3Msg{kk}.Point.Z==j
                partition_bar(3,j+1).Position=[uav3Msg{kk}.Point.Y-p.building_block_length/10,...
                    (j)*p.building_block_height,p.building_block_length/5,p.building_block_height];       
            end
        end
        if uav1Msg{kk+1}.Point.Z~=uav1Msg{kk}.Point.Z
            base_value=base_value+2.5*30/div;
        end
    end
    drawnow
    Frame2(ii) = getframe(gcf);
    pause(0.001)
end

  % create the video writer with 1 fps
  writerObj = VideoWriter('rosvisualize.avi');
  writerObj.FrameRate = 30/div;
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


function posx=positionxfrompixel(pixel,videoxmin,videox_1_2,videox_2_4,videox_3_6)
    if pixel<videoxmin
        posx=0-(videoxmin-pixel)/videoxmin*0.3;
    elseif pixel>=videoxmin && pixel <videox_1_2
        posx=(pixel-videoxmin)/(videox_1_2-videoxmin)*1.2;
    elseif pixel>=videox_1_2 && pixel <videox_2_4
        posx=1.2+(pixel-videox_1_2)/(videox_2_4-videox_1_2)*1.2;
    elseif pixel>=videox_2_4 && pixel <videox_3_6
        posx=2.4+(pixel-videox_2_4)/(videox_3_6-videox_2_4)*1.2;
    else
        posx=3.6+(pixel-videox_3_6)/(1920-videox_3_6)*0.3;
    end
end