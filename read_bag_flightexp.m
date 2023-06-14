close all
clear
bag{1} = rosbag('/home/iot/msc_bags/2022-05-05-20-13-12.bag');
bag{2} = rosbag('/home/iot/msc_bags/2022-05-05-20-13-36.bag');
bag{3} = rosbag('/home/iot/msc_bags/2022-05-05-20-13-50.bag');

bagselect{1} = select(bag{1}, 'Topic', '/uav1/sweepinginfo');    
msg{1}  = readMessages(bagselect{1});
bagselect{2} = select(bag{2}, 'Topic', '/uav2/sweepinginfo');    
msg{2}  = readMessages(bagselect{2});
bagselect{3} = select(bag{3}, 'Topic', '/uav3/sweepinginfo');    
msg{3}  = readMessages(bagselect{3});

start_time = 146.41;
end_time = 159.6;
for jj=1:3
    prev_value = -100;
    count =1;
    %start_time = min(start_time, msg{jj}{1}.Header.Stamp.Sec + msg{jj}{1}.Header.Stamp.Nsec/10e8);
    for i=1:size(msg{jj},1)
        value = msg{jj}{i}.MIJplus1;
        time = msg{jj}{i}.Header.Stamp.Sec + msg{jj}{i}.Header.Stamp.Nsec/10e8;
        if value ~= prev_value && time <end_time && time >start_time
            p{jj}.value(count) = value;
            p{jj}.time(count) = time;
            count = count+1;
        end
        prev_value = value;
    end
end

figure
xlabel('Time (s)','FontSize',15);
ylabel('Workload','FontSize',15)
hold on
for jj=1:3
    plot (p{jj}.time-start_time,p{jj}.value,'-','LineWidth',1.5)
end
    legend({'$\mu^2_1(k)$','$\mu^2_2(k)$','$\mu^2_3(k)$'},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
