clear;
sampling_time = 0.5;
bag = rosbag('sweep4.bag');
bagselect1 = select(bag, 'Topic', '/uav1/workload');
bagselect2 = select(bag, 'Topic', '/uav2/workload');
bagselect3 = select(bag, 'Topic', '/uav3/workload');
uav1Msg = readMessages(bagselect1);
uav2Msg = readMessages(bagselect2);
uav3Msg = readMessages(bagselect3);
timeStart = uav1Msg{1}.Header.Stamp.Sec+ uav1Msg{1}.Header.Stamp.Nsec/1e9;

start_time=52
end_time=107
workloadinfo=[];
timeStart = uav1Msg{start_time}.Header.Stamp.Sec+ uav1Msg{start_time}.Header.Stamp.Nsec/1e9;
for i=1:end_time-start_time
    workloadinfo(i,1) = uav1Msg{start_time+i-1}.Header.Stamp.Sec+ uav1Msg{start_time+i-1}.Header.Stamp.Nsec/1e9 - timeStart;
    workloadinfo(i,2) = uav1Msg{start_time+i-1}.Point.X;
    workloadinfo(i,3) = uav2Msg{start_time+i-1}.Point.X;
    workloadinfo(i,4) = uav3Msg{start_time+i-1}.Point.X;
end

str3=["$\psi^q_1(k)$";"$\psi^q_2(k)$";"$\psi^q_3(k)$";"$\psi^q_4(k)$";"$\psi^q_5(k)$";"$\psi^q_6(k)$";"$\psi^q_7(k)$";"$\psi^q_8(k)$"];

% odomPt=[];
% j=0;
% for i=1:100:length(odomMsg)
%     j=j+1;
%     odomPt(j,1)=odomMsg{i}.Header.Stamp.Sec+ odomMsg{i}.Header.Stamp.Nsec/1e9 - timeStart;
%     odomPt(j,2)=odomMsg{i}.Pose.Pose.Position.X;
%     odomPt(j,3)=odomMsg{i}.Pose.Pose.Position.Y;
% end
% 
% wall=[];
% for i=2:length(wallMsg)
%     wall(i,1)=wallMsg{i}.Header.Stamp.Sec+ wallMsg{i}.Header.Stamp.Nsec/1e9 - timeStart;
%     wall(i,2)=wallMsg{i}.Pose.Position.X-wallMsg{i-1}.Pose.Position.X;
%     wall(i,3)=wallMsg{i}.Pose.Position.Y-wallMsg{i-1}.Pose.Position.Y;
%     wall(i,4)=wallMsg{i}.Pose.Orientation.X;
%     wall(i,5)=wallMsg{i}.Pose.Orientation.Y;   
%     wall(i,6)=wallMsg{i}.Pose.Orientation.Z;
%     wall(i,7)=wallMsg{i}.Pose.Orientation.W;
%     wall(i,8)=wallMsg{i}.Pose.Position.Z-wallMsg{i-1}.Pose.Position.Z;
% 
% end

close all;

plot(workloadinfo(:,1)/sampling_time,workloadinfo(:,2),'LineWidth',2);
hold on
plot(workloadinfo(:,1)/sampling_time,workloadinfo(:,3),'LineWidth',2);
plot(workloadinfo(:,1)/sampling_time,workloadinfo(:,4),'LineWidth',2);
xlabel('Sampling Time k','FontSize',14);
ylabel('Operation cycle','FontSize',14)
legend({str3(1),str3(2),str3(3)},'Interpreter','latex','FontSize',16, 'Orientation','vertical')
%legend('p_x(k)','\hat{p}_{r,x}(k|k)');
% figure
% plot(odomPt(:,1)/sampling_time,odomPt(:,3),'LineWidth',2);
% hold on
% plot(trajPt(:,1)/sampling_time,trajPt(:,3),'LineWidth',2);
% xlabel('Sampling Time k');
% ylabel('Position y')
% legend({'$p_y(k)$','$\hat{p}_{r,y}(k|k)$'},'Interpreter','latex','FontSize',14)
% 
% figure
% plot(wall(:,1)/sampling_time,wall(:,4),wall(:,1)/sampling_time,wall(:,5),wall(:,1)/sampling_time,wall(:,6),wall(:,1)/sampling_time,wall(:,7),'LineWidth',2);
% xlabel('Sampling Time k');
% ylabel('Facade Coefficients Estimate')
% legend({'$\hat{a}(k)$','$\hat{b}(k)$','$\hat{c}(k)$','$\hat{d}(k)$'},'Interpreter','latex','FontSize',14, 'Orientation','horizontal')
% figure
% plot(wall(:,1)/sampling_time,wall(:,2),wall(:,1)/sampling_time,wall(:,3),wall(:,1)/sampling_time,wall(:,8),'LineWidth',2);
% xlabel('Sampling Time k');
% ylabel('Tracking Output Estimate')
% legend({'$\delta{z}_{r,x}(j|k)$','$\delta{z}_{r,y}(j|k)$','$\delta{z}_{r,z}(j|k)$'},'Interpreter','latex','FontSize',14)
