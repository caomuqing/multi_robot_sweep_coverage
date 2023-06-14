% fig=figure(10);
% if ishandle(fig)
% clf(fig);
% end
colours1=["#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];

figure

partition_bar=gobjects(n,p.num_vert);
for i=1:n
    for j=1:p.num_vert
        partition_bar(i,j)=rectangle('EdgeColor','b','FaceColor',colours1(i));
        hold on
    end
end
for i=1:n
    for j=1:p.num_vert
        %partition_bar(i,1)=rectangle('EdgeColor','g','FaceColor','g');
        partition_bar(i,j).Position=[partition(i,j,kmax(j)),(j-1)*p.building_block_height,...
            partition(i+1,j,kmax(j))-partition(i,j,kmax(j)),p.building_block_height];
    end
end
drawnow
