colours1=["#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#0072BD";"#D95319"];
fig1=figure(q)
if ishandle(fig1)
clf(fig1)
end
fig1
axis([0  la    0  lb+0.5])
fullregion=rectangle('EdgeColor','k');
fullregion.Position=[0 0    la  lb];
drawnow
%h = rectangle('EdgeColor','k','FaceColor','b');
kk=1;
for j=1:q
%     for i=1:n
%         final_workload(i)=rectangle('EdgeColor','k');
%         final_workload(i).Position=[mark(i,j,k_max1(j)), (j-1)*d, mark(i+1,j,k_max1(j))-mark(i,j,k_max1(j)),d];
%     end
%     drawnow
    for k=k_max2(j):k_max2(j)
        for li=1:n
               %partition_progress(li)=rectangle('EdgeColor','b','FaceColor',colours1(li)); 
               partition_progress(li)=rectangle('EdgeColor','b','FaceColor',colours1(li));
               partition_progress(li).Position=[mark2(li,j,k), (j-1)*d, mark2(li+1,j,k)-mark2(li,j,k),d];
        end
        kk=kk+1;
        drawnow
        %pause(0.001)
    %drawnow
    end
end
        xlabel('Length x');
        ylabel('Width y')