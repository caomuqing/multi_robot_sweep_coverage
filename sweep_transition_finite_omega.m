clear variables
close all
finite_omega_values=[1, 2, 3]*3;
n= 5;
p.sampling_period=1; %has to be one
p.transition_speed= 1.2;
p.sweeping_speed= 1;
p.building_block_length = 20;
p.building_block_height= 3;
p.num_horiz= 10;
p.num_vert= 10;
kappa= 0.05; %can be adjusted
p.building=zeros(p.num_vert, p.num_horiz);

%-----------calculate lambda2------------------
v=p.sweeping_speed/p.transition_speed;
lambda2_min=sweep_transition_lambda2_min(n,finite_omega_values,v);
%-----------calculate lambda_L-------------------
L_full=-eye(n);
for i=1:n-1
    L_full(i,i+1)=1;
end
L=L_full(1:n-1,:)*(-1);
eigenL=sort(eig(L'*L));
lambda_L=eigenL(n);
%-----------calculate eta-------------------------
omega_max=max(finite_omega_values);
eta=2*sqrt(n-1)/p.sweeping_speed^2*(4*omega_max^2+4*v*omega_max+(2-3/n)*v^2);
%-----------calculate kappa-------------------------
kappa_theoretical=4*lambda2_min/(p.sweeping_speed*eta*lambda_L);

for run=1:1
    clear m pos pos_x pos_y partition finish_sweep_time finish_transition_time;
    %-----------generate random workload profile-------------------------
    n_order=7; %create a bounded polynomial
    p1=rand(n_order+1,1)*p.building_block_length;
    p2=rand(n_order+1,1)*p.building_block_length;
    t_tot=p.num_vert*p.building_block_height;    
    idx=1;
    for tt=0:0.1:t_tot
        boundary(idx,:)=[tt,0,0];
        for i = 0:n_order
            basis_p = nchoosek(n_order, i) * (tt/t_tot)^i * (1-(tt/t_tot))^(n_order-i);
            boundary(idx,2) = boundary(idx,2) + basis_p * p1(i+1);
            boundary(idx,3) = boundary(idx,3) + basis_p * p2(i+1);
        end
        boundary(idx,3)=boundary(idx,3)+(p.num_horiz-1)*p.building_block_length;
        idx=idx+1;
    end
    for j=1:p.num_vert
        x_0(j,1)=min(boundary((j-1)*p.building_block_height/0.1+1:j*p.building_block_height/0.1+1,2));
        x_0(j,2)=max(boundary((j-1)*p.building_block_height/0.1+1:j*p.building_block_height/0.1+1,2));
        
        x_n(j,1)=min(boundary((j-1)*p.building_block_height/0.1+1:j*p.building_block_height/0.1+1,3));
        x_n(j,2)=max(boundary((j-1)*p.building_block_height/0.1+1:j*p.building_block_height/0.1+1,3));
    end
    
    for i=1:p.num_vert
        for j=1:p.num_horiz
            p.building(i,j)=randomFromArray(finite_omega_values);
        end
    end
    %---------------------zhai chao method------------------------------
    
    m=[];
    pos=[];
    partition=[];

    partition(1,1,1)=x_0(1,1);
    for i=1:n
        partition(i+1,1,1)=p.building_block_length*p.num_horiz/n*i;
    end
    partition(n+1,1,1)=x_n(1,2);
    diff=1;
    for i=1:n
        m(i,1,1)=calculateWorkload(partition(i,1,1),partition(i+1,1,1),1,p)/...
            (p.sweeping_speed*p.sampling_period);
    end
    while (diff>0.0001)
        diff=0;
        for i=1:n-1
            partition(i+1,1,1)=partition(i+1,1,1)+kappa*(m(i+1,1,1)-m(i,1,1));
        end
        for i=1:n
            m(i,1,1)=calculateWorkload(partition(i,1,1),partition(i+1,1,1),1,p)/...
                (p.sweeping_speed*p.sampling_period);
        end
        for i=1:n-1
            diff=diff+abs(m(i,1,1)-m(i+1,1,1));
        end
    end
    pos(1:n,1,1)=partition(1:n,1,1);
    
    

    kmax=ones(j);    
    finish_transition=ones(n);
    finish_sweep=zeros(n);
    finish_sweep_time=zeros(n,p.num_vert);
    finish_transition_time=zeros(n,p.num_vert);
    finish_transition_time(:,1)=ones(n,1);
    for j=1:p.num_vert
        k=1;
        if j~=1
            pos(1:n,j,1)=pos(1:n,j-1,kmax(j));
        end
        if j~=p.num_vert
            partition(1:n+1,j+1,1)=partition(1:n+1,j,kmax(j));
            partition(1,j+1,1)=x_0(j+1,1);
            partition(n+1,j+1,1)=x_n(j+1,2);
            for i=1:n %calculate workload for first time
                if rem(j,2)==1
                    m(i,j+1,k)=...
                        calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                        (p.sweeping_speed*p.sampling_period);
                else
                    m(i,j+1,k)=...
                        calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                        (p.sweeping_speed*p.sampling_period);  
                end
            end   
        end
        for i=1:n
            states(i,j,k)=1; %transition
        end
        while (prod(finish_sweep)~=1) %sweeping/transitioning/update loop
            k=k+1;
            for i=1:n
                posChange=0;
                if ~finish_transition(i)
                    states(i,j,k)=1;
                    if (rem(j,2)==1)
                        if abs(pos(i,j,k-1)-partition(i,j,kmax(j)))>p.sampling_period*p.transition_speed
                            posChange=sign(partition(i,j,kmax(j))-pos(i,j,k-1))*...
                                p.sampling_period*p.transition_speed;
                        else
                            finish_transition(i)=true;
                            [posChange1,finish_sweep(i)]=calculatePosChange(partition(i,j,kmax(j)),...
                                (p.sampling_period*p.transition_speed-abs(pos(i,j,k-1)-partition(i,j,kmax(j))))...
                                /p.transition_speed,j,p,partition(i+1,j,kmax(j)));
                            posChange=(partition(i,j,kmax(j))-pos(i,j,k-1))+posChange1;
                        end
                    else
                        if abs(pos(i,j,k-1)-partition(i+1,j,kmax(j)))>p.sampling_period*p.transition_speed
                            posChange=sign(partition(i+1,j,kmax(j))-pos(i,j,k-1))*...
                                p.sampling_period*p.transition_speed;
                        else
                            finish_transition(i)=true;
                            [posChange1,finish_sweep(i)]=calculatePosChange(partition(i+1,j,kmax(j)),...
                                (p.sampling_period*p.transition_speed-abs(pos(i,j,k-1)-partition(i+1,j,kmax(j))))...
                                /p.transition_speed,j,p,partition(i,j,kmax(j)));
                            posChange=(partition(i+1,j,kmax(j))-pos(i,j,k-1))+posChange1;
                        end
                    end
                    if finish_transition(i)
                        finish_transition_time(i,j)=k;
                        states(i,j,k)=2;
                    end
                elseif ~finish_sweep(i)
                    states(i,j,k)=2;
                    if rem(j,2)==1
                        [posChange,finish_sweep(i)]=calculatePosChange(pos(i,j,k-1),...
                            p.sampling_period,j,p,partition(i+1,j,kmax(j)));
                    else
                        [posChange,finish_sweep(i)]=calculatePosChange(pos(i,j,k-1),...
                            p.sampling_period,j,p,partition(i,j,kmax(j)));
                    end
                    if finish_sweep(i)
                        finish_sweep_time(i,j)=k;
                        states(i,j,k)=3;
                    end
                else
                    states(i,j,k)=3;
                end
                pos(i,j,k)=pos(i,j,k-1)+posChange;
                if j==p.num_vert&&states(i,j,k)==3
                	states(i,j,k)=5;
                end
            end
            if j~=p.num_vert %as long as it is not the last stripe
        %             if rem(j,2)==1
                    for i=1:n-1 %update partition mark
                        partition(i+1,j+1,k)=partition(i+1,j+1,k-1)+kappa*(m(i+1,j+1,k-1)-m(i,j+1,k-1));
                    end
                    partition(1,j+1,k)=partition(1,j+1,k-1);
                    partition(n+1,j+1,k)=partition(n+1,j+1,k-1);
        %             else
        %                 for i=2:n
        %                     partition(i,j+1,k)=partition(i,j+1,k-1)+kappa*(m(i,j+1,k-1)-m(i-1,j+1,k-1));
        %                 end
        %             end
                for i=1:n %calculate workload
                    if rem(j,2)==1
                        m(i,j+1,k)=...
                            calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                            (p.sweeping_speed*p.sampling_period);
                    else
                        m(i,j+1,k)=...
                            calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                            (p.sweeping_speed*p.sampling_period);    
                    end
                end        
            end
        end
        kmax(j+1)=k;
        finish_transition=zeros(n);
        finish_sweep=zeros(n);
    end

    k=0;
    k_vert=ceil(p.building_block_height/(p.transition_speed*p.sampling_period));   
    for j=1:p.num_vert
        for kk=1:kmax(j+1)
            pos_x(1:n,k+kk)=pos(1:n,j,kk);
            pos_y(1:n,k+kk)=(j-0.5)*p.building_block_height;
            states_all(1:n,k+kk)=states(1:n,j,kk);
        end
        if j~=p.num_vert
            for kk=1:k_vert-1 %ascending
                pos_x(1:n,k+kmax(j+1)+kk)=pos_x(1:n,k+kmax(j+1));
                pos_y(1:n,k+kmax(j+1)+kk)=(j-0.5)*p.building_block_height+kk*p.transition_speed*p.sampling_period; 
                states_all(1:n,k+kmax(j+1)+kk)=4;
                if kk==k_vert
                    pos_y(1:n,k+kmax(j+1)+kk)=(j+0.5)*p.building_block_height;
                end
            end
            k=k+kmax(j+1)+k_vert-1;
        else
            k=k+kmax(j+1);
        end
    end

    % time_used=zeros(p.num_vert);
    % for j=1:p.num_vert
    %     time_used(j)=ceil(max(m(1:n,j,kmax(j)))/(p.sampling_period*p.sweeping_speed));
    % end
    % time_used_proposed=sum(time_used(1:p.num_vert))+(p.num_vert-1)*k_vert; %vertical movement time
    %fprintf('Total time steps used is %d\n',k);
    time_previous(run)=k;
    
    %---------------------proposed method------------------------------
    clear pos_x pos_y
    m=[];
    pos=[];
    partition=[];

    partition(1,1,1)=x_0(1,1);
    for i=1:n
        partition(i+1,1,1)=p.building_block_length*p.num_horiz/n*i;
    end
    partition(n+1,1,1)=x_n(1,2);
    diff=1;
    for i=1:n
        m(i,1,1)=calculateWorkload(partition(i,1,1),partition(i+1,1,1),1,p)/...
            (p.sweeping_speed*p.sampling_period);
    end
    while (diff>0.0001)
        diff=0;
        for i=1:n-1
            partition(i+1,1,1)=partition(i+1,1,1)+kappa*(m(i+1,1,1)-m(i,1,1));
        end
        for i=1:n
            m(i,1,1)=calculateWorkload(partition(i,1,1),partition(i+1,1,1),1,p)/...
                (p.sweeping_speed*p.sampling_period);
        end
        for i=1:n-1
            diff=diff+abs(m(i,1,1)-m(i+1,1,1));
        end
    end
    pos(1:n,1,1)=partition(1:n,1,1);
    
    

    kmax=ones(j);    
    finish_transition=ones(n);
    finish_sweep=zeros(n);
    finish_sweep_time=zeros(n,p.num_vert);
    finish_transition_time=zeros(n,p.num_vert);
    finish_transition_time(:,1)=ones(n,1);
    for j=1:p.num_vert
        k=1;
        if j~=1
            pos(1:n,j,1)=pos(1:n,j-1,kmax(j));
        end
        if j~=p.num_vert
            partition(1:n+1,j+1,1)=partition(1:n+1,j,kmax(j));
            partition(1,j+1,1)=x_0(j+1,1);
            partition(n+1,j+1,1)=x_n(j+1,2);
            for i=1:n %calculate workload for first time
                if rem(j,2)==1
                    m(i,j+1,k)=abs(partition(i+1,j+1,k)-partition(i+1,j,kmax(j)))/...
                        (p.transition_speed*p.sampling_period)+...
                        calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                        (p.sweeping_speed*p.sampling_period);
                else
                    m(i,j+1,k)=abs(partition(i,j+1,k)-partition(i,j,kmax(j)))/...
                        (p.transition_speed*p.sampling_period)+...
                        calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                        (p.sweeping_speed*p.sampling_period);  
                end
            end   
        end
        for i=1:n
            states(i,j,k)=1; %transition
        end
        while (prod(finish_sweep)~=1) %sweeping/transitioning/update loop
            k=k+1;
            for i=1:n
                posChange=0;
                if ~finish_transition(i)
                    states(i,j,k)=1;
                    if (rem(j,2)==1)
                        if abs(pos(i,j,k-1)-partition(i,j,kmax(j)))>p.sampling_period*p.transition_speed
                            posChange=sign(partition(i,j,kmax(j))-pos(i,j,k-1))*...
                                p.sampling_period*p.transition_speed;
                        else
                            finish_transition(i)=true;
                            [posChange1,finish_sweep(i)]=calculatePosChange(partition(i,j,kmax(j)),...
                                (p.sampling_period*p.transition_speed-abs(pos(i,j,k-1)-partition(i,j,kmax(j))))...
                                /p.transition_speed,j,p,partition(i+1,j,kmax(j)));
                            posChange=(partition(i,j,kmax(j))-pos(i,j,k-1))+posChange1;
                        end
                    else
                        if abs(pos(i,j,k-1)-partition(i+1,j,kmax(j)))>p.sampling_period*p.transition_speed
                            posChange=sign(partition(i+1,j,kmax(j))-pos(i,j,k-1))*...
                                p.sampling_period*p.transition_speed;
                        else
                            finish_transition(i)=true;
                            [posChange1,finish_sweep(i)]=calculatePosChange(partition(i+1,j,kmax(j)),...
                                (p.sampling_period*p.transition_speed-abs(pos(i,j,k-1)-partition(i+1,j,kmax(j))))...
                                /p.transition_speed,j,p,partition(i,j,kmax(j)));
                            posChange=(partition(i+1,j,kmax(j))-pos(i,j,k-1))+posChange1;
                        end
                    end
                    if finish_transition(i)
                        finish_transition_time(i,j)=k;
                        states(i,j,k)=2;
                    end
                elseif ~finish_sweep(i)
                    states(i,j,k)=2;
                    if rem(j,2)==1
                        [posChange,finish_sweep(i)]=calculatePosChange(pos(i,j,k-1),...
                            p.sampling_period,j,p,partition(i+1,j,kmax(j)));
                    else
                        [posChange,finish_sweep(i)]=calculatePosChange(pos(i,j,k-1),...
                            p.sampling_period,j,p,partition(i,j,kmax(j)));
                    end
                    if finish_sweep(i)
                        finish_sweep_time(i,j)=k;
                        states(i,j,k)=3;
                    end
                else
                    states(i,j,k)=3;
                end
                pos(i,j,k)=pos(i,j,k-1)+posChange;
                if j==p.num_vert&&states(i,j,k)==3
                	states(i,j,k)=5;
                end
            end
            if j~=p.num_vert %as long as it is not the last stripe
        %             if rem(j,2)==1
                    for i=1:n-1 %update partition mark
                        partition(i+1,j+1,k)=partition(i+1,j+1,k-1)+kappa*(m(i+1,j+1,k-1)-m(i,j+1,k-1));
                    end
                    partition(1,j+1,k)=partition(1,j+1,k-1);
                    partition(n+1,j+1,k)=partition(n+1,j+1,k-1);
        %             else
        %                 for i=2:n
        %                     partition(i,j+1,k)=partition(i,j+1,k-1)+kappa*(m(i,j+1,k-1)-m(i-1,j+1,k-1));
        %                 end
        %             end
                for i=1:n %calculate workload
                    if rem(j,2)==1
                        m(i,j+1,k)=abs(partition(i+1,j+1,k)-partition(i+1,j,kmax(j)))/...
                            (p.transition_speed*p.sampling_period)+...
                            calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                            (p.sweeping_speed*p.sampling_period);
                    else
                        m(i,j+1,k)=abs(partition(i,j+1,k)-partition(i,j,kmax(j)))/...
                            (p.transition_speed*p.sampling_period)+...
                            calculateWorkload(partition(i,j+1,k),partition(i+1,j+1,k),j+1,p)/...
                            (p.sweeping_speed*p.sampling_period);    
                    end
                end        
            end
        end
        kmax(j+1)=k;
        finish_transition=zeros(n);
        finish_sweep=zeros(n);
    end

    k=0;
    k_vert=ceil(p.building_block_height/(p.transition_speed*p.sampling_period*0.5));   
    for j=1:p.num_vert
        for kk=1:kmax(j+1)
            pos_x(1:n,k+kk)=pos(1:n,j,kk);
            pos_y(1:n,k+kk)=(j-0.5)*p.building_block_height;
            states_all(1:n,k+kk)=states(1:n,j,kk);
        end
        if j~=p.num_vert
            for kk=1:k_vert-1 %ascending
                pos_x(1:n,k+kmax(j+1)+kk)=pos_x(1:n,k+kmax(j+1));
                pos_y(1:n,k+kmax(j+1)+kk)=(j-0.5)*p.building_block_height+kk*p.transition_speed*p.sampling_period*0.5; 
                states_all(1:n,k+kmax(j+1)+kk)=4;
                if kk==k_vert
                    pos_y(1:n,k+kmax(j+1)+kk)=(j+0.5)*p.building_block_height;
                end
            end
            k=k+kmax(j+1)+k_vert-1;
        else
            k=k+kmax(j+1);
        end
    end

    % time_used=zeros(p.num_vert);
    % for j=1:p.num_vert
    %     time_used(j)=ceil(max(m(1:n,j,kmax(j)))/(p.sampling_period*p.sweeping_speed));
    % end
    % time_used_proposed=sum(time_used(1:p.num_vert))+(p.num_vert-1)*k_vert; %vertical movement time
    %fprintf('Total time steps used is %d\n',k);
    time_proposed(run)=k;
    
    %--------------------calculate time for other methods------------------
    num_turns=ceil(p.num_vert/n);
    width=p.num_vert*p.building_block_height/(num_turns*n);
            
    work_new_bound=ones(num_turns*n,2)*p.num_horiz*p.building_block_length/2;
    for j=1:num_turns*n
        work_new(j)=0;
        current_height=(j-1)*width;
        current_j=floor(current_height/p.building_block_height);
        for h=current_j+1:p.num_vert
                if j*width<=h*p.building_block_height+0.00001
                    work_new(j)=work_new(j)+(j*width-current_height)/p.building_block_height*...
                        calculateWorkload(x_0(h,1),x_n(h,2),h,p);
                    if work_new_bound(j,1)>=x_0(h,1)
                        work_new_bound(j,1)=x_0(h,1);
                    end
                    if work_new_bound(j,2)<=x_n(h,2);
                        work_new_bound(j,2)=x_n(h,2);
                    end
                    break;
                else
                    work_new(j)=work_new(j)+(h*p.building_block_height-current_height)/p.building_block_height*...
                        calculateWorkload(x_0(h,1),x_n(h,2),h,p);
                    current_height=h*p.building_block_height;
                    if work_new_bound(j,1)>=x_0(h,1)
                        work_new_bound(j,1)=x_0(h,1);
                    end
                    if work_new_bound(j,2)<=x_n(h,2);
                        work_new_bound(j,2)=x_n(h,2);
                    end
                end
        end
    end
    for i=1:p.num_vert
        work(i)=calculateWorkload(x_0(i,1),x_n(i,2),i,p);
    end
    if sum(work_new)<sum(work)-0.000001||sum(work_new)>sum(work)+0.000001
        fprintf('Something is wrong.\n');
        fprintf('work_new sum: %f.\n',sum(work_new));
        fprintf('work total: %f.\n',sum(work));
    end
        
    k_vert_new=ceil(width/(p.transition_speed*p.sampling_period));
    for i=1:n
        time_used_cell(i)=0;
        for ii=1:num_turns
            time_used_cell(i)= time_used_cell(i)+ ceil(work_new(num_turns*(i-1)+ii)...
                /(p.sampling_period*p.sweeping_speed));
            if ii~=1
                time_used_cell(i)= time_used_cell(i)+ ...
                    abs(work_new_bound(num_turns*(i-1)+ii,2)-work_new_bound(num_turns*(i-1)+ii-1,2))/...
                    (p.transition_speed*p.sampling_period);
            end
        end
        time_used_cell(i)= time_used_cell(i)+k_vert_new*(num_turns-1);
    end
    %fprintf('Total time steps used using equal cell decomposition %d\n',max(time_used_cell));    
    time_cell_decomp(run)=max(time_used_cell);
    for i=1:n
        time_used_line(i)=0;
        for ii=1:num_turns
            time_used_line(i)= time_used_line(i)+ ceil(work_new(n*(ii-1)+i)...
                /(p.sampling_period*p.sweeping_speed));
            if ii~=1
                time_used_cell(i)= time_used_cell(i)+ ...
                    abs(work_new_bound(n*(ii-1)+i,2)-work_new_bound(n*(ii-2)+i,2))/...
                    (p.transition_speed*p.sampling_period);
            end            
        end
        time_used_line(i)= time_used_line(i)+k_vert_new*(num_turns-1)*n;
%         time_used_line(i)=ceil(work_new(i)...
%             /(p.sampling_period*p.sweeping_speed))+...
%             ceil(work_new(i+n)...
%             /(p.sampling_period*p.sweeping_speed))+n*k_vert_new;
    end
    %fprintf('Total time steps used using line formation %d\n',max(time_used_line)); 
    time_line_formation(run)=max(time_used_line);
    for j=1:1
        sumofworkload=calculateWorkload(x_0(j,1),x_n(j,2),j,p)/...
            (p.sweeping_speed*p.sampling_period);
        if sum(m(1:n,j,kmax(j)))<sumofworkload-0.000001||sum(m(1:n,j,kmax(j)))>sumofworkload+0.000001
            fprintf('workload sum does not amount to the total. Something is wrong.\n');
            fprintf('workload sum: %f.\n',sum(m(1:n,j,kmax(j))));
            fprintf('workload total: %f.\n',calculateWorkload(0,p.building_block_length*p.num_horiz,j,p)/...
                (p.sweeping_speed*p.sampling_period));
        end
    end
end
 final_result=[time_proposed; time_previous; time_cell_decomp; time_line_formation]';
 
fprintf('Avg time used for proposed methods is %f\n',mean(time_proposed));
fprintf('Avg time used for previous methods is %f\n',mean(time_previous));

fprintf('Avg time steps used using cell decomposition %f\n',mean(time_cell_decomp));    
fprintf('Avg time steps used using line formation %f\n',mean(time_line_formation)); 

function [posChange,finish_current]=calculatePosChange(currentPos, time, j, p, bar)
    if j>p.num_vert
        posChange=0;
        finish_current=true;
        return;
    end
    initialPos=currentPos;
    remaining_effort=time*p.sweeping_speed;
    num=constrain(floor(currentPos/p.building_block_length), 0, p.num_horiz-1);
    if (rem(j,2)==1)
        while remaining_effort>0
            if ((num+1)*p.building_block_length-currentPos)*p.building(p.num_vert-j+1, num+1)...
                    >remaining_effort
                currentPos = currentPos + remaining_effort/p.building(p.num_vert-j+1, num+1);
                remaining_effort=0;
            else
                remaining_effort=remaining_effort-...
                    ((num+1)*p.building_block_length-currentPos)*p.building(p.num_vert-j+1, num+1);
                currentPos=(num+1)*p.building_block_length;
            end
            if currentPos>=bar
                currentPos=bar;
                finish_current=true;
                break
            else
                finish_current=false;
            end
            num=num+1;
        end
    else
        while remaining_effort>0
            if (currentPos-(num*p.building_block_length))*p.building(p.num_vert-j+1, num+1)...
                    >remaining_effort
                currentPos = currentPos - remaining_effort/p.building(p.num_vert-j+1, num+1);
                remaining_effort=0;
            else
                remaining_effort=remaining_effort-...
                    (currentPos-num*p.building_block_length)*p.building(p.num_vert-j+1, num+1);
                currentPos=num*p.building_block_length;
            end
            if currentPos<=bar
                currentPos=bar;
                finish_current=true;
                break
            else
                finish_current=false;
            end
            num=num-1;
        end
    end
    posChange=currentPos-initialPos;                
end



function workload=calculateWorkload(xa,xb,j,p)
    blk_num_left=constrain(floor(xa/p.building_block_length), 0, p.num_horiz-1);
    blk_num_right=constrain(floor(xb/p.building_block_length), 0, p.num_horiz-1);
    workload=0;
    for num=blk_num_left:blk_num_right
        if xb>(num+1)*p.building_block_length
            workload=workload+((num+1)*p.building_block_length-xa)*p.building(p.num_vert-j+1, num+1);
            xa=(num+1)*p.building_block_length;
        else
            workload=workload+(xb-xa)*p.building(p.num_vert-j+1, num+1);
            break;
        end
    end
end
            
function y=constrain(x, a, b)
    if x<a
        y=a;
    elseif x>b
        y=b;
    else
        y=x;
    end
end

function rand =randomFromArray(x)
    msize = numel(x);
    idx = randperm(msize);
    rand=x(idx(1));
end