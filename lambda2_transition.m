clear all
close all
n=3;
step_size=5;
omega_low=5;
omega_high=10;
omega_transition=4.9;

steps=ceil((omega_high-omega_low)/step_size)+1;


% only the directed chain graph is tested
L_full=-eye(n);
for i=1:n-1
    L_full(i,i+1)=1;
end
L=L_full(1:n-1,:)*(-1);
valid_graphs(:,:,1)=L;
%for testing only
% valid_graphs(:,:,1) = [1,-1,0,0,0;
%                         0,1,-1,0,0;
%                         0,0,1,-1,0;
%                         0,0,0,1,-1];
% valid_graphs(:,:,1) = [1,0,-1,0,0;
%                         0,1,-1,0,0;
%                         -1,-1,3,-1,0;
%                         0,0,-1,2,-1];
% valid_graphs(:,:,1) = [1,0,0,0,-1;
%                         0,1,0,0,-1;
%                         0,0,1,0,-1;
%                         0,0,0,1,-1];
num_valid_graph=1;

Omg_bare_temp=eye(n,n);    %define Omega(k)
for j=1:n-1
    Omg_bare_temp(j+1,j)=-1;
end
Omg_bare=Omg_bare_temp(:,1:n-1);        
      
omega_matrix=zeros(n-1,n-1);
for counter=1:(steps)^(n-1)
    quotient=counter;
    for j=1:n-1
        bitmask=rem(quotient,steps);
        quotient=fix(quotient/steps);
        omega_matrix(j,j)=omega_low+bitmask*step_size;
        if omega_matrix(j,j)>omega_high
            omega_matrix(j,j)=omega_high;
        end
    end
    omega_matrix_big(:,:,counter)=omega_matrix;
end
num_workload_dist=counter;

I_matrix=zeros(n,n-1);
for counter=1:2^(n-1)
    quotient=counter;
    for j=1:n-1
        bitmask=rem(quotient,2);
        quotient=fix(quotient/2);
        if bitmask==1
            I_matrix(j,j)=omega_transition;
        else
            I_matrix(j,j)=-omega_transition;
        end
    end
    I_matrix_big(:,:,counter)=I_matrix;
end
num_I_matrix=counter;

counter_lambda2_not_larger_zero=0;
lambda2_min=10;
lambda2_max_vec=ones(num_workload_dist,1)*0;
lambda2_max_vec_index=ones(num_workload_dist,1);

good_graph_counter=0;
for i=1:num_valid_graph
    this_graph_is_good=1;
    for j=1:num_workload_dist
        %graph_matrix(:,:)=valid_graphs(i,:,:);
        %omega_matrix=omega_matrix_big(j,:,:);
        for k=1:num_I_matrix
            Omega=(Omg_bare*omega_matrix_big(:,:,j)+I_matrix_big(:,:,k))*valid_graphs(:,:,i);
            eigen_omega=sort(eig(0.5*Omega.'+0.5*Omega));
            lambda2_temp=eigen_omega(2);
            if lambda2_temp<lambda2_min
                lambda2_min=lambda2_temp;
                lambda2_min_index=[i,j];
            end
            if lambda2_temp>lambda2_max_vec(j)
                lambda2_max_vec(j)=lambda2_temp;
                lambda2_max_vec_index(j)=i;
            end
            if lambda2_temp<=0&&isreal(lambda2_temp)
                counter_lambda2_not_larger_zero=counter_lambda2_not_larger_zero+1;
                lambda2_index(:,:,counter_lambda2_not_larger_zero)=[i,j];
                this_graph_is_good=0;
                Omega
                return
            end
        end
    end
end

% a_graph=     [1     0    -1     0     0;
%            0     1    -1     0     0;
%            0     0     1    -1     0;
%            0     0     0     1    -1];
% syms ome1 ome2 ome3 ome4 ome5 ome6 ome7 ome8 ome9
%        omega_vector=[ome1 0 0 0
%            0 ome2 0 0
%            0 0 ome3 0
%            0 0 0 ome4];
%        omega_vector=[10 0 0 0
%            0 1 0 0
%            0 0 3 0
%            0 0 0 2];
%       0.5*Omg_bare* omega_vector*a_graph+0.5*a_graph'* omega_vector*Omg_bare'