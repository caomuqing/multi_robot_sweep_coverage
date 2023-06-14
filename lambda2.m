clear all
close all
n=5;
finite_density_value=true;
step_size=1.5;
omega_low=1;
omega_high=10;
finite_omega_values=[1, 1.5, 2];

if finite_density_value
    steps=size(finite_omega_values, 2)
else
    steps=ceil((omega_high-omega_low)/step_size)+1;
end

for i=1:2^(n-1)-1
    quotient=i;
    for j=1:n-1
        bit(j)=rem(quotient,2);
        quotient=fix(quotient/2);
    end
    binary(i,:)=-bit;
end

num_valid_graph=1;
for counter=1:(2^(n-1)-1)^(n-1)
    quotient=counter;
    for j=1:n-1
        bitmask=rem(quotient,2^(n-1)-1)+1;
        quotient=fix(quotient/(2^(n-1)-1));
        
        ii=1;
        for i=1:n
            if i~=j
                graph(j,i)=binary(bitmask,ii);
                ii=ii+1;
            end
        end
        graph(j,j)=-sum(binary(bitmask,:));
    end
    if rank(graph)==n-1
        valid_graphs(:,:,num_valid_graph)=graph;
        num_valid_graph=num_valid_graph+1;
    end
end
num_valid_graph=num_valid_graph-1;
valid_graphs_reshape=reshape(valid_graphs,n-1,n,num_valid_graph);

%for testing only
valid_graphs(:,:,1) = [1,-1,0,0,0;
                        0,1,-1,0,0;
                        0,0,1,-1,0;
                        0,0,0,1,-1];
% valid_graphs(:,:,1) = [1,0,-1,0,0;
%                         0,1,-1,0,0;
%                         -1,-1,3,-1,0;
%                         0,0,-1,2,-1];
% valid_graphs(:,:,1) = [1,0,0,0,-1;
%                         0,1,0,0,-1;
%                         0,0,1,0,-1;
%                         0,0,0,1,-1];
num_valid_graph=1; %only one graph is defined for this test

Omg_bare_temp=eye(n,n);    %define Omega(k)
for j=1:n-1
    Omg_bare_temp(j+1,j)=-1;
end
Omg_bare=Omg_bare_temp(:,1:n-1);        
      
omega_matrix=zeros(n-1,n-1);

if finite_density_value
    for counter=1:(steps)^(n-1)
        quotient=counter;
        for j=1:n-1
            bitmask=rem(quotient,steps);
            quotient=fix(quotient/steps);
            omega_matrix(j,j)=finite_omega_values(bitmask+1);
        end
        omega_matrix_big(:,:,counter)=omega_matrix;
    end   
else
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
end
num_workload_dist=counter;

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
        Omega=Omg_bare*omega_matrix_big(:,:,j)*valid_graphs(:,:,i);
        eigen_omega=sort(eig(0.5*Omega'+0.5*Omega));
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
        end
    end
    if this_graph_is_good
        good_graph_counter=good_graph_counter+1;
        good_graph_index(good_graph_counter)=i;
    end
end

if this_graph_is_good
    fprintf('This is a good graph!')
else
    fprintf('This graph does not guarantee lambda2 > 0')
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