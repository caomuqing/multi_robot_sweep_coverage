clear all
close all
n=3;

syms o1 o2 o3 o4 o5 o6 o7 o8 o9;
ome=[o1 o2 o3 o4 o5 o6 o7 o8 o9];
syms v a1 a2 a3 a4 a5 a6 a7 a8 a9
A_opt=[a1 a2 a3 a4 a5 a6 a7 a8 a9];
%for testing only


Omg_bare_temp=eye(n,n);    %define Omega(k)
for j=1:n-1
    Omg_bare_temp(j+1,j)=-1;
end
Omg_bare=Omg_bare_temp(:,1:n-1);        
      
% omega_matrix=zeros(n-1,n-1);
omega_matrix=sym(zeros(n-1,n-1));
for i=1:n-1
    for j=1:n-1
        if i~=j
            omega_matrix(i,j)=0;
        else
            omega_matrix(j,j)=ome(j);
        end
    end
end

A_matrix=sym(zeros(n,n-1));
for i=1:n-1
    A_matrix(i,i)=A_opt(i);
end
% print result for original algorithm
% Omg_bare.'*omega_matrix.'*omega_matrix*Omg_bare

I_matrix=eye(n)-1/n*ones(n,1)*ones(1,n);
omega_A=Omg_bare*omega_matrix+v*A_matrix;
omega_A.'*I_matrix*omega_A
