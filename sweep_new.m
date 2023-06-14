clear all
close all

plot_figure=false;
run_incremental_algo=true;
run_normal_algo=true;
pos_plot=true;
%common for both methods
n=5;    %number of agents
la=6   %sweep region width
lb=20;   %sweep region length
d=1;    %agent effective width
q=ceil(lb/d);   %number of stripes
rou_low=1; %lowest workload density
rou_up=5   % max workload density
v=0.2 %sweeping rate (workload/time)
Ts=0.3 % sampling period

kappa_init=0.5;
kappa_end=0.5;
kappa=0.1;

Omega=[1,0,0,0;-1,1,0,0;0,-1,1,0;0,0,-1,1;0,0,0,-1]*d*rou_low;
Omega_G=[1,0,0,0;-1,1,0,0;0,-1,1,0;0,0,-1,1;0,0,0,-1]*d*rou_up;
L=[-1,1,0,0,0;0,-1,1,0,0;0,0,-1,1,0;0,0,0,-1,1]*(-1); % communication graph, depends on
%number of agents
%L=[-3,1,1,1;1,-3,1,1;1,1,-3,1]*(-1); % communication graph
%L=[-2,1,1,0;0,-2,1,1;0,1,1,-2]*(-1); % communication graph

P=Omega*L;
Pu=(P+P')/2;
eigen=eig(Pu);
lambda2=eigen(2);
alpha=rou_up/rou_low;
max_time=q*d*rou_up*la/(Ts*v);
%specific for proposed method
eigenL=eig(L'*L);
lambda_L=eigenL(4);
lx=la;
omega_up=d*rou_up;
omega_low=d*rou_low;
etas(1)=2*norm(Omega'*Omega);
etas(2)=2*norm(Omega'*Omega_G);
etas(3)=2*norm(Omega_G'*Omega);
etas(4)=2*norm(Omega_G'*Omega_G);
eta=8*sqrt(n-1)*omega_up^2; %max lipschitz constant based on matrix norm
%inequality
%eta=max(etas); %max lipschitz constant based on sampling
K_compare=zeros(100,2);
k=1;

reached_final=zeros(n);
pos=zeros(n,10^5);
pos_level=zeros(n,10^5);
path_left=zeros(n,10^5);
pos_level(:,1)=ones(n,1);
bars=zeros(n+1,q,10^5);

for i=1:n     %set initial bars position
    for j=1:q
        bars(i+1,j,1)= i*la/n;
        bars(1,j,1)=0;
    end
    pos(i,1)=(i-1)*la/n;
end

syms x y
%rou=@(x,y) (rou_up-rou_low)/2*sin(x+y)+(rou_up+rou_low)/2; 
rou=@(x,y) (rou_up-rou_low)/2*sin(x.^2/10+(y).^2/10)+(rou_up+rou_low)/2;

while (prod(reached_final)~=1)
    k=k+1; 
    for i=1:n
       x2=sym('x2');
%        f1=(rou_up+rou_low)*(pos(i,k-1)-x2)*(-d)+(rou_up-rou_low)*...
%           (-sin(pos(i,k-1)+(pos_level(i,k-1)-1)*d)+sin(x2+(pos_level(i,k-1)-1)*d)+...
%           sin(pos(i,k-1)+pos_level(i,k-1)*d)-sin(x2+pos_level(i,k-1)*d))==2*Ts*v;
       f1=int(int(rou,x,pos(i,k-1),x2),y,(pos_level(i,k-1)-1)*d,pos_level(i,k-1)*d)==Ts*v;
       pos(i,k)=vpasolve(f1,x2,pos(i,k-1)+Ts*v*2/((rou_up+rou_low)*d));
       if pos(i,k)>=bars(i+1,pos_level(i,k-1),k-1)
           if (pos_level(i,k-1))~=q
               pos_level(i,k)=pos_level(i,k-1)+1;
               pos(i,k)=bars(i,pos_level(i,k),k-1);    %uplevel
           else %last stripe already
               pos(i,k)=bars(i+1,pos_level(i,k-1),k-1);
               pos_level(i,k)=pos_level(i,k-1);
               reached_final(i)=1;
           end
       else
           pos_level(i,k)=pos_level(i,k-1);
       end

    end
    
    for i=1:n
        for j=pos_level(i,k):q
            if j>pos_level(i,k)
                path_left(i,k)=path_left(i,k)+bars(i+1,j,k-1)-bars(i,j,k-1);
            else
                path_left(i,k)=bars(i+1,pos_level(i,k),k-1)-pos(i,k);
            end
        end
    end
    
    bars(:,:,k)=bars(:,:,k-1);
    for i=1:n-1
        adjust_level=max(pos_level(i,k),pos_level(i+1,k))+1;
        if adjust_level<=q
            bars(i+1,adjust_level,k)=bars(i+1,adjust_level,k)+kappa*...
                (path_left(i+1,k)-path_left(i,k));
            if bars(i+1,adjust_level,k)<bars(i,adjust_level,k)
                bars(i+1,adjust_level,k)=bars(i,adjust_level,k);
            end
        end
    end
   
    parfor i=1:n
        actual_workload(i,k)=0;
        for j=1:q
            actual_workload(i,k)=actual_workload(i,k)+...
                integral2(rou, bars(i,j,k), bars(i+1,j,k), (j-1)*d, j*d);
        end
    end
            
end
    
k_max=k;
str=["$m_1(k)$";"$m_2(k)$";"$m_3(k)$";"$m_4(k)$";"$m_5(k)$";"$m_6(k)$";"$m_7(k)$";"$m_8(k)$"];

figure
xlabel('Sampling Time k');
ylabel('Workload')
for i=1:n
    plot(1:k_max,squeeze(actual_workload(i,1:k_max)),'LineWidth',2)
    legend({str(i)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
    hold on
end

figure
xlabel('Sampling Time k');
ylabel('path left')
for i=1:n
    plot(1:k_max,squeeze(path_left(i,1:k_max)),'LineWidth',2)
    hold on
    legend({str(1),str(2),str(3),str(4),str(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')

end

