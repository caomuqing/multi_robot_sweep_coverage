clear all
close all

plot_figure=true;
run_incremental_algo=true;
run_normal_algo=true;
normal_algo_pos_plot=false;
%common for both methods
n=5;    %number of agents
la=6   %sweep region width
lb=6;   %sweep region length
d=1;    %agent effective width
q=ceil(lb/d);   %number of stripes
rou_low=1; %lowest workload density
rou_up=1.5   % max workload density
v=0.5 %sweeping rate (workload/time)
Ts=0.01 % sampling period

kappa_init=8;
kappa_end=8;

N=eye(n)-ones(n,n)/n;  %average_matrix
Omg_bare_temp=eye(n,n);    %define Omega(k)
for j=1:n-1
    Omg_bare_temp(j+1,j)=-1;
end
Omg_bare=Omg_bare_temp(:,1:n-1);
Omega=Omg_bare*d*rou_low;
Omega_G=Omg_bare*d*rou_up;

Phi_bare=zeros(n,n); %construct Phi(k)
Phi_bare(n,n)=1;
for j=1:n-1
    Phi_bare(j+1,j)=1;
end
Phi_low=Phi_bare*d*rou_low;
Phi_high=Phi_bare*d*rou_up;

Gamma_bare=[Omg_bare Phi_bare];
Gamma_high=[Omega_G Phi_high];
Gamma_low=[Omega Phi_low];
eigen_gamma=sort(eig(Gamma_high*Gamma_high'));
lambda_gamma=eigen_gamma(n);

syms rho
syms rho1 rho2 rho3 rho4 rho5 rho6 rho7 rho8 rho9 rho10 rho11 rho12 rho13 rho14 rho15 rho16 rho17
rho_x=[rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13,rho14,rho15,rho16,rho17];
syms rhop1 rhop2 rhop3 rhop4 rhop5 rhop6 rhop7 rhop8 rhop9 rhop10 rhop11 rhop12 rhop13 rhop14 rhop15 rhop16 rhop17
rho_p=[rhop1,rhop2,rhop3,rhop4,rhop5,rhop6,rhop7,rhop8,rhop9,rhop10,rhop11,rhop12,rhop13,rhop14,rhop15,rhop16,rhop17];

%Gamma_rho=Gamma_bare;%*rho;
for i=1:n-1
    Gamma_rho(:,i)=Gamma_bare(:,i)*rho_x(i);
end
for i=1:n
    Gamma_rho(:,n-1+i)=Gamma_bare(:,n-1+i)*rho_p(i);
end
Gamma_rho.'*Gamma_rho;
Gamma_rho.'*N*Gamma_rho;
Gamma_rho.'*Gamma_rho;
%Omega=[1,0,0,0;-1,1,0,0;0,-1,1,0;0,0,-1,1;0,0,0,-1]*d*rou_low;
%Omega_G=[1,0,0,0;-1,1,0,0;0,-1,1,0;0,0,-1,1;0,0,0,-1]*d*rou_up;
L_full=-eye(n);
for i=1:n-1
    L_full(i,i+1)=1;
end
L=L_full(1:n-1,:)*(-1);

% L = [1,-1,0,0,0;
%     -1,2,-1,0,0;
%     0,-1,2,-1,0;
%     0,0,0,1,-1];
%L=[-1,1,0,0,0;0,-1,1,0,0;0,0,-1,1,0;0,0,0,-1,1]*(-1); % communication graph, depends on
%number of agents
%L=[-3,1,1,1;1,-3,1,1;1,1,-3,1]*(-1); % communication graph
%L=[-2,1,1,0;0,-2,1,1;0,1,1,-2]*(-1); % communication graph
omega_new_up=eye(n,n)*d*rou_up;     %%for analyzing the incremental algorithm
omega_new_low=eye(n,n)*d*rou_low;
    for j=1:n-1
        omega_new_up(j+1,j)=d*(rou_up-rou_low);
        omega_new_low(j+1,j)=d*(rou_low-rou_up);
    end

P=Omega*L;
Pu=(P+P')/2;
eigen=sort(eig(Pu));
lambda2=eigen(2);
alpha=rou_up/rou_low;
max_time=q*d*rou_up*la/(Ts*v);
%specific for proposed method
eigenL=sort(eig(L'*L));
lambda_L=eigenL(n);
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
sigma=6*sqrt(2*n-1)*omega_up^2;
b=1/(2*lambda2)+100;
a=2*sigma*lambda_L/(1/b-2*lambda2)^2+100;
kappa_incremental=-1*(1/b-2*lambda2)/(eta*lambda_L);%(2*lambda2)/(sigma*lambda_L)%
% kappa_incremental=0.2
conv_rate=-2*kappa_incremental*lambda2+1/2*sigma*kappa_incremental^2*lambda_L+kappa_incremental/b+1/a;
g=-conv_rate;
bb=a*lambda_gamma+sigma/2+(b*kappa_incremental*sigma^2*lambda_L)/4;
h=bb;
aa=-conv_rate;
% delta_p_infty=v*Ts/omega_low;
% delta_K4=1/(Ts*v)*sqrt((n-1)/n)*(q-1)*sqrt(bb/aa)*delta_p_infty+q
delta_K_LSR=sqrt((n-1)/n)*(q-1)/omega_low*sqrt(h*(2*n-1)/g)+q

K_compare=zeros(100,2);
k=1;

%%%%%--------------------estimate error bound--------------%%%%%%%%%%%%%%

for scale=kappa_init:0.1:kappa_end
    
    kappa=scale*lambda2/(eta*lambda_L);% use the proposed method to determin kappa
    kappa1=kappa/Ts;
    % find upper bound for method 1
    beta=rou_up^4/rou_low^2-rou_low^2;
    t_low=(d*rou_low*la)/(n*v);
    gamma=exp((4*(n-1)*kappa1*(d^2)*(rou_up^2)*la)/v);

    f2_alternate=0;
    for i=1:q-1
        f1_alternate=0;
        for j=0:q-1-i
            f1_alternate=f1_alternate+alpha^(j)*exp(-lambda2*kappa1*(j+1)*t_low);
        end
        f2_alternate=f2_alternate+f1_alternate;
    end
    CT_error=f2_alternate;
    
    syms j
    f=j*alpha^(q-1-j)*exp(-lambda2*kappa1*(q-j)*t_low);
    double(subs(symsum(f, j, 1, q-1)));
    truncation_error=(n+2)*kappa1*d^2*rou_up^2*la/(4*v)*((1-gamma^q)/(1-gamma)-q)*Ts;
    Deata_T=d*la/v*sqrt(beta/n)*symsum(f, j, 1, q-1)+truncation_error;
    delta_T=double(subs(Deata_T));
    delta_K05=(CT_error+truncation_error)/Ts;
    delta_K=double(subs(Deata_T/Ts))
    K_compare(k,1)=delta_K;
    
    %find upper bound for proposed method
    beta=lx^2/n*(omega_up^4/omega_low^2-omega_low^2);
    zeta=1-2*kappa*lambda2+0.5*eta*kappa^2*lambda_L;
    f2=0;
    for j=1:q-1
        f1=0;
        for i=0:q-1-j
            f1=f1+alpha^i*zeta^((i+1)/2*(lx*omega_low/(n*Ts*v)));
        end
        f2=f2+f1;
    end
    delta_K2=sqrt(beta)/(Ts*v)*sqrt((n-1)/n)*f2+q
%     syms i j;
%     f1=symsum(alpha^i*zeta^((i+1)/2*(lx*omega_low/(n*Ts*v))), i, 0, q-1-j);
%     f2=sqrt(beta)/(Ts*v)*sqrt((n-1)/n)*symsum(f1,j,1,q-1)+q;
%     delta_K21=double((f2))
    K_compare(k,2)=delta_K2;
    k=k+1;
end

s=svds(L);
syms i
fzeta=zeta^(i/2);
k_max=lx*rou_up*d/(Ts*v);
k_min=lx*rou_low*d/(n*Ts*v);
%fprintf('max displacement possible per stripe is');
max_disp=double(subs(kappa*s(1)*sqrt((n-1)/n)*sqrt((n-1)/n)*...
        (lx*rou_up*d)*symsum(fzeta,i,0,ceil(k_max)-1)));%max displacement possible per stripe
fprintf('max time difference based on max initial disorder is')
delta_K3=...%(1/(Ts*v)*sqrt((n-1)/n)*sqrt((n-1)/n)*(lx*rou_up*d)+1)+...
           (q-1)*(1/(Ts*v)*sqrt((n-1)/n)*sqrt((n-1)/n)*(lx*rou_up*d)*zeta^(k_min/2))+q

syms x y
rou=@(x,y) (rou_up-rou_low)/2*sin(x+y)+(rou_up+rou_low)/2;  %adjust workload density function here
%rou=@(x,y) (rou_up-rou_low)/2*sin(x.^2/5+y.^2/5)+(rou_up+rou_low)/2;

%mark=zeros(n+1,q);
%mark2=zeros(n+1,q);
mark2_move=zeros(n+1,q+1);
mark(1,1,1)=0;
mark(n+1,1,1)=lx;
k_max1(1)=1;
%m=zeros(n,q,100);
k_total=0;
for i=1:n-1
    mark(i+1,1,1)=i/n*lx;
end
for j=1:5000    % this is to get same initial workload
    for i=1:n
        initial_workload(i)=integral2(rou, mark(i,1,1), mark(i+1,1,1), 0, (d));
    end
    for i=1:n-1  %update partition mark
        mark(i+1,1,1)=mark(i+1,1,1)+kappa*(initial_workload(i+1)-initial_workload(i));
    end
end
mark2=mark;
for i=1:n
    m(i,1,1)=integral2(rou, mark(i,1,1), mark(i+1,1,1), 0, d);
    m2(i,1,1)=integral2(rou, mark(i,1,1), mark(i+1,1,1), 0, d);
end

%%%%%----------------------DTSC method------------------------%%%%%%%%%%%%%%

if (run_normal_algo)
    for j=1:q
        k=1;
        mark(:,j+1,k)=mark(:,j,k_max1(j));
        pos1(1:n,j,k)=zeros(n,1,1);
        temp=zeros(n);
        while k<=ceil(max(m(:,j,k_max1(j)))/(Ts*v))
            mark(1,j+1,k)=0;
            mark(n+1,j+1,k)=lx;
            if (normal_algo_pos_plot)
                if k>1 temp=pos1(:,j,k-1); end
                if (rem(j,2)==1) %odd stripe
                    parfor i=1:n
                        if (k==1)||(temp(i)<mark(i+1,j,k_max1(j)))
                            x2=sym('x2');
    %                         f1=1/2*((rou_up+rou_low)*(mark(i,j,k_max1(j))-x2)*(-d)+(rou_up-rou_low)*...
    %                                 (-sin(mark(i,j,k_max1(j))+(j-1)*d)+sin(x2+(j-1)*d)+sin(mark(i,j,k_max1(j))+j*d)-sin(x2+j*d)))==k*Ts*v;
                            f1=int(int(rou,x,mark(i,j,k_max1(j)),x2),y,(j-1)*d,j*d)==k*Ts*v;

                            pos1(i,j,k)=vpasolve(f1,x2,temp(i)+Ts*v*2/((rou_up+rou_low)*d));
                            if pos1(i,j,k)>=mark(i+1,j,k_max1(j))
                                pos1(i,j,k)=mark(i+1,j,k_max1(j));
                            end
                        else
                            pos1(i,j,k)=mark(i+1,j,k_max1(j));
                        end
                    end
                else
                    parfor i=1:n
                        if (k==1)||(temp(i)>mark(i,j,k_max1(j)))
                            x2=sym('x2');
                            f1=int(int(rou,x,x2,mark(i+1,j,k_max1(j))),y,(j-1)*d,j*d)==k*Ts*v;

                            pos1(i,j,k)=vpasolve(f1,x2,temp(i)-Ts*v*2/((rou_up+rou_low)*d));
                            if pos1(i,j,k)<=mark(i,j,k_max1(j))
                                pos1(i,j,k)=mark(i,j,k_max1(j));
                            end
                        else
                            pos1(i,j,k)=mark(i,j,k_max1(j));
                        end
                    end                    
                end
            end
            for i=1:n
                m(i,j+1,k)=integral2(rou, mark(i,j+1,k), mark(i+1,j+1,k), (j)*d, (j+1)*d);
            end
            for i=1:n-1  %update partition mark
                mark(i+1,j+1,k+1)=mark(i+1,j+1,k)+kappa*(m(i+1,j+1,k)-m(i,j+1,k));
                %mark(i+1,j+1)=mark(i+1,j+1)+kappa*(sum(m(1:n,j+1,k))-n*m(i,j+1,k));
            end
            k=k+1;
        end
        k_simul_opt_compare(j,1)=k-1;
        k_max1(j+1)=k-1;
        k_total=k_total+k_max1(j+1);
    end
    fprintf('actual number of steps in simulation is')
    k_total
end

%find optimum number of steps
k_opt=0;
% for i=1:q
%     k_opt=k_opt+ceil(integral2(rou, 0, lx, (i-1)*d, (i)*d)/(n*Ts*v));
%     k_simul_opt_compare(i,2)=ceil(integral2(rou, 0, lx, (i-1)*d, (i)*d)/(n*Ts*v));
% end
k_opt = ceil(integral2(rou, 0, lx, 0, lb)/(n*Ts*v));
fprintf('optimal possible number of steps is')
k_opt

%%%%%--------------------proposed method--------------%%%%%%%%%%%%%%
% a practical way of updating the partition mark miu
if (run_incremental_algo)
    k_total2=0;
    k_max2(1)=1;
    %pos2=zeros(n);
    pos_prev=zeros(n);
    for j=1:q
    %for j=1:2 %adjust here
        k=1;
        if (rem(j,2)==1)
            mark2(1,j+1,k)=0;
            mark2(2:n+1,j+1,k)=mark2(1:n,j,k_max2(j));
            %mark2(n+1,j+1)=lx;
            %pos2(1:n,j,k)=zeros(n,1,1);
            pos_prev=mark2(1:n,j,k_max2(j));
            pos2(1:n,j,k)=pos_prev;
            while k<=ceil(max(m2(:,j,k_max2(j)))/(Ts*v))   %adjust here
            %while k<=ceil(max(m2(:,j,k_max2(j)))/(Ts*v))+1500   

                % compute current workload for each agent
                extra_load_for_neighbour=zeros(n);
                for i=1:n
                    extra_load_for_neighbour(i+1)=integral2(rou, mark2(i+1,j+1,k), pos2(i,j,k), (j)*d, (j+1)*d);

                    m2(i,j+1,k)=integral2(rou, mark2(i,j,k_max2(j)), mark2(i+1,j+1,k), (j)*d, (j+1)*d)+...
                                extra_load_for_neighbour(i);
                end

                %update law  
                for i=1:n-1          
                    mark2_move(i+1,j+1)=mark2_move(i+1,j+1)+kappa_incremental*(m2(i+1,j+1,k)-m2(i,j+1,k));            
                    mark2(i+1,j+1,k)=pos2(i,j,k)+mark2_move(i+1,j+1);
                    if mark2(i+1,j+1,k)<mark2(i,j,k_max2(j))
                        mark2(i+1,j+1,k)=mark2(i,j,k_max2(j));
                    end
                end

                %update the position change
                temp=mark2(2:n+1,j+1,k);
                for i=1:n
                    if pos_prev(i)<mark2(i+1,j,k_max2(j))
                        x2=sym('x2');
                        f1=int(int(rou,x,mark2(i,j,k_max2(j)),x2),y,(j-1)*d,j*d)==k*Ts*v;
                        pos2(i,j,k+1)=vpasolve(f1,x2,pos_prev(i)+Ts*v*2/((rou_up+rou_low)*d));
                        if pos2(i,j,k+1)>=mark2(i+1,j,k_max2(j))
                            pos2(i,j,k+1)=mark2(i+1,j,k_max2(j));
                        end
                    else
                        pos2(i,j,k+1)=mark2(i+1,j,k_max2(j));
                    end
                    temp(i)=temp(i)+pos2(i,j,k+1)-pos_prev(i);
                    pos_prev(i)=pos2(i,j,k+1);
                end 
                mark2(2:n+1,j+1,k+1)=temp(1:n);        
                mark2(1,j+1,k+1)=0;
                k=k+1;
            end
            % compute current workload for each agent, again
            extra_load_for_neighbour=zeros(n);
            for i=1:n
                extra_load_for_neighbour(i+1)=integral2(rou, mark2(i+1,j+1,k), pos2(i,j,k), (j)*d, (j+1)*d);
                m2(i,j+1,k)=integral2(rou, mark2(i,j,k_max2(j)), mark2(i+1,j+1,k), (j)*d, (j+1)*d)+...
                            extra_load_for_neighbour(i);
            end
        else
            mark2(n+1,j+1,k)=la;
            mark2(1:n,j+1,k)=mark2(2:n+1,j,k_max2(j));
            pos_prev=mark2(2:n+1,j,k_max2(j));
            pos2(1:n,j,k)=pos_prev;
            while k<=ceil(max(m2(:,j,k_max2(j)))/(Ts*v))   %adjust here
            %while k<=ceil(max(m2(:,j,k_max2(j)))/(Ts*v))+1500   
                % compute current workload for each agent
                extra_load_for_neighbour=zeros(n);
%                 for i=2:n+1
%                     extra_load_for_neighbour(i-1)=integral2(rou, pos2(i,j,k), mark2(i,j+1,k), (j)*d, (j+1)*d);
% 
%                     m2(i,j+1,k)=integral2(rou, mark2(i,j,k_max2(j)), mark2(i+1,j+1,k), (j)*d, (j+1)*d)+...
%                                 extra_load_for_neighbour(i);
%                 end
                for i=n:-1:1
                    if i~=1
                        extra_load_for_neighbour(i-1)=integral2(rou, pos2(i,j,k), mark2(i,j+1,k), (j)*d, (j+1)*d);
                    end
                    m2(i,j+1,k)=integral2(rou, mark2(i,j+1,k), mark2(i+1,j,k_max2(j)), (j)*d, (j+1)*d)+...
                                extra_load_for_neighbour(i);
                end 
                %update law  
                for i=2:n          
                    mark2_move(i,j+1)=mark2_move(i,j+1)-kappa_incremental*(m2(i-1,j+1,k)-m2(i,j+1,k));            
                    mark2(i,j+1,k)=pos2(i,j,k)+mark2_move(i,j+1);
                    if mark2(i,j+1,k)>mark2(i+1,j,k_max2(j))
                        mark2(i,j+1,k)=mark2(i+1,j,k_max2(j));
                    end
                end

                %update the position change
                temp=mark2(1:n,j+1,k);
                for i=1:n
                    if pos_prev(i)>mark2(i,j,k_max2(j))
                        x2=sym('x2');
                        f1=int(int(rou,x,x2,mark2(i+1,j,k_max2(j))),y,(j-1)*d,j*d)==k*Ts*v;
                        pos2(i,j,k+1)=vpasolve(f1,x2,pos_prev(i)-Ts*v*2/((rou_up+rou_low)*d));
                        if pos2(i,j,k+1)<=mark2(i,j,k_max2(j))
                            pos2(i,j,k+1)=mark2(i,j,k_max2(j));
                        end
                    else
                        pos2(i,j,k+1)=mark2(i,j,k_max2(j));
                    end
                    temp(i)=temp(i)+pos2(i,j,k+1)-pos_prev(i);
                    pos_prev(i)=pos2(i,j,k+1);
                end 
                mark2(1:n,j+1,k+1)=temp(1:n);        
                 mark2(n+1,j+1,k+1)=la;

                k=k+1;
            end
            % compute current workload for each agent, again
            extra_load_for_neighbour=zeros(n);
            for i=n:-1:1
                if i~=1
                    extra_load_for_neighbour(i-1)=integral2(rou, pos2(i,j,k), mark2(i,j+1,k), (j)*d, (j+1)*d);
                end
                m2(i,j+1,k)=integral2(rou, mark2(i,j+1,k), mark2(i+1,j,k_max2(j)), (j)*d, (j+1)*d)+...
                            extra_load_for_neighbour(i);
            end            
        end
        
        k_simul_opt_compare(j,3)=k-1;
        k_max2(j+1)=k-1;
        k_total2=k_total2+k_max2(j+1)-1;
    end
    fprintf('Number of steps when using incremental adjustment is')
    k_total2
end

str=["$m^q_1(k)$";"$m^q_2(k)$";"$m^q_3(k)$";"$m^q_4(k)$";"$m^q_5(k)$";"$m^q_6(k)$";"$m^q_7(k)$";"$m^q_8(k)$"];
str2=["$\mu^q_1(k)$";"$\mu^q_2(k)$";"$\mu^q_3(k)$";"$\mu^q_4(k)$";"$\mu^q_5(k)$";"$\mu^q_6(k)$";"$\mu^q_7(k)$";"$\mu^q_8(k)$"];

%%%%%-------------------------------------plot--------------%%%%%%%%%%%%%%

if (plot_figure)
    for j=2:q
        figure
        xlabel('Sampling Time k','FontSize',15);
        ylabel('Workload','FontSize',15)
        hold on
        for i=1:n
            if (run_incremental_algo)
                plot(1:k_simul_opt_compare(j-1,3),squeeze(m2(i,j,1:k_simul_opt_compare(j-1,3))),'LineWidth',2)
                legend({str(i)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')

            end
            hold on
            if (run_normal_algo)
                plot(1:k_simul_opt_compare(j-1,1),squeeze(m(i,j,1:k_simul_opt_compare(j-1,1))),'LineWidth',2)
               
            end
        end
         %legend({str(1),str(2),str(3),str(4),str(5),str2(1),str2(2),str2(3),str2(4),str2(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
        legend({str2(1),str2(2),str2(3),str2(4),str2(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
    end
end
set(gcf, 'Position',  [100, 100, 550, 280])        

begin_lb = 0;
for i=1:n
    workload_boustro(i)=integral2(rou, 0, la, begin_lb, begin_lb+lb/n);
    begin_lb = begin_lb+lb/n;
    workload_boustro(i) = ceil(workload_boustro(i)/(Ts*v));
end

Delta_K_boustro = max(workload_boustro)-ceil(integral2(rou, 0, la, 0, lb)/(n*Ts*v))

%calculate for 5 robots specifically, each robot get a same area
workload_boustro(1)=integral2(rou, 0, la/2, 0, lb/5*2);
workload_boustro(2)=integral2(rou, la/2, la, 0, lb/5*2);

workload_boustro(3)=integral2(rou, 0, la/3, lb/5*2, lb);
workload_boustro(4)=integral2(rou, la/3, la/3*2, lb/5*2, lb);
workload_boustro(5)=integral2(rou, la/3*2, la, lb/5*2, lb);
for i=1:n
    workload_boustro(i) = ceil(workload_boustro(i)/(Ts*v));
end
Delta_K_boustro = max(workload_boustro)-ceil(integral2(rou, 0, la, 0, lb)/(n*Ts*v))

% Delta_K_formation = max(workload_formation)-ceil(integral2(rou, 0, la, 0, lb)/(n*Ts*v))

% tss=[10, 1, 0.1, 0.05, 0.01, 0.001];
% for iii=1:size(tss,2)
%     Ts = tss(iii);
%     pos_y_line = 0;
%     k_line = 0;
%     while pos_y_line<lb
%         pos_i = 2*[lb, lb, lb, lb, lb, lb, lb];
%         for i=1:n
%             y2=sym('y2'); 
%             if i==1
%                 f1=int(int(rou,x,(i-1)*la/n,i*la/n),y,pos_y_line,y2)==Ts*v;
%                 pos_i(i) = vpasolve(f1,y2, pos_y_line+Ts*v/(((rou_up+rou_low)/2)*la/n));
%             elseif int(int(rou,x,(i-1)*la/n,i*la/n),y,pos_y_line,min(pos_i)) > Ts*v
%                 f1=int(int(rou,x,(i-1)*la/n,i*la/n),y,pos_y_line,y2)==Ts*v;
%                 pos_i(i) = vpasolve(f1,y2, pos_y_line+Ts*v/(((rou_up+rou_low)/2)*la/n));
%             else
%                 pos_i(i) = min(pos_i);
%             end                
%         end
%         k_line = k_line+1;
%         pos_y_line = min(pos_i);
%     end
%     fprintf('for Ts = %f, Delta_K_line is', Ts);
%     Delta_K_line = k_line-ceil(integral2(rou, 0, la, 0, lb)/(n*Ts*v))
% end
