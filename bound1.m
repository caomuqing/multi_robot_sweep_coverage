clear all
close all

plot_figure=true;
run_incremental_algo=true;
run_normal_algo=true;
pos_plot=true;
%common for both methods
n=4;    %number of agents
la=4   %sweep region width
lb=6;   %sweep region length
d=1;    %agent effective width
q=ceil(lb/d);   %number of stripes
rou_low=1; %lowest workload density
rou_up=2   % max workload density
v=60 %sweeping rate (workload/time)
Ts=1e-90 % sampling period

kappa_init=0.5;
kappa_end=0.5;

Omega=[1,0,0;-1,1,0;0,-1,1;0,0,-1]*d*rou_low;
Omega_G=[1,0,0;-1,1,0;0,-1,1;0,0,-1]*d*rou_up;
L=[-1,1,0,0;0,-1,1,0;0,0,-1,1]*(-1); % communication graph, depends on
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

kappa=5;% use the proposed method to determin kappa
% find upper bound for method 1
beta=rou_up^4/rou_low^2-rou_low^2;
t_low=(d*rou_low*la)/(n*v);
gamma=exp((4*(n-1)*kappa*d^2*rou_up^2*la)/v);


syms j
f=j*alpha^(q-1-j)*exp(-lambda2*kappa*(q-j)*t_low);
double(subs(symsum(f, j, 1, q-1)));
CT_ERROR=double(subs(d*la/v*sqrt(beta/n)*symsum(f, j, 1, q-1)))
Deata_T=d*la/v*sqrt(beta/n)*symsum(f, j, 1, q-1)+...
        (n+2)*kappa*d^2*rou_up^2*la/(4*v)*((1-gamma^q)/(1-gamma)-q)*Ts;
delta_T=double(subs(Deata_T))
delta_K=double(subs(Deata_T/Ts))
K_compare(k,1)=delta_K;