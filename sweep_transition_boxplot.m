load('result_transition_420.mat')
final_result_4=final_result;
clearvars -except final_result_4
load('result_transition_520.mat')
final_result_5=final_result;
clearvars -except final_result_4 final_result_5
load('result_transition_620.mat')
final_result_6=final_result;
% data = {final_result_4, final_result_5, final_result_6}; 
% boxplotGroup(data, 'PrimaryLabels', {'a' 'b' 'c'}, ...
%   'SecondaryLabels',{'Group1', 'Group2', 'Group3'}, 'InterGroupSpace', 2)


% prepare data
data=cell(4,3);
for ii=1:size(data,1)
Ac{ii}=final_result_6(:,ii);
Bc{ii}=final_result_5(:,ii);
Cc{ii}=final_result_4(:,ii);
end
data=vertcat(Cc,Bc,Ac);

xlab={'4','5','6'};
col=[...
    102,255,255, 200;
51,153,255, 200;
0, 0, 255, 200;
0,255,128, 200];
col=col/255;

multiple_boxplot(data,xlab,{'Proposed', 'Equal Workload','Cell Decomp', 'Line Formation'},col')
xlabel('Number of robots') 
ylabel('Completion Time') 
yticks(1600:300:3600);