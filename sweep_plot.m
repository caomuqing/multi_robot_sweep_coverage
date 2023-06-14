close all
if (true)
    for j=2:q
        figure
        xlabel('Sampling Time k','FontSize',14);
        ylabel('Workload','FontSize',14)
        hold on
        for i=1:n
            if (true) %to plot for proposed algo
                plot(1:k_simul_opt_compare(j-1,3),squeeze(m2(i,j,1:k_simul_opt_compare(j-1,3))),'LineWidth',2)
                legend({str(i)},'Interpreter','latex','FontSize',14, 'Orientation','vertical','Box','off')

            end
            hold on
            if (false) %to plot for DTSC algo
                plot(1:k_simul_opt_compare(j-1,1),squeeze(m(i,j,1:k_simul_opt_compare(j-1,1))),'LineWidth',2)
               
            end
        end
         %legend({str(1),str(2),str(3),str(4),str(5),str2(1),str2(2),str2(3),str2(4),str2(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical','NumColumns',2)
        legend({str2(1),str2(2),str2(3),str2(4),str2(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical','NumColumns',2, 'Box', 'On')
    end
end
        
set(gcf, 'Position',  [100, 100, 550, 290])