close all hidden

str=["$m^q_1(k)$";"$m^q_2(k)$";"$m^q_3(k)$";"$m^q_4(k)$";"$m^q_5(k)$";"$m^q_6(k)$";"$m^q_7(k)$";"$m^q_8(k)$"];
str2=["$\mu^q_1(k)$";"$\mu^q_2(k)$";"$\mu^q_3(k)$";"$\mu^q_4(k)$";"$\mu^q_5(k)$";"$\mu^q_6(k)$";"$\mu^q_7(k)$";"$\mu^q_8(k)$"];
str3=["$\psi^q_1(k)$";"$\psi^q_2(k)$";"$\psi^q_3(k)$";"$\psi^q_4(k)$";"$\psi^q_5(k)$";"$\psi^q_6(k)$";"$\psi^q_7(k)$";"$\psi^q_8(k)$"];

for j=2:p.num_vert
    figure
    xlabel('Sampling Time k');
    ylabel('Operation cycle')
    hold on
    for i=1:n
            plot(1:kmax(j),squeeze(m(i,j,1:kmax(j))),'LineWidth',2)
            legend({str(i)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
        hold on
    end
     %legend({str(1),str(2),str(3),str(4),str(5),str2(1),str2(2),str2(3),str2(4),str2(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
    legend({str3(1),str3(2),str3(3),str3(4),str3(5)},'Interpreter','latex','FontSize',14, 'Orientation','vertical')
end