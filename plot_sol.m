function [] = plot_sol(time, t_begin, t_end, intervals, sol, solExactRef)

fig =figure;
hold on;
for i=1:size(sol,1)
    txt = ['sol = ',num2str(i)];
    plot(time, sol(i,:), '*b', 'DisplayName', txt), hold on
    
    txt = ['exact = ',num2str(i)];
    plot(time, solExactRef(i,:), '-og', 'DisplayName', txt), hold on
    
end

h=legend('Location', 'best');
set(h, 'FontSize', 12)

vertLine = linspace(t_begin, t_end, intervals+1);
for j=1:intervals+1
    line([vertLine(j) vertLine(j)],get(gca, 'ylim'), 'LineStyle', ':', 'Color', [0 0 0], 'HandleVisibility','off')
end
xlabel('time', 'FontSize', 16);
ylabel('u(t)','FontSize', 16);

hold on
set(gca,'FontSize',14)
hold on
end