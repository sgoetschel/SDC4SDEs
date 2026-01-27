function [] = plot_err(step_size, errWeak, errStrong, errL2, errT)

%creates error plot

figErrStrong = figure(2);
loglog(step_size, errStrong, 'ob')
xlabel('time step', 'FontSize', 16);
ylabel('strong convergence error','FontSize', 16);
hold on
set(gca,'FontSize',14)
set(gca, 'YScale', 'log')

figErrWeak = figure(3);
loglog(step_size, errWeak, 'ob')
xlabel('time step', 'FontSize', 16);
ylabel('weak convergence error','FontSize', 16);
hold on
set(gca,'FontSize',14)
set(gca, 'YScale', 'log')

figErrL2 = figure(4);
loglog(step_size, errL2, 'ob')
xlabel('time step', 'FontSize', 16);
ylabel('expected L2-error','FontSize', 16);
hold on
set(gca,'FontSize',14)
set(gca, 'YScale', 'log')

figErrT = figure(5);
loglog(step_size, errT, 'ob')
xlabel('time step', 'FontSize', 16);
ylabel('final time error','FontSize', 16);
hold on
set(gca,'FontSize',14)
set(gca, 'YScale', 'log')
end