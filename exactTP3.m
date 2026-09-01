%% Stratonovich solution to problem case d='TP3'
% this represents an approximation of the exact solution
function [sol] = exactTP3(~, ~, t_begin, t_end, step_size, X0, eta)
%Stratonovich-SDE
    % time = t_begin:step_size:t_end;
    % int = 0;
    % for i=1:length(time)-1
    %     sol1(i) = X0*exp(0.5*time(i)+eta(1,i))./(1+X0*int);
    %     int = int + (time(i+1)-time(i)).*exp(0.5.*time(i+1) + eta(1,i+1));
    % end
    % sol1(end+1) =  X0*exp(0.5*time(end)+eta(1,end))./(1+X0*int);

    % fixed analytic sol approx
    nsteps = (t_end-t_begin)/step_size;
    t  = linspace(t_begin,t_end,nsteps+1);
    %W = zeros(1,nsteps+1);
    W = eta;
    integrand = exp(W+t/2);
    % Cumulative trapezoidal integration
    I = cumtrapz(t,integrand);

    % Analytical solution evaluated numerically
    sol = X0 * exp(W + t/2) ./ (1 + X0*I);

    % figure;
    % plot(t,sol,'LineWidth',1.5);
    % hold on; 
    % plot(t,sol1,'LineWidth',1.5);
    % xlabel('t');
    % ylabel('X(t)');
    % title('Analytical solution of the stochastic logistic equation');
    % legend('new', 'old')
    % grid on;
end