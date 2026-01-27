%% Stratonovich solution to problem case d='TP3'
% this represents an approximation of the exact solution
function [sol] = exactTP3(~, ~, t_begin, t_end, step_size, ~, eta)
%Stratonovich-SDE
    time = t_begin:step_size:t_end;
    int = 0;
    for i=1:length(time)-1
    sol(i) = exp(0.5*time(i)+eta(1,i))./(2+int);
    int = int + (time(i+1)-time(i)).*exp(0.5.*time(i) + eta(1,i));
    end
    sol(end+1) =  exp(0.5*time(end)+eta(1,end))./(2+int);
end