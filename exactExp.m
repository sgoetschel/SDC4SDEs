%Stratonovich solution to problem case d='exp'
function [exact] = exactExp(lambda, beta, t_begin, t_end, step_size, initial, eta)
%Stratonovich-SDE
    time = t_begin:step_size:t_end;
    exact = (initial).*exp( (lambda-0.5*beta^2).*time+beta.*eta);
end