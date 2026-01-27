%% Stratonovich solution to problem case d='TP2'
function [exact] = exactTP2(lambda, beta, t_begin, t_end, step_size, initial, eta)
%Stratonovich-SDE
    time = t_begin:step_size:t_end;
    exact = ( (1+initial).*exp(-2*lambda.*time+2*beta.*eta) -1 )./( (1+initial).*exp(-2*lambda.*time+2*beta.*eta) +1 );
end