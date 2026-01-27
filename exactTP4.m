%% Stratonovich solution to problem case d='TP4', 2-dim
function [sol] = exactTP4(A, B, t_begin, t_end, step_size, initial, eta)
%Stratonovich-SDE
time = t_begin:step_size:t_end;
for k=1:length(time)
    sol(:,k) = expm((A-0.5*(B{1}*B{1}+B{2}*B{2})).*time(k) + B{1}.*eta(1,k) + B{2}.*eta(2,k))*initial;
end
end