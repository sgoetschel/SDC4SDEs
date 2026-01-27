%Stratonovich solution to problem case d='exp'
function [exact] = exactConst(~, ~, t_begin, t_end, step_size, initial, bM)
%Stratonovich-SDE
    time = t_begin:step_size:t_end;
    exact = (initial)+time.*bM;   
end