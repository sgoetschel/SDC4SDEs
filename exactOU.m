%Stratonovich solution to problem case OU -- dummy, return 0, it is not needed
function [exact] = exactOU(~, ~, t_begin, t_end, step_size, initial, bM)
    time = t_begin:step_size:t_end;
    exact = ((initial)+time.*bM).*0;   
end
