%Stratonovich solution to problem case OU
function [sol] = exactOU(lambdaVec, sigma, t_begin, t_end, step_size, X0, eta)
    % lambdaVec = [lambda, mu]
    lambda = lambdaVec(1);
    mu = lambdaVec(2);
    
    nSteps = (t_end-t_begin)/step_size;
    sol = zeros(nSteps+1,1);
    sol(1) = X0;
    %time = t_begin:step_size:t_end;
    for kk = 2:nSteps+1
        dW = eta(kk) - eta(kk-1);             
        e  = exp(-lambda*step_size);
        m  = mu + (sol(kk-1) - mu) * e;         
        v  = (sigma^2 / (2*lambda)) * (1 - exp(-2*lambda*step_size));
        sol(kk) = m + sqrt(max(v,0)) * dW/sqrt(step_size);
    end
end
