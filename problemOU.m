function[rhs, stochRhs, J, RhsIto, exact] = problemOU(lambdaVec, beta)

% lambdaVec = [lambda, mu]
lambda = lambdaVec(1);
mu = lambdaVec(2);

rhs =@(x)lambda*(mu-x);
stochRhs ={@(x) beta*ones(1,length(x))};
J = {@(x) 0};
RhsIto =@(x) lambda*(mu-x);

%exact sol of the Stratonovich SDE
exact = @exactOU;
end
