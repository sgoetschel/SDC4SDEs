function[rhs, stochRhs, J, RhsIto, exact] = problem(~, beta)

rhs =@(x)zeros(1,length(x));
stochRhs ={@(x) beta*ones(1,length(x))};
J = {@(x) 0};
RhsIto =@(x) zeros(1, length(x));

%exact sol of the Stratonovich SDE
exact = @exactConst;
end