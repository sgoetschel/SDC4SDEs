function[rhs, stochRhs, J, RhsIto, exact] = problem2(lambda, beta)

rhs =@(x) -lambda.*(1-x.^2);
stochRhs ={@(x) beta.*(1-x.^2)};
J = {@(x) -2*beta.*x};
RhsIto =@(x) -(lambda + beta^2.*x)*(1-x.^2);

exact =@exactTP2;
end