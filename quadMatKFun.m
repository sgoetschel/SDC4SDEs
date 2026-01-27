%% compute quadrature matrix for each expansion term
% based on Lagrange polynomial coeffs, time points (t_nodes), step size (deltaT), time grid (nodes),# expansion terms (m)
%by Lisa Fischer
function [quadMatK_c] = quadMatKFun(t_nodes, deltaT, nodes, m)
len = length(t_nodes);
quadMatK_c = zeros(len-1, len, m);
%quadMatK_s = zeros(len-1, len, m);

%compute quadrature matrix for each expansion term
for k=1:m
    %compute each matrix entry
    
    %each column consists of different lagrange polynomials
    %while keeping the integration boundaries fixed
    for col=1:len
        c = coeff(t_nodes, nodes, col);
        c = c/c(end);
        c(end) = [];
        
        switch mod(len, 2)
            case 0
                
                for i=1:len
                    c(i) = (-1)^(i) *c(i);
                end
                
            otherwise
                for i=1:len
                    c(i) = (-1)^(i+1)*c(i);
                end
        end
        
        c(1:end) = fliplr(c(1:end));
        %c(1)... coefficient of highest polynomial degree
        %c(end)... coefficient of lowest polynomial degree
        
        %change integral boundaries using row
        %each row has different time integrals
        for row=1:len-1
            quadMatK_c(row,col,k) = lagrange_cos_integrate( t_nodes, row+1, deltaT, c, k);
            %quadMatK_s(row,col,k) = lagrange_sin_integrate( t_nodes, row+1, deltaT, c, k);
        end
    end
end
end